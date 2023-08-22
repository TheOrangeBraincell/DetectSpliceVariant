# -*- coding: utf-8 -*-
"""
Date: Mon Aug 21 12:51:16 2023
File Name: Run_Pipeline.py
Author: Mirjam Karlsson-Müller

Description:
    Framework for the pipeline. Reads in inputs, allocates resources and runs 
    different parts of the pipeline.
    
List of Functions:
    Pipeline(gene): Runs entire pipeline for one gene.    
    
Useage:
    #Eventually run on Bianca with sbatch script.
    
    #But if run directly from console use:
    
python ../gitrepo/Run_Pipeline.py -s sample_file.txt -g genes.txt -o Outputs_Pipeline/ -gc ~/MasterProject/Database/hg38_GENCODE39_all.tsv -rs ~/MasterProject/Database/hg38_NCBI_all.tsv -vcf ~/Sample_Data -bam ~/Sample_Data -gff ~/MasterProject/Database/gencode.v39.annotation.gff3 --cores 6
    
Possible Bugs:
"""

#%% 0.0 Imports

import argparse
import math
import time
import multiprocessing as mp
import subprocess
import os
from shutil import which

#%% 0.1 Argparse

parser = argparse.ArgumentParser(prog='Score Alternative Splicing',
                                 usage='',
                                 description="""Creates
                                 a table with the PSI scores supporting said event
                                 per sample in sample folder.""")

parser.add_argument('--samples', '-s', required=True,
                    help="""file containing the names of all samples that the pipeline is to be run for.
                    One line per sample.""")
parser.add_argument('--genes', '-g', required=True,
                    help="File containing all HGNC gene names, the pipeline should run for. One gene name per line.")
parser.add_argument('--out', '-o', required=True,
                    help="""Output directory, containing pipeline outputs in subfolders.""")
parser.add_argument('--gencode', '-gc', required=True,
                    help="""tsv file containing bed file information on 
                    annotated exons from GENCODE39 as well as gene names. Instructions on
                    how to find the tables, can be found on github.""")
parser.add_argument('--refseq', '-rs', required=True,
                    help="""tsv file containing bed file information on 
                    annotated exons from RefSeq as well as gene names. Instructions on
                    how to find the tables, can be found on github.""")
parser.add_argument('--vcf', '-vcf', required=True,
                    help="""Folder containing vcf files for all samples.""")
parser.add_argument('--bam', '-bam', required=True,
                    help="""Folder containing bam files for all samples.""")
parser.add_argument('--gff', '-gff', required=True,
                    help="""gff3 annotation file from gencode (here used v39).""")
parser.add_argument('--cores', '-n', required=True,
                    help="""how many cores youd like the pipeline to use.""")



args = parser.parse_args()



#%% 0.2 Functions

def Pipeline(gene):
    #Initiate log file for the gene
    subprocess.call(["touch Log_Files_genes/"+gene+"_log.txt"], shell=True)
    
    #this one needs 2 cores, cause refseq and gencode go parallel.
    subprocess.call(["python Scripts/Identify_AS.py -o "+ args.out+"/AS_Events/AS_events_"+gene+".tsv -n "+ gene+ " -g "+ args.gencode+" -r " +args.refseq+ " >Log_Files_genes/"+gene+"_log.txt"], shell=True)
    
    #PSI, needs 3 cores to run 3 samples at a time. 
    subprocess.call(["python Scripts/PsiScores.py -i "+args.out+"/AS_Events/AS_events_"+gene+".tsv"+ " -o "+ args.out+"/PSI_Tables/"+gene+"_PSI.tsv"+ " -s bam_file_list.txt -is "+ parameter_string+ " >Log_Files_genes/"+gene+"_log.txt"], shell=True)
    
    #vcf parser (1 core)
    subprocess.call(["python Scripts/vcf_location_table.py -s vcf_file_list.txt -o "+ args.out+"/Variant_Locations/"+gene+"_locations.tsv -n "+ gene + " -r gene_ranges.bed >Log_Files_genes/"+gene+"_log.txt"], shell=True)

    #genotype (1 core)
    subprocess.call(["python Scripts/genotype.py -s bam_file_list.txt -i " +args.out+"/Variant_Locations/locations_"+gene+".tsv -o " + args.out+"/Genotype_Tables/"+gene+"genotypes.tsv >Log_Files_genes/"+gene+"_log.txt"], shell=True)
    #Write update into logfile.
    print(("Pipeline complete for gene "+ gene+"\n"))
    print((str(genes.index(gene)+1) + "/" +str(len(genes))+" genes.\n"))
    
    
#%% 0.3 Time Start

start_time=time.time()

#%% 0.4 Check if software requirements are given.

#samtools and pysam need to be installed.
#pysam. It would throw an error when its not installed. But we want a nicer error.

try:
    import pysam
except:
    print("The module pysam was not found. Install pysam v0.19.0 over conda with conda install -c bioconda pysam=0.19.0 ")
    quit()
#if it didnt quit, it managed to load pysam
    
#samtools
if which("samtools") is None:
    print("samtools has not been found. Install samtools v.1.15 over conda with conda install -c bioconda samtools=1.15")
    quit()

#%% 1. Preparations for pipeline

"Check if input files exists"
if not os.path.isfile(args.samples):
    print("The sample file does not exist. Please input a valid file containing the sample names.")
    quit()

if not os.path.isfile(args.genes):
    print("The gene file does not exist. Please input a file containing the names of the genes the pipeline should run for.")
    quit()

if not os.path.isfile(args.refseq):
    print("The RefSeq file does not exist.")
    quit()

if not os.path.isfile(args.gencode):
    print("The GENCODE file does not exist.")
    quit()
    

"Check for output directories"
if not os.path.isdir(args.out):
    #Then we need to create the directory and all subdirectories.
    os.mkdir(args.out)
#check for sub directories needed for results
if not os.path.isdir(args.out+"/AS_Events"):
    os.mkdir(args.out+"/AS_Events")
if not os.path.isdir(args.out+"/PSI_Tables"):
    os.mkdir(args.out+"/PSI_Tables")
if not os.path.isdir(args.out+"/Variant_Locations"):
    os.mkdir(args.out+"/Variant_Locations")
if not os.path.isdir(args.out+"/Genotype_Tables"):
    os.mkdir(args.out+"/Genotype_Tables")
#We also need a folder for the logfiles for the pipeline run for each gene.
if not os.path.isdir("Log_Files_genes"):
    os.mkdir("Log_Files_genes")
    
#Now the structure to save outputs is given.

"Create necessary files."
#BAM file list (we dont check if its already there, as it might lead to faulty outputs if the sample list changes)
#subprocess.call(["find "+ args.bam+ " -name alignment.bam >bam_file_list.txt"], shell=True)
#VCF list (same as for bam)
#subprocess.call(["find "+ args.vcf + " -name variants-annotated.vcf >vcf_file_list.txt"], shell=True)
#Gene ranges file (same as bam and vcf but for genes.)
#This requires the python file gene_ranges to run.
#subprocess.run(["python Scripts/gene_ranges.py -g "+ args.gencode+ " -r "+ args.refseq+ " -o gene_ranges.bed -i "+ args.genes], shell=True)
#Mean insert size and standard deviation
gff_file=str(args.gff)
print(gff_file)
parameter_string=subprocess.check_output(["/bin/bash" ,"Scripts/Insert_Size_Bamfiles.sh" , gff_file], shell=True)
print(parameter_string)
quit()
#%% 2. Read in sample and gene names.

sample_names=[]
with open(args.samples) as sample_file:
    for line in sample_file:
        sample_names.append(line.strip("\n"))

genes=[]
with open(args.genes) as gene_file:
    for line in gene_file:
        genes.append(line.strip("\n"))
        

#%% 3. Start Pipeline for every gene separately, allocate resources.

if __name__ == '__main__':
    with mp.Pool(math.floor(int(args.n)/3)) as p:
        p.map(Pipeline, genes)
        
#%% Time End

print("Run time: {:.2f} seconds.".format(time.time()-start_time))
