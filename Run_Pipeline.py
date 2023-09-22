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
    
    #run directly from console use:
    
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
    #Run bash script with all subprocesses for the entire pipeline in it.
    subprocess.run(["./Scripts/Run_Pipeline.sh " +gene+ " " +args.out +" "+ args.refseq + " " + args.gencode +" " + parameter_string])
    """
    #Initiate log file for the gene
    header="#Log file for "+gene
    subprocess.run(["echo '"+ header+ "' > Log_Files_genes/"+gene+"_log.txt"], shell=True)
    
    #this one needs 2 cores, cause refseq and gencode go parallel.
    subprocess.run(["python Scripts/Identify_AS.py -o "+ args.out+"AS_Events/"+gene+"_AS_events.tsv -n "+ gene+ " -g "+ args.gencode+" -r " +args.refseq+ " >>Log_Files_genes/"+gene+"_log.txt"], shell=True)
    #Because AS needs only 2 cores, and vcf needs 1, we can run those simultaneously.
    #vcf parser (1 core)
    subprocess.run(["python Scripts/vcf_location_table.py -s vcf_file_list.txt -o "+ args.out+"Variant_Locations/"+gene+"_locations.tsv -n "+ gene + " -r gene_ranges.bed >>Log_Files_genes/"+gene+"_log.txt"], shell=True)
    
    #PSI, needs 3 cores to run 3 samples at a time. 
    subprocess.run(["python Scripts/PsiScores.py -i "+args.out+"AS_Events/"+gene+"_AS_events.tsv"+ " -o "+ args.out+"PSI_Tables/"+gene+"_PSI.tsv"+ " -s bam_file_list.txt -is '"+ parameter_string+ "' >>Log_Files_genes/"+gene+"_log.txt"], shell=True)
    
    #genotype (3 cores)
    subprocess.run(["./Scripts/Read_Depth.sh " + gene])
    #wait, then read read depth table off.
    subprocess.run(["python genotype.py -i " + gene+ "_read_depth.tsv -o " + args.out+"Genotype_Tables/"+gene+"_genotypes.tsv"])
    #subprocess.run(["python Scripts/genotype.py -s bam_file_list.txt -i " +args.out+"Variant_Locations/"+gene+"_locations.tsv -o " + args.out+"Genotype_Tables/"+gene+"_genotypes.tsv >>Log_Files_genes/"+gene+"_log.txt"], shell=True)
    """
    #Write update into logfile.
    print(("Pipeline complete for gene "+ gene+"\n"))
    print((str(genes.index(gene)+1) + "/" +str(len(genes))+" genes.\n"))
    
    
#%% 0.3 Time Start

start_time=time.time()
print("Starting Pipeline...")

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

"Check if Scripts folder is here and all scripts are present"
if not os.path.isdir("/Scripts"):
    print("The Scripts directory with the parts of the pipeline is missing.")
    quit()
else:
    #Check for every script needed.
    #gene_ranges.py
    if not os.path.isfile("/Scripts/gene_ranges.py"):
        print("The gene_ranges.py script is missing from the Scripts directory. Please make sure all scripts are present before starting the Pipeline.")
    #Identify_AS.py
    if not os.path.isfile("/Scripts/Identify_AS.py"):
        print("Identify_AS.py script is missing from the Scripts directory. Please make sure all scripts are present before starting the Pipeline.")
    
    #PsiScores.py
    if not os.path.isfile("/Scripts/PsiScores.py"):
        print("The PsiScores.py script is missing from the Scripts directory. Please make sure all scripts are present before starting the Pipeline.")
    
    #genotype.py
    if not os.path.isfile("/Scripts/genotype.py"):
        print("The genotype.py script is missing from the Scripts directory. Please make sure all scripts are present before starting the Pipeline.")
    
    #Run_Pipeline.sh
    if not os.path.isfile("/Scripts/Run_Pipeline.sh"):
        print("The Run_Pipeline.sh script is missing from the Scripts directory. Please make sure all scripts are present before starting the Pipeline.")
    
    #Read_Depth.sh
    if not os.path.isfile("/Scripts/Read_Depth.sh"):
        print("The Read_Depth.sh script is missing from the Scripts directory. Please make sure all scripts are present before starting the Pipeline.")
    
    #vcf_location_table.py
    if not os.path.isfile("/Scripts/vcf_location_table.py"):
        print("The vcf_location_table.py script is missing from the Scripts directory. Please make sure all scripts are present before starting the Pipeline.")
    

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
if not os.path.isdir(args.out+"/Read_Depths"):
    os.mkdir(args.out+"/Read_Depths")
#We also need a folder for the logfiles for the pipeline run for each gene.
if not os.path.isdir("Log_Files_genes"):
    os.mkdir("Log_Files_genes")
#And for Database for insert size
if not os.path.isdir("Database"):
    os.mkdir("Database")

print("Directories created if needed, input files checked.")
#Now the structure to save outputs is given

#%% 2. Read in sample and gene names.

sample_names=[]
with open(args.samples) as sample_file:
    for line in sample_file:
        sample_names.append(line.strip("\n"))

genes=[]
with open(args.genes) as gene_file:
    for line in gene_file:
        genes.append(line.strip("\n"))
        
print("Sample names are read in, genes are read in. ")

#%% 3. Creating needed files.

"Create necessary files."
#BAM file list (we dont check if its already there, as it might lead to faulty outputs if the sample list changes)
#find all 
#subprocess.run(["find "+ args.bam+ " -name alignment.bam >all_bam_file_list.txt"], shell=True)
#Exclude whats not on the list
out=open("bam_file_list.txt", "w")
with open("all_bam_file_list.txt", "r") as all_bam:
    for line in all_bam:
        sections=line.split("/")
        for i in sections:
            if i.startswith("S00"):
                #thats the sample name
                if i in sample_names:
                    out.write(line)
                    
out.close()
#VCF list (same as for vcf)
#subprocess.run(["find "+ args.vcf + " -name variants-annotated.vcf >all_vcf_file_list.txt"], shell=True)
out=open("vcf_file_list.txt", "w")
with open("all_vcf_file_list.txt", "r") as all_vcf:
    for line in all_vcf:
        sections=line.split("/")
        for i in sections:
            if i.startswith("S00"):
                #thats the sample name
                if i in sample_names:
                    out.write(line)
                    break
out.close()
#Gene ranges file (same as bam and vcf but for genes.)
#This requires the python file gene_ranges to run.
#subprocess.run(["python Scripts/gene_ranges.py -g "+ args.gencode+ " -r "+ args.refseq+ " -o gene_ranges.bed -i "+ args.genes], shell=True)
#Mean insert size and standard deviation, used to be a script. Now its here.
#subprocess.run("cat "+args.gff+""" | grep "three_prime_UTR" | awk 'BEGIN{OFS="\t"} ($5-$4>1000) {print($1 OFS $4 OFS $5)}' > Database/3_UTR.bed""", shell=True)
#subprocess.run("""cat bam_file_list.txt | while read i; do samtools view $i -L Database/3_UTR.bed -h -f PROPER_PAIR -F UNMAP,SECONDARY,QCFAIL --subsample 0.25 | head -n 500000 | samtools stats | grep -A 1 "insert size average"; done > average_insert.txt""", shell=True)
#subprocess.run("""awk 'BEGIN{FS="\t"} (NR %2 == 1) {s+=$3} (NR %2 == 0) {sd+=$3^2} END{print("Mean", s/(NR/2) , "Standard Deviation", sqrt(sd/(NR/2)))}' average_insert.txt > parameter_string.txt""", shell=True)

with open("parameter_string.txt", "r") as param_file:
    for line in param_file:
        if line.startswith("Mean"):
            parameter_string=line.strip("\n")
print("Additional Files needed for all genes created.")


#%% 3. Start Pipeline for every gene separately, allocate resources.

if __name__ == '__main__':
    with mp.Pool(math.floor(int(args.cores)/3)) as p:
        p.map(Pipeline, genes)
        
#%% Time End
print("Pipeline Run completed!")
print("Run time: {:.2f} seconds.".format(time.time()-start_time))
