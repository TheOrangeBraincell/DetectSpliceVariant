# -*- coding: utf-8 -*-
"""
Date: Thu Jan 19 10:24:12 2023
File Name: genotype.py
Author: Mirjam Karlsson-Müller

Description:
    Assigns genotypes NE (not expressed) or 0/0 (homozygous reference) to location table entries marked with "-"
    
List of Functions:
    none
    
Procedure: 
    1. Reads in gene range file created by gene_ranges.py
    2. Iterates through genes and location table at the same time (as they are sorted by coordinates thats ok)
    3. If variant in gene, assign genotype based on fpkm (>=10 0/0, otherwise NE), otherwise assign genotype NE (not expressed)

Input: 
    - Sorted Location table as created by vcf_location_table.py and sorted by Sort_Locations.sh
    - sample folder containing gene.tsv


Useage:
    python ../gitrepo/genotype.py -s bam_file_list.txt -i location_table_ESR1.tsv -o genotype_table_ESR1.tsv 
    
    
Possible Bugs:
    location tables need to be created by vcf_location_table.py
    gene_ranges.tsv needs to be created with gene_ranges.py
"""

#%% Imports
import argparse
import time
from subprocess import check_output

#%% 0.1 Argparse

parser = argparse.ArgumentParser(prog='assigning genotypes',
                                 usage='%(prog)s -s INPUT-FOLDER -o OUTPUT \
                                     -i LOCATION-TABLE [-c] "chrX:XXXXXX-XXXXXX"\
                                         -r RANGE-TSV -f FPKM-TABLE',
                                 description="""Creates a genotype table out of
                                 a location table. Containing genotypes for location x sample.""")

parser.add_argument('--input', '-i', required=True,
                    help="""Input file containing location table.""")
parser.add_argument('--out', '-o', required=True,
                    help="""Output file containing genotypes.""")
parser.add_argument('--samples', '-s', required=True,
                    help='file containing the paths to the bam files of all samples.')
                    

args = parser.parse_args()


#Server option: Make sure no outputs are accidentally deleted.
# if os.path.isfile(args.out):
#     print("The output file already exists, we assume it is complete.")
#     quit()

#%% 0.2 Functions

#%% 0.3 Start Timer 

start_time=time.time()

#%% 1. Read in list of bam files. 

print("Reading in BAM files...", end="\r")

bam_file_list=[]
with open(args.samples, "r") as infile:
    for line in infile:
        bam_file_list.append(line.strip())

#Assuming that each bam.file is in its own sample folder with the sample name as folder name. 
sample_names=sorted(list(set([i.split("/")[-4] for i in bam_file_list])))

#make that a dictionary.
files=dict()
for sample in sample_names:
    for file in bam_file_list:
        if sample in file:
            files[sample]=file

print("Reading in BAM files: Done! \n", end="\r")

#%% 2. Read in location table

out=open(args.out, "w")
var_dict=dict()

with open(args.input, "r") as infile:
    for line in infile:
        if line.startswith("#"):
            #headers
            out.write(line)
            #extract gene name
            if line.startswith("#Variant Location"):
                gene=line.strip("\n").split(" ")[7]
            continue
        if line.startswith("Location"):
            #Write header for output file + column for gene name!
            out.write(line.split("\t")[0]+ "\tGene\t"+"\t".join(sample_names)+"\n")
            continue
        
        #Else its a variant entry. Save into dictionary.
        location=line.split("\t")[0]
        var_dict[location]=dict()
        for i in range(0,len(sample_names)):
            var_dict[location][sample_names[i]]=line.strip("\n").split("\t")[i+1]


#%% 3. Go through bam file list, call sequence depth for every variant location.


for sample in sample_names:
    #iterate through all locations and extract seq depth.
    for key in var_dict:
        if var_dict[key][sample]!="-":
            #We already have a genotype there, no need to check count.
            continue
        
        reads_counted=[]
        chrom=key.split("_")[0]
        position=int(key.split("_")[1])
        bam_file=files[sample]
        #samtools takes 1-based coordinates. And the end is inclusive.
        region_string=chrom+":"+str(position+1)+"-"+str(position+1)
        output=check_output("samtools depth -s -r "+region_string+ " " + bam_file, shell=True)
        #samtools depth returns nothing if the read count is 0.
        if str(output)=="b''":
            count=0
        else:
            count=int(str(output).split("\\t")[2][0:-3])
        
        #if we had a readcount threshold it would look smth like this.
        
        if count>=23:
            var_dict[key][sample]="0/0"
        else:
            var_dict[key][sample]="NE"
        """
        #For now we return read count for these cases.
        var_dict[key][sample]=str(count)
        """
        

#%% 4. Write Genotype output file

for location in var_dict:
    new_line=location
    for sample in var_dict[location]:
        new_line+="\t"+var_dict[location][sample]
        
    out.write(new_line+"\n")



#%% Close output file
out.close()
