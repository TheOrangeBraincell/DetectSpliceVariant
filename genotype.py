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
import subprocess
import multiprocessing as mp

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
print("Starting Genotyping script!")
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

print("Reading in Location Table:...", end="\r")
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
            sample_names=[i for i in line.strip("\n").split("\t")[1:]]
            out.write(line.split("\t")[0]+ "\tGene\t"+"\t".join(sample_names)+"\n")
            
            continue
        
        #Else its a variant entry. Save into dictionary.
        location=line.split("\t")[0]
        var_dict[location]=dict()
        for i in range(0,len(sample_names)):
            var_dict[location][sample_names[i]]=line.strip("\n").split("\t")[i+1]

print("Reading in Location Table: Done! \n", end="\r")
#%% 3. Go through bam file list, call sequence depth for every variant location.

def Genotype(sample):
    #create output dictionary
    sample_dict=dict()
    sample_dict[sample]=dict()
    print("Finding Genotypes for sample ", sample_names.index(sample)+1, "/", len(sample_names), "...")
    bam_file=files[sample]
    #iterate through all locations and extract locations we need to get with samtools
    with open("temp_out_"+sample+".txt", "o") as out:
        for key in var_dict:
            if var_dict[key][sample]=="-":
                out.write(var_dict[key].split("_")[0]+"\t"+str(int(var_dict[key].split("_")[1])+1)+"\t"+str(int(var_dict[key].split("_")[1])+1)+"\n")
        
        
    
    subprocess.run("cat temp_out_"+sample+".txt | while read coordinate; do samtools depth -a -s -r $coordinate " + bam_file+ " > counts_"+sample+".txt", shell=True)
    
    #Check output
    #Put it into sample dict.
    #Then merge after multiprocessing according to psi
    #dont forget to delete temp out file and counts file 
                
    count=int(str(output).split("\\t")[2][0:-3])
    if count>=23:
        result[key][sample]="0/0"
    else:
        result[key][sample]="NE"

    return sample_dict
# def Genotype(sample):
#     print("Finding Genotypes for sample ", sample_names.index(sample)+1, "/", len(sample_names), "...")
#     bam_file=files[sample]
#     #iterate through all locations and extract seq depth.
#     for key in var_dict:
#         if var_dict[key][sample]!="-":
#             #We already have a genotype there, no need to check count.
#             continue
        
#         chrom=key.split("_")[0]
#         position=int(key.split("_")[1])
#         #samtools takes 1-based coordinates. And the end is inclusive. So to get one nucleotide, the start and stop are the same.
#         region_string=chrom+":"+str(position+1)+"-"+str(position+1)
#         output=check_output("samtools depth -s -r "+region_string+ " " + bam_file, shell=True)
#         #samtools depth returns nothing if the read count is 0.
#         if str(output)=="b''":
#             count=0
#         else:
#             count=int(str(output).split("\\t")[2][0:-3])
        
#         #if we had a readcount threshold it would look smth like this.
        
#         if count>=23:
#             result[key][sample]="0/0"
#         else:
#             result[key][sample]="NE"
        
#         """
#         #For now we return read count for these cases.
#         var_dict[key][sample]=str(count)
#         """
  

# result=mp.Manager().dict(var_dict)
 with mp.Pool(3) as pool:
     result=pool.map(Genotype, sample_names)

# print(result)

print("All genotypes have been assigned!")

#%% 4. Write Genotype output file

for location in var_dict:
    new_line=location
    for sample in var_dict[location]:
        new_line+="\t"+var_dict[location][sample]
        
    out.write(new_line+"\n")

print("Genotype script complete!")

#%% Close output file
out.close()

#%% End time

print("Run time: {:.2f} seconds.".format(time.time()-start_time))  