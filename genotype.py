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
import pysam

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

def check_read(read):
    filter_pass=True
    
    if read.is_secondary:
        filter_pass=False
    
    if not read.is_mapped:
        filter_pass=False
    
    if read.is_duplicate:
        filter_pass=False
    
    if read.is_qcfail:
        filter_pass=False
    
    if read.query_name in reads_counted:
        filter_pass=False
        
    if filter_pass==True:
        reads_counted.append(read.query_name)
        
    #This also counts read that have an N or a D or I at this position
    #we only want to count reads that have an M at the variant position.
    
    #So either we implement that here. Or we skip that and just use 
    #subprocess + github
    
    
    return filter_pass

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
        bam_file=pysam.AlignmentFile(files[sample], 'rb', index_filename=files[sample][0:-1]+"i")
        
        count=bam_file.count(chrom, position, position+1, read_callback=check_read)
        #print(reads_counted)
        """
        #if we had a readcount threshold it would look smth like this.
        
        if count>=5:
            var_dict[key][sample]="0/0"
        else:
            var_dict[key][sample]="NE"
        """
        #For now we return read count for these cases.
        var_dict[key][sample]=str(count)
        
        

#%% 4. Write Genotype output file

for location in var_dict:
    new_line=location
    for sample in var_dict[location]:
        new_line+="\t"+var_dict[location][sample]
        
    out.write(new_line+"\n")



#%% Close output file
out.close()


# #%% Open location file and output genotype table

# locations=open(args.input, "r")
# genotypes=open(args.out, "w")

# genotypes.write("#Genotype Table\n")
# #%% Assign Genotypes


# #Iterate through lines in location table as well as genes simultaneously..
# for line in locations:
#     if line.startswith("#"):
#         #headers
#         genotypes.write(line)
#         #extract gene name
#         if line.startswith("#Variant Location"):
#             gene=line.strip("\n").split(" ")[7]
#         continue
#     if line.startswith("Location"):
#         #Thats the column names. We write those and extract sample names.
#         sample_names=line.strip("\n").split("\t")[1:]
#         #Write header for output file + column for gene name!
#         genotypes.write(line.split("\t")[0]+ "\tGene\t"+"\t".join(sample_names)+"\n")
#         continue
#     #Now all that is left is entries. shape: chr_position_ref_(alt)\t samples
#     #print(line.split("\t")[0].split("_")[0:2])
#     variant_chrom,variant_position = line.split("\t")[0].split("_")[0:2]
#     #infostring to use in first column of output.
#     infostring=line.split("\t")[0]
#     entries=line.strip("\n").split("\t")[1:]
    
#     #Since we checked for the vcf location, to only take variants within the genes range
#     #we do not need to check for this again.

#     for i in range(0, len(entries)):
#         if entries[i]=="-":
#             sample= sample_names[i]
#             #here goes the read_threshold
#             entries[i]="WIP"
#     #write it
#     genotypes.write(infostring+"\t"+gene+"\t"+"\t".join(entries)+"\n")


# print("Genotyping Done!                 \n", end="\r")


# #%% Close location file and output file.

# locations.close()
# genotypes.close()

# print("Files closed!")

#%% End time

print("Run time: {:.2f} seconds.".format(time.time()-start_time))   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
