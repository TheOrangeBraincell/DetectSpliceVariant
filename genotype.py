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
    for ESR1 f.e.
    python variants_in_AS_Pipeline/genotype.py -f Database/fpkm_table.tsv -i location_table_ESR1.tsv -o genotype_table_ESR1.tsv -r Database/gene_ranges.tsv -c "chr6:151656691-152129619"
    
    
Possible Bugs:
    location tables need to be created by vcf_location_table.py
    gene_ranges.tsv needs to be created with gene_ranges.py
"""

#%% Imports
import argparse
import time


#%% Start Timer 

start_time=time.time()

#%% Argparse

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
parser.add_argument('--ranges','-r', required=True, help= "File containing \
                    gene ranges created with refseq and gencode.")
                    

args = parser.parse_args()


#Server option: Make sure no outputs are accidentally deleted.
# if os.path.isfile(args.out):
#     print("The output file already exists, we assume it is complete.")
#     quit()
    
#%% Open location file and output genotype table

locations=open(args.input, "r")
genotypes=open(args.out, "w")

genotypes.write("#Genotype Table\n")
#%% Assign Genotypes


#Iterate through lines in location table as well as genes simultaneously..
for line in locations:
    if line.startswith("#"):
        #headers
        genotypes.write(line)
        #extract gene name
        if line.startswith("#Variant Location"):
            gene=line.strip("\n").split(" ")[7]
        continue
    if line.startswith("Location"):
        #Thats the column names. We write those and extract sample names.
        sample_names=line.strip("\n").split("\t")[1:]
        #Write header for output file + column for gene name!
        genotypes.write(line.split("\t")[0]+ "\tGene\t"+"\t".join(sample_names)+"\n")
        continue
    #Now all that is left is entries. shape: chr_position_ref_(alt)\t samples
    #print(line.split("\t")[0].split("_")[0:2])
    variant_chrom,variant_position = line.split("\t")[0].split("_")[0:2]
    #infostring to use in first column of output.
    infostring=line.split("\t")[0]
    entries=line.strip("\n").split("\t")[1:]
    
    #Since we checked for the vcf location, to only take variants within the genes range
    #we do not need to check for this again.

    for i in range(0, len(entries)):
        if entries[i]=="-":
            sample= sample_names[i]
            #here goes the read_threshold
            entries[i]="WIP"
    #write it
    genotypes.write(infostring+"\t"+gene+"\t"+"\t".join(entries)+"\n")


print("Genotyping Done!                 \n", end="\r")


#%% Close location file and output file.

locations.close()
genotypes.close()

print("Files closed!")

#%% End time

print("Run time: {:.2f} seconds.".format(time.time()-start_time))   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
