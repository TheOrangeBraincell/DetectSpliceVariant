# -*- coding: utf-8 -*-
"""
Date: Thu Sep 28 10:39:07 2023
File Name: genotype.py
Author: Mirjam Müller

Description:
    Reads in the variant locations and determines genotype for samples
    that have a space holder ("-") instead of a genotype at a giving location,
    based on read depth found in Read_Depth directory.
    
List of Functions:
    
Useage:
    python ../gitrepo/genotype.py -v ESR1_locations.tsv -r Read_Depth/ESR1_read_depth.tsv -g ESR1 -o ESR1_genotypes.tsv
    
Possible Bugs:
"""
#%% 0.0 Imports
import argparse
import time
import re
#%% 0.1 Argparse

"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='',
                                 usage='',
                                 description=""" """)

parser.add_argument('--variant', '-v', required=True,
                    help='input file containing variant locations.')
parser.add_argument('--read_depth', '-r', required=True,
                    help='input file containing read depths')
parser.add_argument('--output', '-o', required=True,
                    help='output file')
parser.add_argument('--gene', '-g', required=True,
                    help='output file')

args = parser.parse_args()

#%% 0.2 Functions

#%% 0.3 Start Timer

start_time=time.time()
print("Running Script genotype.py.")

#%% 1. Read in Read depth table.

read_depths=dict()

with open(args.read_depth, "r") as depths:
    for line in depths:
        if not line.startswith("chr"):
            #then thats the bam file names, i.e. way to get the sample.
            p=re.compile('(S0\d{5})')
            sample_names=p.findall(line)
        else:
            read_depths[line.split("\t")[3]]=dict()
            #add all samples to inner dictionary
            counts=line.strip("\n").split("\t")[6:]
            
            if len(sample_names)!=len(counts):
                print("There is an error with the number of columns!")
                quit()
            for i in range(0, len(sample_names)):
                read_depths[line.split("\t")[3]][sample_names[i]]=counts[i]

#%% 2. Read through variant location table, check read depth table, write output.


with open(args.variant, "r") as variants, open(args.output, "w") as out:
    for line in variants:
        if line.startswith("#Var"):
            out.write("#Genotype table for variants in gene "+ args.gene)
        elif line.startswith("#"):
            out.write(line)
        elif line.startswith("Location"):
            out.write(line)
            sample_names=list(line.strip("\n").split("\t")[1:])
        else:
            #Entries
            variant_ID=line.split("\t")[0]
            genotypes=line.strip("\n").split("\t")[1:]
            
            for index, entry in enumerate(genotypes):
                if entry=="-":
                    #genotype needs to be determined.
                    sample=sample_names[index]
                    for location in read_depths:
                        if variant_ID.startswith(location):
                            #Found match
                            count=int(read_depths[location][sample])
                            #if count > 23
                            if count > 23:
                                genotypes[index]= "0/0"
                            else:
                                genotypes[index]="NE"
            out.write(variant_ID+"\t"+"\t".join(genotypes)+"\n")

#%% End Timer

print("Run time genotype.py: {:.2f} seconds.".format(time.time()-start_time))  