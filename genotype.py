# -*- coding: utf-8 -*-
"""
Date: Thu Sep 28 10:39:07 2023
File Name: genotype.py
Author: Mirjam Karlsson-Müller

Description:
    Reads in the variant locations and determines genotype for samples
    that have a space holder ("-") instead of a genotype at a giving location,
    based on read depth found in Read_Depth directory.
    
List of Functions:
    
Procedure: 
    1.
    2.
    3.
    
Useage:
    python ../gitrepo/genotype.py -v ESR1_locations.tsv -r Read_Depth/ESR1_read_depth.tsv -g ESR1 -o ESR1_genotypes.tsv
    
Possible Bugs:
"""
#%% 0.0 Imports
import argparse
import time

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

#%% 1. Read in Read depth table.

read_depths=dict()

with open(args.read_depth, "r") as depths:
    for line in depths:
        if line.startswith("#Read Depth"):
            continue
        else:
            if line.split("\t")[1] not in read_depths:
                read_depths[line.split("\t")[1]]={line.split("\t")[0]: line.split("\t")[2].strip("\n")}
            else: 
                read_depths[line.split("\t")[1]][line.split("\t")[0]] = line.split("\t")[2].strip("\n")
            
            
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

print("Run time: {:.2f} seconds.".format(time.time()-start_time))  