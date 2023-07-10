# -*- coding: utf-8 -*-
"""
Date: Mon Jul 10 13:04:13 2023
File Name: gene_ranges.py
Author: Mirjam Karlsson-Müller

Description:
    Takes a list of gene names as input file, returns an output table containing
    location information on the genes in the input.
    
    
List of Functions:
    
Procedure: 
    1. Read Input file into list
    2. Extract information on genes from databases RefSeq and GENCODE
    3. Write output file.
    
Useage:
    
    
Possible Bugs:
"""

#%% 0.0 imports

import time
import argparse
import re

#%% 0.1 argparse

"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='Identify Alternative Splicing',
                                 usage='%(prog)s -o OUTPUT-FILE \
                                     -g GENCODE-FILE -r REFSEQ-FILE \
                                         [-c] "chrX:XXXXXX-XXXXXX" -as AS-TYPE',
                                 description="""Per AS event of interest, creates
                                 a table with the PSI scores supporting said event
                                 per sample in sample folder.""")

parser.add_argument('--out', '-o', required=True,
                    help="""Output file, bed format, containing info for each input gene.""")
parser.add_argument('--gencode', '-g', required=True,
                    help="""tsv file containing bed file information on 
                    annotated exons from GENCODE39 as well as gene names.""")
parser.add_argument('--refseq', '-r', required=True,
                    help="""tsv file containing bed file information on 
                    annotated exons from RefSeq as well as gene names.""")
parser.add_argument('--input', '-i', required=True,
                    help="Table containing gene names in first column.")

args = parser.parse_args()


#%% 0.2 Functions

#%% 0.3 Start Timer

start_time=time.time()

#%% 1. Read in input file into list.

gene_list=[]
#Read file
with open(args.input, "r") as infile:
    for line in infile:
        #skip potential header lines.
        if line.startswith("#"):
            continue
        gene_list.append(line.strip())

#%% 2. Read in database files (RefSeq and GENCODE) and save information.

gene_ranges=dict()

for file in [args.gencode, args.refseq]:
    with open(file, "r") as infile:
        for line in infile:
            # To exclude potential title lines/empty lines, formatting mistakes
            # Only takes chr[] and chr[]_random lines, in accordance with bam.
            if re.search(r"(.+)\t(?:([a-z]{3}[X,M,Y]?\d*)|([a-z]{3}[X,M,Y]?\d*)"
                         r".+_random)\t(\-?\+?)\t(\d+)\t(\d+)\t(\d+)\t([\d,]+)"
                         r"\t([\d,]+)\t(.+)", line):
                # specify groups.
                entry = re.search(r"(.+)\t([a-z]{3}[X,M,Y]?\d*).*\t"
                             r"(\-?\+?)\t(\d+)\t(\d+)\t(\d+)\t([\d,]+)"
                             r"\t([\d,]+)\t(.+)", line)
                #Assign variables to groups
                chrom=entry.group(2)
                strand=entry.group(3)
                start= entry.group(4)
                stop=entry.group(5)
                gene_name=entry.group(9)
                
                #If the gene name is not in the gene list, we dont need it.
                if gene_name not in gene_list:
                    continue
                
                if gene_name not in gene_ranges:
                    gene_ranges[gene_name]=[chrom, start, stop, gene_name, "0", strand]

                else:
                    #replace start if this transcripts start is smaller
                    if gene_ranges[gene_name][1]> start:
                        gene_ranges[gene_name][1]=start
                    #Replace stop if this transcripts start is bigger.
                    if gene_ranges[gene_name][2]<stop:
                        gene_ranges[gene_name][2]=stop
                    
                    """Note that it does not matter which strand we are on at this point,
                    as the start and stop just refer to smallest and biggest coordinate found
                    in the database file. Since we want to include everything, we are simply
                    looking for the biggest distance between smallest and biggest,
                    regardless if they are the start or stop coordinates."""


#%% 3. Write output file 

with open(args.out, "w") as outfile:
    #A bed file traditionally does not have a header, but i think it would help.
    outfile.write("chrom\tstart\tstop\tgene\tscore\tstrand\n")
    #write entries into file
    for gene in gene_list:
        outfile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(*gene_ranges[gene]))

#%% Time

print("Run time: {:.2f} seconds.".format(time.time()-start_time))  