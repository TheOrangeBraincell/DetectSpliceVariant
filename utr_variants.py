# -*- coding: utf-8 -*-
"""
Date: Tue Oct  3 14:20:53 2023
File Name:
Author: Mirjam Karlsson-Müller

Description:
    
    
List of Functions:
    
Procedure: 
    1.
    2.
    3.
    
Useage:
    
    
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

parser.add_argument('--input', '-i', required=True,
                    help='input file containing variant locations for a gene')
parser.add_argument('--output', '-o', required=True,
                    help='variant locations that are not in the 3 prime utr')
parser.add_argument('--gene', '-g', required=True,
                    help='name of the gene we are running for.')
parser.add_argument('--gff_annotation', '-gff', required=True,
                    help='gff annotation file, containing annotated exons.')


args = parser.parse_args()

#%% 0.2 Functions

#%% 0.3 Start Timer

start_time=time.time()

#%% 1. Read in gene annotation for gene.
annotation=dict()

with open(args.gff_annotation, "r") as gff:
    for line in gff:
        if line.startswith("#"):
            #header
            continue
        
        if args.gene in line:
            #right gene.
            annotation[line.split("\t")[2]+"_"+line.split("\t")[3]+"_"+ line.split("\t")[4]]=[line.split("\t")[2],line.split("\t")[3], line.split("\t")[4]]
            

#%% 2. Go through the variant locations, remove those that are ONLY in a 3' UTR

with open(args.input, "r") as locations, open(args.output, "w") as out:
    for line in locations:
        if line.startswith("#"):
            #descriptive header
            out.write(line)
        elif line.startswith("Location"):
            #column names
            out.write(line)
        else:
            #this is a variant entry.
            coordinate=int(line.split("\t")[0].split("_")[1])
            #check if its past the limit!
            
            flag="not_found"
            break_flag=False
            for region in annotation:
                if flag=="found":
                    print("stopped loop for ", coordinate)
                    break
                if coordinate >= int(annotation[region][1]) and coordinate<int(annotation[region][2]):
                    #found a region gene is in.
                    if annotation[region][0]=="three_prime_UTR" and flag=="not_found":
                        flag="3_UTR"
                    elif annotation[region][0]=="exon":
                        #We found the variant in a region that is not an UTR. Keep it.
                        out.write(line)
                        flag="found"
                        break
                    elif annotation[region][0]=="CDS":
                        #We found the variant in a region that is not an UTR. Keep it.
                        out.write(line)
                        flag="found"
                        break
                
            #at the end of the loop check flag:
            #if location still not found, keep it. But raise an error.
            #if location is UTR, dont keep it.
            if flag=="not_found":
                out.write(line)
            
#%% End Timer

print("Run time: {:.2f} seconds.".format(time.time()-start_time))  