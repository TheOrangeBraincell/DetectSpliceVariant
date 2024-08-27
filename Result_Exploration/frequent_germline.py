# 17.07.24
# Mirjam Müller
# frequent_germline.py
#
# Compares frequent variants (>1%) to dbsnp variants and returns a table with the variants that are found in both.
#
#
#########################################################################################################
#
# IMPORTS
#

import sys
import os
import subprocess

#########################################################################################################
#
# INPUT
#

germline = sys.argv[1]
frequent = sys.argv[2]
output = sys.argv[3]


#read in germline variant, save ids in list
dbsnp=[]
with open(germline, "r") as infile:
    for line in infile:
        dbsnp.append(line.strip("\n").split("\t")[-1])

#Read in frequent variants, if they are on the dbsnp list, write into output.

with open(frequent, "r") as infile, open(output, "w") as outfile:
    for line in infile:
        if line.startswith("#"):
            continue
        elif line.startswith("Gene"):
            #column names
            samples=line.strip("\n").split("\t")[2:]
            outfile.write("Gene\tLocation\t"+ "\t".join(samples)+"\n")
            continue

        variant_id=line.split("\t")[1]

        if variant_id in dbsnp:
            outfile.write(line)