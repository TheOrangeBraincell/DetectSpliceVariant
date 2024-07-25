# 10.06.24
# Mirjam Müller
# frequent_variants.py
#
# Reads in all genotype files, saves the variants with frequency >= 1% and then returns them as a table.
#
#
#########################################################################################################
#
# IMPORTS
#

import sys
import os

#########################################################################################################
#
# INPUT
#

input_folder = sys.argv[1]
files = os.listdir(input_folder)

out=open(sys.argv[2], "w")

out.write("# Variants with frequency >= 1% across all samples.\n")
counter=0
header_written=False
for file in files:
    gene=file.split("_")[0]
    counter+=1
    print(counter,  end="\r")
    with open(input_folder + file, "r") as genotypes:
        for line in genotypes:
            if line.startswith("Location"):
                #This is the header
                if header_written==False:
                    out.write("Gene\t"+line)
                    header_written=True
                samples=line.split("\t")[1:]
                continue

            #For everything else, we count the number of proper genotypes per line.
            var_id=line.split("\t")[0]
            entries=line.split("\t")[1:]

            count=0
            for e in entries:
                if e not in ("NE", "ND"):
                    count+=1
            
            #If they are more than 1 % of samples, we put it into the output file.
            if count>= len(entries)/100:
                out.write(gene+"\t"+line)
out.close()
