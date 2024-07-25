# 05.07.24
# Mirjam Müller
# germline_variants.py
#
# Reads in all genotype files, compares them to dbsnp and returns the variants it finds in dbsnp as well, assuming they are germline variants.
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

input_folder = sys.argv[1]
files = os.listdir(input_folder)

output_file = sys.argv[2]


out=open(sys.argv[2], "w")
out.write("Chromosome\tPosition\tdbSNP\tvar_ID\n")

for file in files:
    with open(input_folder+file, "r") as table:
        for line in table:
            #we dont need the header
            if line.startswith("Location"):
                continue
            
            var_id = line.split("\t")[0]
            chrom = var_id.split("_")[0].strip("chr")
            location = int(var_id.split("_")[1])+1

            #Check if in dbsnp
            region = chrom+":"+str(location)+"-"+str(location)
            query_result = subprocess.run("tabix 00-All.vcf.gz "+ region, shell=True, capture_output=True, text=True).stdout
            #grab output
            if query_result!="":
                out.write(chrom+"\t"+str(location)+"\t"+query_result.split("\t")[2]+"\t"+var_id+"\n")
            
