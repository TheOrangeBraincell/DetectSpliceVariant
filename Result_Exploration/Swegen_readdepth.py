# 17.07.24
# Mirjam Müller
# Swegen_readdepth.py
#
# Calculates average read depth from our output files for Swegen/SCANB matching variant locations
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

variant_file =sys.arv[2]

out=open(sys.argv[3], "w")

locations=[]
swegen_dict=dict()
with open(variant_file, "r") as variants:
    for line in variants:
        if line.startswith("chrom"):
            #header
            out.write(line.strip("\n")+"\tReadDepth\n")
            continue

        locations.append(line.split("\t")[0]+"_"+line.split("\t")[1])
        swegen_dict[line.split("\t")[0]+"_"+line.split("\t")[1]]=line.strip("\n")


for file in files:
    with open(file, "r") as read_depth:
        for line in file:
            if line.startswith("\home"):
                #header
                continue
            
            if line.split("\t")[0]+"_"+line.split("\t")[2] in locations:
                sum=0
                for i in line.split("\t")[4:]:
                    sum+=int(i)
                
                mean=sum/len(line.split("\t")[4:])

                out.write(swegen_dict[line.split("\t")[0]+"_"+line.split("\t")[2]]+"\t"+str(mean)+"\n")



            