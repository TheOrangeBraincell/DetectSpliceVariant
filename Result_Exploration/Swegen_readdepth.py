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

variant_file =sys.argv[2]

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

print(len(locations))
counter=0

for file in files:
    with open(input_folder+file, "r") as read_depth:
        for line in read_depth:
            if "mirjam" in line:
                #header
                continue
            
            if line.split("\t")[0]+"_"+line.split("\t")[2] in locations:
                counter+=1
                print(counter, end="\r")
                sum=0
                for i in line.strip("\n").split("\t")[6:]:
                    sum+=int(i)
                
                mean=round(sum/len(line.split("\t")[4:]), ndigits=2)
                
                out.write(swegen_dict[line.split("\t")[0]+"_"+line.split("\t")[2]]+"\t"+str(mean)+"\n")


out.close()