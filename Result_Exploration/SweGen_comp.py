# 12.06.24
# Mirjam Müller
#
# SweGen_comp.py
#
#
# This script reads in the variants found in SweGen, saves their locations. 
# Then it reads through the genotype files and calculates AF for the locations found in SweGen as well. 
# It then returns a table with my allele frequencies, SweGen allele frequencies and variant ID.
#
#
#
######################################################################################################
# IMPORT

import sys
import os

# Handle input

swegen = sys.argv[1]
genotype_folder = sys.argv[2]

genotype_files = os.listdir(genotype_folder)

# Read in Swegen information
swegen_dict=dict()
with open(swegen, "r") as swefile:
    for line in swefile:
        if line.startswith("#"):
            #header
            continue

        #Remove variants that did not get a pass on quality
        if line.split("\t")[6]!="PASS":
            continue

        # We only use single nucleotide variants. so if they are indels, skip.
        if len(line.split("\t")[3])>1: #more than one reference base
            continue
        # More than one alternative base, we also remove non nucleotide characters
        alt_bases=line.split("\t")[4].split(",")
        current=""
        for base in alt_bases:
            if base in ["A", "C", "G", "T"]:
                current+=base
  
        alt = "".join(sorted(current)) 
        
        var_id = "_".join([line.split("\t")[0], line.split("\t")[1], line.split("\t")[3]])

        #save entry in dictionary
        swegen_dict[var_id]=[alt, line.split("\t")[7].split(";")[1].split("=")[1]]

print("Read in Swegen file.")

#open output file
out = open("swegen_scanb345_variants.tsv", "w")
#Write column names
out.write("chrom\tposition\tref\talt_SW\talt_SC\tAF_SW\tAF_SC\tMajority_VarAllele\n")

# Now we read in the genotype files and save those variants that have an equivalent in SweGen.
scanb_dict=dict()
for file in genotype_files:
    with open("../../Results_170524/Genotypes_in_Exons/"+ file, "r") as genotypes:
        for line in genotypes:
            if line.startswith("Location"):
                #header
                continue
            
            #Adjust coordinates so they are 1 based like Swegens. (and like variants usually are annotated.)
            position=str(int(line.split("\t")[0].split("_")[1])+1)
            #match by chrom, position and reference base
            new_id=line.split("\t")[0].split("_")[0]+"_"+position+"_"+line.split("\t")[0].split("_")[2]

            if new_id in swegen_dict:
                # calculate allele frequency.
                homR=0
                het=0
                homA=0
                for gt in line.split("\t")[1:]:
                    if gt=="ND" or gt =="NE":
                        continue
                    if gt=="0/0":
                        homR +=1
                    elif gt in ["1/1", "2/2", "3/3"]:
                        homA+=1
                    elif gt in ["0/1", "1/2", "0/2", "1/3", "2/3", "0/3"]:
                        het+=1
                if homR + het+ homA <345:
                    stop=True
                    continue    
                AF = (het+2*homA)/(2*(het+homA+homR))
                if homA > het:
                    majority = "homA"
                elif het > homA:
                    majority = "het"
                else:
                    majority = "both"


                alt=line.split("\t")[0].split("_")[3].strip(")").strip("(").strip(",")
                # Write output line
                #print(new_id.split("_")[0]+"\t"+new_id.split("_")[1]+"\t"+new_id.split("_")[2]+"\t"+swegen_dict[new_id][0]+"\t"+alt+swegen_dict[new_id][1]+str(AF)+"\n")
                out.write(new_id.split("_")[0]+"\t"+new_id.split("_")[1]+"\t"+new_id.split("_")[2]+"\t"+swegen_dict[new_id][0]+"\t"+alt+"\t"+swegen_dict[new_id][1]+"\t"+str(AF)+"\t"+majority+"\n")




out.close()

