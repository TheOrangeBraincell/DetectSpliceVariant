# Genotype_AF.py
# 27.03.24
# Mirjam Müller
#
# Calculates the allele frequencies for a set of variants.
#
# imports
import sys
import os

input_folder = sys.argv[1]
files = os.listdir(input_folder)

variant_file = sys.argv[2]

outfile = sys.argv[3]

locations=[]
with open(variant_file, "r") as variants:
    for line in variants:
        if line.startswith("Location"):
            #header
            continue
        locations.append(line.split("\t")[1])

skipped=0
with open(outfile, "w") as out:
    out.write("Location\tgene\tAllele_Frequency_Observed\tNumber_Samples\n")
    for file in files:
        with open(input_folder + file, "r") as genotype:
            for line in genotype:
                if line.startswith("#"):
                    skipped+=1
                    continue #thats header information
                elif line.startswith("Gene"):
                    skipped+=1
                    continue #we will change column names anyways
                elif line.startswith("Location"):
                    skipped+=1
                    continue
                else:
                    infostring, gene=line.split("\t")[0:2]
                    getypes=line.strip("\n").split("\t")[2:]

                    if infostring in locations:
                        counts={"0/0":0, "0/1":0, "1/1":0}
                        for type in getypes:
                            if type =="NE":
                                #not expressed, move on.
                                continue
                            elif type=="ND":
                                #filtered variant, move on
                                continue
                            elif type=="0/0":
                                counts["0/0"]+=1
                            elif type in ["1/1", "2/2", "3/3"]:
                                counts["1/1"]+=1
                            else:
                                #heterozgous
                                counts["0/1"]+=1

                        #Now that we have the counts for this location, calculate observed allele frequency
                        total=counts["1/1"]+counts["0/0"]+counts["0/1"]
                        observed=(2*counts["1/1"]+counts["0/1"])/(2*total)

                        #Now write it into output file
                        out.write(infostring+"\t"+gene+"\t"+str(observed)+"\t"+str(total)+"\n")


print(skipped)
