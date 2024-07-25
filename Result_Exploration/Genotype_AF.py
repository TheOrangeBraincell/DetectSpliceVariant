# Make_PSI_Table.py
# 27.03.24
# Mirjam Müller
#
# To do an all genes allele frequency comparison between SWEGEN and SCANB, i need to have AF for all genes.
# However, R refuses to read in all the genotype files. this will help with that.
#
# imports
import argparse

#argparse
parser = argparse.ArgumentParser(prog='Make AF summary table',
                                 usage='',
                                 description="""Parses genotype output tables for each gene and puts AF into one table.""")

parser.add_argument('--files', '-f', required=True,
                    help="""file containing all the genotype table names""")
parser.add_argument('--out', '-o', required=True,
                    help="""output file.""")

args = parser.parse_args()

skipped=0
with open(args.files, "r") as infile, open(args.out, "w") as out:
    out.write("Location\tgene\tAllele_Frequency_Observed\tNumber_Samples\n")
    for file in infile:
        with open(file.strip("\n"), "r") as genotype:
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
