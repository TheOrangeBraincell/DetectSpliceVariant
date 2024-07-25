# Make_PSI_Table.py
# 19.03.24
# Mirjam Müller
#
# To do summary statistics, we need to process all the PSI files for all the genes we have run for.
# This turned out to be too much for R. So this script processes all the PSI files and creates one table with it.
# Which is then still not possible to be read into R.
#
#
# imports
import argparse
import os



#argparse
parser = argparse.ArgumentParser(prog='Make PSI summary table',
                                 usage='',
                                 description="""Parses PSI output tables for each gene and puts them into one table.""")

parser.add_argument('--folder', '-f', required=True,
                    help="""input folder of PSI tables""")
parser.add_argument('--out', '-o', required=True,
                    help="""output file.""")

args = parser.parse_args()

files = os.listdir(args.folder)

temp_dict = dict()
counter = 0
with open(args.out, "w") as out:
    out.write("Gene\tEvent\tSample\tPSI\n")
    for file in files:
        with open(args.folder+file, "r") as infile:
            counter += 1
            gene=file.split("_")[0]
            temp_dict[gene] = dict()
            for line in infile:
                if line.startswith("Event"):
                    sample_names = line.strip("\n").split("\t")[1:]
                    continue
                event_id = line.split("\t")[0]
                #add event to nested dictionary
                temp_dict[gene][event_id] = dict()
                psi_scores = line.strip("\n").split("\t")[1:]

                for i in range(0, len(psi_scores)):
                    temp_dict[gene][event_id][sample_names[i]]=psi_scores[i]


        # Write all events of one gene into file
        for event in temp_dict[gene]:
            for sample in temp_dict[gene][event]:
                if temp_dict[gene][event][sample]!="NAN":
                    out.write(gene+"\t"+event+"\t"+sample+"\t"+temp_dict[gene][event][sample]+"\n")

        #Reset dictionary
        temp_dict=dict()
        #Update progress
        print("Processing gene number " + str(counter), end="\r")




