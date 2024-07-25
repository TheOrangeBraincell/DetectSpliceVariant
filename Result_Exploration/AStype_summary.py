# AStype_summary.py
# 10.07.24
# Mirjam Müller
#
#
# Description: Parses through PSI_summary file
#               Extracts all Events belonging to one event type
#               Writes them into output file
#
##############################################################################################################3
#
#
#
#
import sys

infile = sys.argv[1]
outfile = sys.argv[2]
as_type = sys.argv[3]

with open(infile, "r") as table, open(outfile, "w") as out:
    for line in table:
        if line.startswith("Gene"):
            #thats the header
            out.write(line)
            continue

        event=line.split("\t")[1].split("_")[0]

        if event == as_type:
            out.write(line)


