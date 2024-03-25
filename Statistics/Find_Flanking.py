# Find_Flanking.py
# 21.03.24
# Mirjam Müller
#
# This is to match variants to each AS event for statistical testing.
# For IR and CE we match to variants in flanking exons. For CE, AD and AA we match variants in same exon.
# (Yes we do both for CE).
########################################################################################################################

# imports

import argparse
import re
import time
from multiprocessing import Pool

# 0.1 Argparse
"Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='Filter Exon Variants',
                                 usage='%(prog)s -o OUTPUT-FILE \
                                     -g GENCODE-FILE -r REFSEQ-FILE',
                                 description="""Per gene, filters the variants to only include
                                 variants found in exons.""")

parser.add_argument('--genotype', '-gt', required=True,
                    help="Input file, containing genotype events for the gene.")
parser.add_argument('--gencode', '-g', required=True,
                    help="""tsv file containing bed file information on 
                    annotated exons from GENCODE39 as well as gene names.""")
parser.add_argument('--refseq', '-r', required=True,
                    help="""tsv file containing bed file information on 
                    annotated exons from RefSeq as well as gene names.""")
parser.add_argument('--psi', '-p', required=True,
                    help="""tsv file containing PSI scores for the gene.""")



args = parser.parse_args()

# 0.2 Functions
def database_read(file):
    gene_chrom = "NA"
    with open(file, "r") as infile:
        for line in infile:
            # To exclude potential title lines/empty lines, formatting mistakes
            # Only takes chr[] and chr[]_random lines, in accordance with bam.
            if re.search(
                    r"(.+)\t([a-z]{3}[X,M,Y]?\d*).*\t(\-?\+?)\t(\d+)\t(\d+)\t\d+\t\d+\t(\d+)\t([\d,]+)\t([\d,]+)\t(.+)",
                    line):
                # specify groups.
                entry = re.search(
                    r"(.+)\t([a-z]{3}[X,M,Y]?\d*).*\t(\-?\+?)\t(\d+)\t(\d+)\t\d+\t\d+\t(\d+)\t([\d,]+)\t([\d,]+)\t(.+)",
                    line)
                # Assign variables to groups
                trans_ID = entry.group(1)
                chrom = entry.group(2)
                strand = entry.group(3)
                gene_name = entry.group(9)

                # Check if transcript belongs to gene in question, if not skip.
                if gene_name != gene:
                    if gene_chrom != "NA":
                        """When we reach a transcript from a different gene that starts 
                        after the "latest" transcript of the gene we look for ended, then we stop looking."""
                        if int(entry.group(4)) > gene_stop:
                            break
                    continue
                else:
                    # Save information on gene, so we can stop iterating through db file when its found.
                    if gene_chrom == "NA":
                        gene_stop = int(entry.group(5))
                        gene_chrom = chrom
                    else:
                        if gene_stop < int(entry.group(5)):
                            gene_stop = int(entry.group(5))

                """To not get caught in exon start/stop, -/+ strand complications,
                coordinates are referred to as bigger and smaller instead."""

                number_exons = int(entry.group(6))
                exon_smaller = entry.group(7).split(",")[0:-1]
                exon_bigger = entry.group(8).split(",")[0:-1]
                gene_name = entry.group(9)


                # make entries for each exon.
                for i in range(0, number_exons):
                    if i == 0:
                        # Theres an exon with bigger coordinates but no smaller
                        if exon_smaller[i]+"_"+exon_bigger[i] not in gene_dict:
                            gene_dict[exon_smaller[i] + "_" + exon_bigger[i]]=[]
                        if number_exons>1:
                            gene_dict[exon_smaller[i]+"_"+exon_bigger[i]].append([exon_smaller[i+1], exon_bigger[i+1]])
                    elif i == number_exons-1:
                        #There is only exons before this one
                        if exon_smaller[i] + "_" + exon_bigger[i] not in gene_dict:
                            gene_dict[exon_smaller[i] + "_" + exon_bigger[i]] = []
                        gene_dict[exon_smaller[i]+"_"+exon_bigger[i]].append([exon_smaller[i-1], exon_bigger[i-1]])
                    else:
                        if exon_smaller[i] + "_" + exon_bigger[i] not in gene_dict:
                            gene_dict[exon_smaller[i] + "_" + exon_bigger[i]] = []
                        #Theres an exon before this on
                        gene_dict[exon_smaller[i] + "_" + exon_bigger[i]].append(
                           [exon_smaller[i - 1], exon_bigger[i - 1]])
                        #and one after
                        gene_dict[exon_smaller[i] + "_" + exon_bigger[i]].append(
                            [exon_smaller[i + 1], exon_bigger[i + 1]])


    return gene_dict


# To merge two nested dictionaries, I wrote a function.
# This atm works for this specific situation but id like it to work for more general as well. Eventually.
def merge_nested_dict(list_of_dicts):
    """
    Merges a list of dicts into a new dictionary.

    Parameters
    ----------
    list_of_dicts : DICT
        List containing all dictionaries to be merged.

    Returns new_dict containing all information from the list of dicts.
    -------

    """
    # Initialize resulting dict
    new_dict = dict()

    for dictionary in list_of_dicts:
        for k, v in dictionary.items():
            # Check if the dictionary is nested:
            # k is the gene name, v is the transcript dictionary.
            if str(type(v)) == "<class 'dict'>":
                # if the nested thing is a dictionary, add them together.
                new_dict[k] = dict()
                for i, j in v.items():
                    new_dict[k][i] = j
            # Just normal 1 level dictionaries. Proceed without complications.
            else:
                new_dict[k] = v
    return new_dict


# 0.3 Timer

start_time=time.time()

# 1. Process databases in parallel.

# Initialize
gene_dict=dict()

# Find gene name
gene = args.genotype.split("/")[-1].split("_")[0]
print(gene)
#print("Creating Database Dictionary...", end="\r")

# Parallelisation so gencode and refseq are read in at the same time.
if __name__=="__main__":
    with Pool(2) as pool:
        # result is a list. i.e. two gene_dicts.
        result=pool.map(database_read, [args.refseq, args.gencode])

# Merge gencode and refseq outputs
gene_dict=merge_nested_dict(result)

#print("Creating Database Dictionary: Done! \n", end="\r")



#print("Finding matches between AS and variants...", end="\r")

# Iterate through PSI file, note that we only need the first column.
AS_dict=dict()
with open(args.psi, "r") as psi:
    for line in psi:
        if line.startswith("Event"):
            continue
        #reset variables
        matching_exons=[]
        infostring= line.split("\t")[0]
        if infostring.startswith("A"):
            AS_type, group, chrom, strand, position = infostring.split("_")
            # Find the exon the AS event is in.
            for exon in gene_dict:
                if int(position) == int(exon.split("_")[0]) or int(position) == int(exon.split("_")[1]):
                    # found exon
                    matching_exons.append([exon.split("_")[0],exon.split("_")[1]])
                    #We dont look in flanking exons for AA and AD.
        else:
            AS_type, chrom, strand, start, stop = infostring.split("_")
            # Find the exon the AS event is in.
            if AS_type == "CE":
                for exon in gene_dict:
                    if start == exon.split("_")[0] and stop == exon.split("_")[1]:
                        #found exon
                        matching_exons.append([start, stop])
                        for matches in gene_dict[exon]:
                            matching_exons.append(matches)
            else:
                for exon in gene_dict:
                    if start==exon.split("_")[1]:
                        #found exon.
                        matching_exons.append([exon.split("_")[0], exon.split("_")[1]])
                        #for IR we only want the exons immediately bordering to the intron.
                        # So given that the start of an intron matched exon end, we only want exons AFTER this exon
                        for matches in gene_dict[exon]:
                            if int(matches[0])>int(start):
                                matching_exons.append(matches)

                    elif stop == exon.split("_")[0]:
                        #found exon
                        matching_exons.append([exon.split("_")[0], exon.split("_")[1]])
                        # only exons BEFORE this exon
                        for matches in gene_dict[exon]:
                            if int(matches[0])<int(stop):
                                matching_exons.append(matches)
        AS_dict[infostring]=matching_exons

#Now that we have the matching exons, we match the variants to the matching exons.
#The file are already only exon variants, so we only need to check which exons they are in :)
#Initialize dictionary to save pairs in.
pairs=dict()
for event in AS_dict:
    pairs[event]=[]
with open(args.genotype, "r") as genotype:
    for line in genotype:
        if line.startswith("#"):
            continue
        elif line.startswith("Location"):
            continue

        locstring = line.split("\t")[0]
        position = int(locstring.split("_")[1])
        for event in AS_dict:
            matching_exons=AS_dict[event]
            for exon in matching_exons:
                if position>=int(exon[0]) and position<=int(exon[1]):
                    #match
                    if event.startswith("A"):
                        for e in AS_dict:
                            if e.startswith(event[0:5]):
                                pairs[e].append(locstring)
                    else:
                        pairs[event].append(locstring)

#Remove duplicates from dictionary
for event in pairs:
    entry = list(set(pairs[event]))
    pairs[event] = entry

#Write into output file
with open("VAR_AS_Matches/"+gene+"_"+"pairs.tsv", "w") as out:
    out.write("AS_event\tvariants\n")
    for event in pairs:
        out.write(event+"\t"+",".join(pairs[event])+"\n")



