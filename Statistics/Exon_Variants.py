# Exon_Variants.py
# 19.03.24
# Mirjam Müller
#
# Based on the SweGen comparison, we saw that our genotype predictions are mostly accurate for
# variants found in exons. So we filter our genotype output tables to only include variants in
# exons for each gene.
#
# This means first reading in the annotation files, finding the exons and saving them in a dictionary.
#
# Then we open the genotype file, go through the variants, and if the variant is in an exon,
# stop iterating through the annotation and write it directly into an output file.
#
# NOTE THIS IS SUPPOSED TO RUN PER GENE, IT WILL BE TOO MUCH FOR ALL GENES AT ONCE.
#
##################################################################################################

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
                if file == args.gencode:
                    db = "G"
                elif file == args.refseq:
                    db = "R"
                else:
                    db = "You got a bug concerning the db check in database_read function."

                # Add transcript ID to gene dictionary.
                gene_dict[trans_ID] = []

                # make entries for each exon.
                for i in range(0, number_exons):
                    if i == 0:
                        if strand == "+":
                            position = "first"
                        else:
                            position = "last"
                    elif i == number_exons - 1:
                        if strand == "+":
                            position = "last"
                        else:
                            position = "first"
                    else:
                        position = "middle"
                    gene_dict[trans_ID].append([chrom, exon_smaller[i],
                                                exon_bigger[i], strand, position, db])
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

print("Starting Identify AS script! ")
# 1. Process databases in parallel.

# Initialize
gene_dict=dict()

# Find gene name
gene = args.genotype.split("/")[-1].split("_")[0]
print("Creating Database Dictionary...", end="\r")

# Parallelisation so gencode and refseq are read in at the same time.
if __name__=="__main__":
    with Pool(2) as pool:
        # result is a list. i.e. two gene_dicts.
        result=pool.map(database_read, [args.refseq, args.gencode])

# Merge gencode and refseq outputs
gene_dict=merge_nested_dict(result)

print("Creating Database Dictionary: Done! \n", end="\r")

# That gives me the exon coordinates for each transcript in each gene.

#variants = open("ESR1_variants.bed", "w")
#exonvariants = open("ESR1_exon_variants.bed", "w")
# Now we read in the variants, check if they are in an exon and if they are, write them into output.

with open(args.genotype, "r") as genotype, open("Exon_Genotypes/"+gene+"_exgeno.tsv","w") as out:
    for line in genotype:
        if line.startswith("#"):
            continue
        if line.startswith("Location"):
            # thats the column headers
            out.write(line)
            continue

        variant_ID = line.split("\t")[0]
        #print(variant_ID)
        chrom, position = variant_ID.split("_")[0:2]
        #variants.write(chrom+"\t"+str(position)+"\t"+str(int(position)+1)+"\t"+"+"+"\n")

        #Check position against annotation dictionary
        position = int(position)

        exon_found=False
        for trans_ID in gene_dict:
            if exon_found:
                break
            for exon in gene_dict[trans_ID]:
                if position >= int(exon[1]) and position < int(exon[2]):
                    #Is in exon
                    out.write(line)
                    #exonvariants.write(chrom + "\t" + str(position) + "\t" + str(int(position) + 1) + "\t" + "+" + "\n")
                    exon_found=True
                    break

#variants.close()
#exonvariants.close()

#%% Timer

print("Run time: {:.2f} seconds.".format(time.time()-start_time))