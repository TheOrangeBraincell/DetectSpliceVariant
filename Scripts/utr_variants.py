# -*- coding: utf-8 -*-
"""
Date: Thu Oct 12 15:11:15 2023
File Name: utr_variants.py
Author: Mirjam Müller

Description:
    
    
List of Functions:
    
Procedure: 
    1.
    2.
    3.
    
Useage:
    python ../gitrepo/Scripts/utr_variants.py -v ESR1_locations.tsv -o no_utr_locations.tsv -n ESR1 -g ../Database/GENCODE39.tsv -r ../Database/RefSeq.tsv
    
Possible Bugs:
"""

#%% 0.0 Imports
import argparse
import time
import re
from multiprocessing import Pool

#%% 0.1 Argparse

"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='Filter three prime UTR variants',
                                 usage='',
                                 description=""" Removes variants which can only be found
                                 in the three prime UTR of a gene.""")

parser.add_argument('--variant', '-v', required=True,
                    help='input file containing variant locations.')
parser.add_argument('--out', '-o', required=True,
                    help="Output file, containing variants not in UTRs.")
parser.add_argument('--gene', '-n', required=True,
                    help="Name of the gene we are filtering variants for")
parser.add_argument('--gencode', '-g', required=True,
                    help="""tsv file containing bed file information on 
                    annotated exons from GENCODE39 as well as gene names.""")
parser.add_argument('--refseq', '-r', required=True,
                    help="""tsv file containing bed file information on 
                    annotated exons from RefSeq as well as gene names.""")


args = parser.parse_args()

#%% 0.2 Functions

def database_read(file):
    gene_chrom="NA"
    with open(file, "r") as infile:
        for line in infile:
            # To exclude potential title lines/empty lines, formatting mistakes
            # Only takes chr[] and chr[]_random lines, in accordance with bam.
            if re.search(r"(.+)\t([a-z]{3}[X,M,Y]?\d*).*\t(\-?\+?)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t([\d,]+)\t([\d,]+)\t(.+)", line):
                # specify groups.
                entry = re.search(r"(.+)\t([a-z]{3}[X,M,Y]?\d*).*\t(\-?\+?)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t([\d,]+)\t([\d,]+)\t(.+)", line)
                #Assign variables to groups
                trans_ID=entry.group(1)
                chrom=entry.group(2)
                gene_name=entry.group(11)
                strand=entry.group(3)
                
                #Check if transcript belongs to gene in question, if not skip. 
                if gene_name!= args.gene:
                    if gene_chrom!="NA":
                        """When we reach a transcript from a different gene that starts 
                        after the "latest" transcript of the gene we look for ended, then we stop looking."""
                        if int(entry.group(4))> gene_stop:
                            break
                    continue
                else:
                    #Save information on gene, so we can stop iterating through db file when its found.
                    if gene_chrom=="NA":
                        gene_stop=int(entry.group(5))
                        gene_chrom=chrom
                    else:
                        if gene_stop<int(entry.group(5)):
                            gene_stop=int(entry.group(5))
                
                
                #Add transcript ID to gene dictionary.
                gene_dict[trans_ID]=dict()
                trans_coord=[int(entry.group(4)), int(entry.group(5))]
                cds_coord=[int(entry.group(6)), int(entry.group(7))]
                if strand =="+":
                    gene_dict[trans_ID]["UTR"]=[[cds_coord[1], trans_coord[1]]]
                else:
                    gene_dict[trans_ID]["UTR"]=[[trans_coord[0], cds_coord[0]]]
               
                """For this we need all exon coordinates per transcript. 
                First regardless of strand."""

                number_exons=int(entry.group(8))
                exon_smaller=entry.group(9).split(",")[0:-1]
                exon_bigger=entry.group(10).split(",")[0:-1]
                
                #make entries for each exon.
                gene_dict[trans_ID]["exons"]=[]
                for i in range(0, number_exons):
                    gene_dict[trans_ID]["exons"].append([exon_smaller[i], 
                                                exon_bigger[i]])
    return gene_dict

#To merge two nested dictionaries, I wrote a function.
#This atm works for this specific situation but id like it to work for more general as well. Eventually.
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
    #Initialize resulting dict
    new_dict=dict()
    
    for dictionary in list_of_dicts:
        for k, v in dictionary.items():
            #Check if the dictionary is nested:
            #k is the gene name, v is the transcript dictionary.
            if str(type(v))=="<class 'dict'>":
                #if the nested thing is a dictionary, add them together.
                new_dict[k]=dict()
                for i, j in v.items():
                    new_dict[k][i]=j
            #Just normal 1 level dictionaries. Proceed without complications.
            else:
                new_dict[k]=v      
    return new_dict

#%% 0.3 Start Timer

start_time=time.time()

#%% 1. Read in Gencode and Refseq. 

"""We require the information on end of CDS and end of last exon. per transcript.
Although to check if theres any exons overlapping with any UTR, i guess we need all exon information too.
"""
#Initialize
gene_dict=dict()
print("Creating Database Dictionary...", end="\r")

#Parallelisation so gencode and refseq are read in at the same time.
if __name__=="__main__":
    with Pool(2) as pool:
        #result is a list. i.e. two gene_dicts.
        result=pool.map(database_read, [args.refseq, args.gencode])


#Merge gencode and refseq outputs
gene_dict=merge_nested_dict(result)

print("Creating Database Dictionary: Done! \n", end="\r")


#%% 2. Now we need to find the UTR region coordinates based on the dictionary we just created.


for trans_ID in gene_dict:
    for second_ID in gene_dict:
        #If they are the same, no comparison needed.
        if trans_ID==second_ID:
            continue
        #If UTR is the same, we dont need to check the exons.
        if gene_dict[trans_ID]["UTR"]==gene_dict[second_ID]["UTR"]:
            continue
        
        to_be_removed=[]
        #Otherwise we have two transcripts with different UTRs. 
        #Then we have to check the exons.
        for exon in gene_dict[second_ID]["exons"]:
            #see if in range of UTR
            for utr in gene_dict[trans_ID]["UTR"]:
                if exon[1]<utr[0] or exon[0]>utr[1]:
                    #Not in range
                    continue
                
                #in range
                if exon[0]>utr[0] and exon[1]<utr[1]:
                    #start and stop in UTR
                
                elif exon[0]>utr[0]:
                    #only start in utr
                
                elif exon[1]<utr[1]:
                    #only stop in utr
                
                elif exon==utr:
                    #they have the same coordinates. if last exon, ignore. 
                    #otherwise this utr needs to be ignored instead.
                    if gene_dict[second_ID]["exons"].index(exon)!=len(gene_dict[second_ID]["exons"])-1:
                        #not the last exon.
                        to_be_removed.append(utr)
                    
        new_list_utr=[]
        for utr in gene_dict[trans_ID]["UTR"]:
            if utr not in to_be_removed:
                new_list_utr.append(utr)
        
        gene_dict[trans_ID]["UTR"]=new_list_utr
                
                
                
                
                
                
    

#%% End Timer

print("Run time: {:.2f} seconds.".format(time.time()-start_time))  