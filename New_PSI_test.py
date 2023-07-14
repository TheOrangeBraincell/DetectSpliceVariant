# -*- coding: utf-8 -*-
"""
Date: Fri Jul 14 08:56:18 2023
File Name:
Author: Mirjam Karlsson-Müller

Description:
    
    
List of Functions:
    
Procedure: 
    1.
    2.
    3.
    
Useage:
    
    
Possible Bugs:
"""

#%% 0.0 Imports

import argparse
import re
import time
import pysam
import math
import multiprocessing as mp

#%% 0.1 argparse

"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='Score Alternative Splicing',
                                 usage='%(prog)s -s SAMPLE-FOLDER -o OUTPUT-FOLDER \
                                     -g GENCODE-FILE -r REFSEQ-FILE \
                                         [-c] "chrX:XXXXXX-XXXXXX" -as AS-TYPE \
                                             -is INSERT-SIZE',
                                 description="""Creates
                                 a table with the PSI scores supporting said event
                                 per sample in sample folder.""")

parser.add_argument('--samples', '-s', required=True,
                    help='file containing the paths to the bam files of all samples.')
parser.add_argument('--input', '-i', required=True,
                    help="Table with all identified AS events, created by Identify_AS.py")
parser.add_argument('--out', '-o', required=True,
                    help="""Output file, containing PSI tables.""")
parser.add_argument('--InsertSize', '-is', type=str,
                    help="""Average Insert size plus standard deviation. Format
                    'Mean X Standard Deviation Y' """)

args = parser.parse_args()


#If we need to calculate scores for intron retention, we require an average insert size between the reads.
if args.InsertSize:
    try:
        insert_mean=float(args.InsertSize.split(" ")[1])
        insert_sd=float(args.InsertSize.split(" ")[4].strip("\n"))
    except:
        raise argparse.ArgumentTypeError("""Insert Sizes are either missing or
                                         of wrong format. Required format is:
                                             'Mean [float] Standard Deviation [float]'""")
    
#%% 0.2 Functions    

def Filter_Reads(read, gene_strand):
    skip=False
    #Exclude non-primary alignments
    if read.is_secondary:
        skip=True
    
    #Exclude reads on wrong strand
    if read.mate_is_reverse and read.is_read1:
        read_strand="-"
    elif read.mate_is_reverse and read.is_read2:
        read_strand="+"
    elif read.mate_is_forward and read.is_read1:
        read_strand="+"
    elif read.mate_is_forward and read.is_read2:
        read_strand="-"
    
    if read_strand!=gene_strand:
        skip=True
        
    return skip


#%% 0.3 Start Timer

start_time=time.time()

#%% 1. Read in AS events

#We will save the identifiers in this dict and then add the regions that should be counted.
AS_events=dict()

with open(args.input, "r") as asfile:
    for line in asfile:
        if line.startswith("Location"):
            #Thats the file header, we ignore it
            continue
        if line.startswith("#"):
            #This is the gene header. Print it to log so we can double check that we run for right gene.
            print(line.strip())
            #list of gene information: gene name, chrom, strand, min, max
            gene=[e.strip(" ") for e in line.strip("\n").split(",")]
            gene[3]=int(gene[3])
            gene[4]=int(gene[4])
            continue
        
        #All thats left now is AS entries.
        Location, AS_Type, Gene_name, Database, Gencode_ID, Refseq_ID= line.strip().split("\t")
        #Initialize dictionary
        AS_events[AS_Type+"_"+Location]=dict()

#%% 2. Read in list of bam files. 

print("Reading in BAM files...", end="\r")

def PSI_for_Sample(sample):
    #Gene information to open bam file.
    chrom=gene[1]
    gene_start=gene[3]
    gene_stop=gene[4]
    gene_strand=gene[2]
    
    #initialize counting dictionary for this sample
    for event in AS_events:
        if event.startswith("CE"):
            AS_events[event]={"ER":{"count":0, "reads":[]}, 
                                         "IR":{"count":0, "reads":[], "count_n":0, 
                                               "normalize":int(event.split("_")[4])-int(event.split("_")[3])}}
    
    #open the bam file, and go through the reads.
    samfile=pysam.AlignmentFile(files[sample], 'rb', index_filename=files[sample][0:-1]+"i")
    reads=samfile.fetch(chrom, gene_start, gene_stop)
    
    #Iterate through reads
    for read in reads:
        #Filter
        if Filter_Reads(read, gene_strand)==True:
            continue
    
        #read information
        read_start=int(read.reference_start)
        read_range=sum([int(i) for i in re.findall(r'\d+', read.cigarstring)])
        read_stop=read_start+read_range
        
        #Now we go through all the events and try to figure out what to count them for.
        for event in AS_events:
            #event=event id later used in first column of output.
            #Extract CE info from event key.
            if event.startswith("CE"):
                CE_start=int(event.split("_")[3])
                CE_stop=int(event.split("_")[4])
                
                #If both start and stop of read are before or after CE, then we can go on.
                if read_start < CE_start and read_stop < CE_start:
                    continue
                if read_start > CE_stop and read_stop > CE_stop:
                    continue
                
                #We are in the right range, lets process the read.
                current_cigar = read.cigarstring
                read_name=read.query_name
                #Spliced or not?
                if re.search(r'\d+M(?:\d+[I,D,S,H])?\d+N(?:\d+[I,D,S,H])?\d+M',read.cigarstring):
                    #This is spliced, but its could have insertions and/or deletions around the intron... 
                    if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
                        #we save excluded weird reads like this one for the logfile.
                        excluded_reads.add(read_name)
                        continue
                    
                    current_start = int(read_start)
                    #Otherwise we process it as spliced
                    while re.search(r'(\d+)M(\d+)N(\d+)M', current_cigar):
                        #assign splice junction variables
                        junction = re.search(r'(\d+)M(\d+)N(\d+)M', current_cigar)
                        
                        #assign variables to groups.
                        exon1 = int(junction.group(1))
                        intron=int(junction.group(2))
                        exon2=int(junction.group(3))
                        exon1_start = current_start
                        exon1_end = exon1_start+exon1+1  #exclusive
                        exon2_start = exon1_end+intron -1 #inclusive
                        exon2_end=exon2_start+exon2+1
                            
                        #skip alignments with less than 3 matching bases in an exon.
                        if exon1<3 or exon2<3:
                            # update cigar string
                            current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
                            current_start= exon2_start
                            continue
                        
                        #Number of spliced reads from previous to CE.
                        if exon1_end< CE_start and CE_stop >= exon2_end>=CE_start:
                            #count!
                            if read_name not in AS_events[event]["IR"]["reads"]:
                                AS_events[event]["IR"]["count"]+=1
                                AS_events[event]["IR"]["reads"].append(read_name)
                            
                        #Number of spliced reads from CE to next.   
                        elif CE_start <= exon1_start<=CE_stop and exon2_start >CE_stop:
                            #count!
                            if read_name not in AS_events[event]["IR"]["reads"]:
                                AS_events[event]["IR"]["count"]+=1
                                AS_events[event]["IR"]["reads"].append(read_name)
                                
                        #Number of spliced reads accross CE.
                        elif exon1_end<CE_start and exon2_start>CE_stop:
                            #count!
                            if read_name not in AS_events[event]["ER"]["reads"]:
                                AS_events[event]["ER"]["count"]+=1
                                AS_events[event]["ER"]["reads"].append(read_name)
                                
                        # update cigar string
                        current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
                        current_start= exon2_start
                
                #Otherwise the read is not spliced and we need to count it if its within the CE.
                else:
                    #read information
                    read_start=int(read.reference_start)
                    read_range=sum([int(i) for i in re.findall(r'\d+', read.cigarstring)])
                    read_stop=read_start+read_range
                    read_name=read.query_name
                    
                    #Since we already excluded all reads in the wrong range and spliced are counted, we can just count these.
                    if read_name not in AS_events[event]["IR"]["reads"]:
                        AS_events[event]["IR"]["count_n"]+=1
                        AS_events[event]["IR"]["reads"].append(read_name)
                        
            elif event.startswith("AA"):
                
                
                
                continue
            elif event.startswith("AD"):
                continue
            elif event.startswith("IR"):
                continue
            
    #Now we have gone through all the reads, we can return the PSI scores.
    sample_PSI={"sample":sample}
    for event in AS_events:
        if event.startswith("CE"):
            print(sample, event, AS_events[event]["IR"]["count"], AS_events[event]["IR"]["count_n"], AS_events[event]["ER"]["count"])
            #PSI is NAN if there is 10 or less reads.
            if AS_events[event]["IR"]["count"]+AS_events[event]["IR"]["count_n"]+AS_events[event]["ER"]["count"]<11:
                PSI="NAN"
            else:
                #sum IR, normalize count_n
                IR=AS_events[event]["IR"]["count"]+(AS_events[event]["IR"]["count_n"]/AS_events[event]["IR"]["normalize"])
                #assign ER
                ER=AS_events[event]["ER"]["count"]
                PSI=str(round(IR/(IR+ER), 3))
        
        elif event.startswith("AA"):
            PSI="NAN"
        elif event.startswith("AD"):
            PSI="NAN"
        elif event.startswith("IR"):
            PSI="NAN"
            
        sample_PSI[event]=PSI

    return sample_PSI

bam_file_list=[]
with open(args.samples, "r") as infile:
    for line in infile:
        bam_file_list.append(line.strip())

#Assuming that each bam.file is in its own sample folder with the sample name as folder name. 
sample_names=sorted(list(set([i.split("/")[-4] for i in bam_file_list])))

#make that a dictionary.
files=dict()
for sample in sample_names:
    for file in bam_file_list:
        if sample in file:
            files[sample]=file

print("Reading in BAM files: Done! \n", end="\r")

#%% 3. Read counts, per sample, allow parallelization

excluded_reads=set()

with mp.Pool(3) as pool:
    #result is a list. i.e. two gene_dicts.
    result=pool.map(PSI_for_Sample, sample_names)

with open(args.out, "w") as outfile:
    #write header
    outfile.write("Event\t"+"\t".join(sample_names)+"\n")
    for event in AS_events:
        new_line=event
        for sample in sample_names:
            #find right scores in PSI_dict
            for dictionary in result:
                if dictionary["sample"]==sample:
                    new_line+="\t"+ dictionary[event]
        outfile.write(new_line+"\n")
    

#%% 4. Write output




