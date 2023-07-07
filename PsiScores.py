# -*- coding: utf-8 -*-
"""
Date: Wed Jul  5 11:06:17 2023
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

#%% Imports

import argparse
import glob
import re
import time
import pysam
import math
from inspect import currentframe, getframeinfo

#%% Time

start_time=time.time()

#%% 0. argparse

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
                    help='folder containing sample folders containing among \
                        others, the vcf, bam and gene.tsv files.')
parser.add_argument('--input', '-i', required=True,
                    help="Table with all identified AS events, created by Identify_AS.py")
parser.add_argument('--out', '-o', required=True,
                    help="""Output file, containing PSI tables.""")
parser.add_argument('--InsertSize', '-is', type=str,
                    help="""Average Insert size plus standard deviation. Format
                    'Mean X Standard Deviation Y' """)
parser.add_argument('--coordinates', '-c', type=str,
                    help="""Start and stop coordinates of region of interest,
                    as well as chromosome. Format: chr[]:start-stop""")
parser.add_argument('--AS', '-as', type=str, required=True,
                    help="""Which type of alternative splicing event we are
                    interested in. "CE" for Casette Exons, "AA" for alternative
                    acceptors, "AD" for alternative donors, "IR" for intron
                    retention and "ALL" for all of the types. Several seperated
                    by ,.""")                    
                    


args = parser.parse_args()

# Extract input coordinates, check their format.
if args.coordinates:
    if re.search(r'[a-z]{3}[MXY]?\d*:\d+-\d+', args.coordinates):
        coord = re.search(r'([a-z]{3}[MXY]?\d*):(\d+)-(\d+)', args.coordinates)
        coord_chrom = coord.group(1)
        #To adjust for different inclusivity of stop/start for plus and minus strand, expand coordinate range by 1.
        coord_start = int(coord.group(2))-1
        coord_stop = int(coord.group(3))+1

    else:
        raise argparse.ArgumentTypeError("""The given coordinates are in the 
                                         wrong format. Input as 
                                         chrX:XXXX-XXXX.""")
        quit()


if args.AS:
    allowed_inputs=["CE", "AA", "AD", "IR", "ALL"]
    inputs=args.AS.split(",")
    inputs=[i.upper() for i in inputs]
    for i in inputs:
        if i not in allowed_inputs:
            raise argparse.ArgumentTypeError("""The allowed abbreviations for
                                             splicing events are "CE" for 
                                             Casette Exons, "AA" for 
                                             alternative
                                             acceptors, "AD" for alternative 
                                             donors, "IR" for intron
                                             retention and "ALL" for all of 
                                             the types. Several seperated
                                             by ,.""")
            quit()
            
    #If input is all events, make sure the code runs through all:
    if inputs[0].upper()=="ALL":
        inputs=["CE", "AA", "AD", "IR"]
    #if IR is in the events, but CE is not, then CE still has to run! Because needed for PSI of IR
    if "IR" in inputs and "CE" not in inputs:
        print("Casette exons (CE) will also be identified as needed for PSI score calculation of intron retention (IR).")
        inputs.append("CE")

#If we need to calculate scores for intron retention, we require an average insert size between the reads.
if args.InsertSize:
    try:
        insert_mean=float(args.InsertSize.split(" ")[1])
        insert_sd=float(args.InsertSize.split(" ")[4].strip("\n"))
    except:
        raise argparse.ArgumentTypeError("""Insert Sizes are either missing or
                                         of wrong format. Required format is:
                                             'Mean [float] Standard Deviation [float]'""")
    
#%% User Defined Functions    

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

def PSI_for_Gene(gene, events, sample):
    """
    Calculates PSI score for all events in gene and one sample.

    Parameters
    ----------
    gene : TYPE
        DESCRIPTION.
    events : TYPE
        DESCRIPTION.
    sample : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    #Iterate through event types, as they require different calculations.
    psi_scores=[]
    #In the end returns PSI scores for all events in this gene for a given sample.
    for AS in events:
        for event in events[AS]:
            if AS=="AA":
                #There is several starts for each AA, so there will be several PSI
                PSI=PSI_AA(gene, sample, events[AS][event])
                for i in PSI:
                    psi_scores.append(i)
            elif AS=="AD":
                PSI=PSI_AD(gene, sample, events[AS][event])
                for i in PSI:
                    psi_scores.append(i)
            elif AS=="CE":
                PSI=PSI_CE(sample, event, gene)
                if PSI=="1.0":
                    CE_PSI[sample].append(event)
                psi_scores.append(PSI)
            else:
                PSI=PSI_IR(sample, event, gene)
                psi_scores.append(PSI)
                
    return psi_scores
    
def PSI_CE(sample, CE, gene):
    """
    

    Parameters
    ----------
    sample : string
        Sample name
    exon : list
        contains chromosome, start, stop, strand information for CE

    Returns
    -------
    str
        PSI score

    """
    #Needs coordinates for CE
    chrom=CE[0]
    start=int(CE[2])
    stop=int(CE[3])
    strand=CE[1]                                                                        
    
    #Find .bam file corresponding to sample name.
    for file in bam_file_list:
        if sample in file:
            bam_file=file
            index_file=file[0:-1]+"i"
            break
    
    #We only want to count each read of a read pair once. So we have a list of reads that have already been counted.
    counted=[]
    #Note that we prioritize counting spliced reads, so those are checked first.
    
    #Initialize opening of file
    samfile=pysam.AlignmentFile(bam_file, 'rb', index_filename=index_file)
    
    #For the spliced reads, we need to fetch a bigger region and then recognise the spliced alignments.
    #We dont want to miss any reads, so we make this region the gene range the exon is in.
    spliced_reads=samfile.fetch(chrom, int(current_gene[3]), int(current_gene[4]))
    
    #initialize counters
    counter_left=0
    counter_right=0
    counter_accross=0
    for read in spliced_reads:
        #Filter
        if Filter_Reads(read, strand)==True:
            continue
        #only spliced reads
        if not re.search(r'\d+M(?:\d+[I,D,S,H])?\d+N(?:\d+[I,D,S,H])?\d+M',read.cigarstring):
            #This is spliced, but its could have insertions and/or deletions around the intron... 
            if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
                excluded_reads.add(read.query_name)
                continue
        #Allow for several splice junctions in one read.
        current_cigar = read.cigarstring
        read_start=int(read.reference_start)
        #sum all numbers in cigar string for read length
        read_range=sum([int(i) for i in re.findall(r'\d+', read.cigarstring)])
        current_start = int(read_start)
        #Do not process spliced reads in the wrong place.
        if read_start < start and read_start+read_range < start:
            continue
        if read_start > stop and read_start+read_range > stop:
            continue
        read_name=read.query_name
    
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
            if exon1_end< start and stop >= exon2_end>=start:
                #count!
                if read_name not in counted:
                    counter_left+=1
                    counted.append(read_name)
                else:
                    pass
            
            #Number of spliced reads from CE to next.   
            elif start <= exon1_start<=stop and exon2_start >stop:
                #count!
                if read_name not in counted:
                    counter_right+=1
                    counted.append(read_name)
                else:
                    pass
                    
            #Number of spliced reads accross CE.
            elif exon1_end<start and exon2_start>stop:
                #count!
                if read_name not in counted:
                    counter_accross+=1
                    counted.append(read_name)
                else:
                    pass
                    
            # update cigar string
            current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
            current_start= exon2_start
    
    #C= number of reads within CE.
    #For this we fetch the part of the alignment file in the exon.
    exon_reads=samfile.fetch(chrom, start, stop)
    #The coordinates do not need to be adjusted as pysam uses also start incl and end excl and 0 based.
    counter=0
    for read in exon_reads:
        #Filter
        if Filter_Reads(read, strand)==True:
            continue
        #no spliced reads
        if re.search(r'\d+M(?:\d+[I,D,S,H])?\d+N(?:\d+[I,D,S,H])?\d+M',read.cigarstring):
            if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
                excluded_reads.add(read.query_name)
            continue
        read_name=read.query_name
        if read_name not in counted:
            counter+=1
    
    #counter needs to be normalized with length of the exon.
    counter_n=counter/(stop-start)
    
    IR= counter_left+counter_right+counter_n
    ER= counter_accross
    #minus normalized, plus not normalized
    if ER+ IR-counter_n+counter >10:
        PSI=str(round(IR/(IR+ER),3))
    else:
        PSI="NAN"
     
    return PSI


def PSI_AA(gene, sample, event):
    chrom = current_gene[1]
    strand=current_gene[2]
    gene_start=int(current_gene[3])
    gene_stop=int(current_gene[4])
    
    #Find .bam file corresponding to sample name.
    for file in bam_file_list:
        if sample in file:
            bam_file=file
            index_file=file[0:-1]+"i"
            break
    #Get all starts
    starts=[]
    for i in event:
        starts.append(int(i[2]))
        
    #Sort starts by coordinate.
    starts=sorted(starts)
    
    #We only want to count each read of a read pair once per event. So we have a list of reads that have already been counted.
    counted=[]
    #Note that we prioritize counting spliced reads, so those are checked first.
    
    #Spliced Counters
    spliced_counters=[0]*len(event)
    
    #Initialize opening of file
    samfile=pysam.AlignmentFile(bam_file, 'rb', index_filename=index_file)
    reads=reads=samfile.fetch(chrom, gene_start, gene_stop)
    for read in reads:
        #Filter
        if Filter_Reads(read, strand)==True:
            continue
        #only want spliced reads for this
        if not re.search(r'\d+M(?:\d+[I,D,S,H])?\d+N(?:\d+[I,D,S,H])?\d+M',read.cigarstring):
            #This is spliced, but its could have insertions and/or deletions around the intron... 
            if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
                excluded_reads.add(read.query_name)
                continue
        
        #The rest can be counted for spliced.
        #Allow for several splice junctions in one read.
        current_cigar = read.cigarstring
        current_start = int(read.reference_start)
        read_name=read.query_name
        
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
                
            #skip alignments with less than 3 matching bases in an exon.
            if exon1<3 or exon2<3:
                # update cigar string
                current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
                current_start= exon2_start
                continue
            
            #Update counters based on matching starts
            if strand=="+":
                #Then the AA starts have to match the exon2_starts
                for i in range(0, len(starts)):
                    if starts[i]==exon2_start:
                        if read_name not in counted:
                            spliced_counters[i]+=1
                            counted.append(read_name)
            else:
                #Then the AA starts have to match the exon1_ends
                for i in range(0, len(starts)):
                    #Same as for AD on the plus strand, the exon end has to be adjusted, because we calculate it to be exclusive.
                    if starts[i]==exon1_end-1:
                        if read_name not in counted:
                            spliced_counters[i]+=1
                            counted.append(read_name)
                        

            # update cigar string
            current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
            current_start= exon2_start
    
    #Now count the difference reads!
    difference_counters=[]
    difference_lengths=[] #For normalization
    for i in range(0, len(starts)-1):
        start1=starts[i]
        start2=starts[i+1]
        length=start2-start1 #since they are sorted.
        difference_lengths.append(length)
        difference_reads=samfile.fetch(chrom, start1, start2)
        counter=0
        
        #Filter reads
        for read in difference_reads:
            #Filter
            if Filter_Reads(read, strand)==True:
                continue
            #no spliced reads
            if re.search(r'\d+M(?:\d+[I,D,S,H])?\d+N(?:\d+[I,D,S,H])?\d+M',read.cigarstring):
                #This is spliced, but its could have insertions and/or deletions around the intron... 
                if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
                    excluded_reads.add(read.query_name)
                continue
            
            #Check if already counted
            read_name=read.query_name
            if read_name not in counted:
                #if not counted, add
                counter+=1
                counted.append(read_name)
            
        #Add difference reads into list
        difference_counters.append(counter)
        difference_lengths.append(length)
    
    #The difference reads need to be normalized
    normalized_diff=[]
    for i in range(0, len(difference_counters)):
        normalized_diff.append(difference_counters[i]/difference_lengths[i])
    
    #Now calculate PSI for all starts with the counts.
    psi=[]
    #give event group a number for AA_counts
    n=len(AA_counts[sample].keys())+1
    AA_counts[sample][n]=dict()
    for s in starts:
        start_number=starts.index(s)
        AA_counts[sample][n][s]=spliced_counters[start_number]
        if strand=="+":
            if start_number!=0:
                total_reads=difference_counters[start_number-1]+spliced_counters[start_number]+sum([x for i,x in enumerate(difference_counters) if i!=(start_number-1)])+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)]) 
                #only numerical PSI scores for total reads>10
                if total_reads<=10:
                    PSI="NAN"
                else:
                    IR=normalized_diff[start_number-1]+spliced_counters[start_number]
                    ER=sum([x for i,x in enumerate(normalized_diff) if i!=(start_number-1)])+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)]) 
                    PSI=str(round(IR/(IR+ER), 3))
            else:
                total_reads=spliced_counters[start_number]+sum(difference_counters)+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)]) 
                if total_reads<=10:
                    PSI="NAN"
                else:
                    IR=spliced_counters[start_number]
                    ER= sum(normalized_diff)+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)])
                    PSI=str(round(IR/(IR+ER), 3))
        else:
            if start_number!=(len(starts)-1):
                total_reads=difference_counters[start_number-1]+spliced_counters[start_number]+sum([x for i,x in enumerate(difference_counters) if i!=(start_number-1)])+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)])
                if total_reads<=10:
                    PSI="NAN"
                else:
                    IR=normalized_diff[start_number-1]+spliced_counters[start_number]
                    ER=sum([x for i,x in enumerate(normalized_diff) if i!=(start_number-1)])+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)])
                    PSI=str(round(IR/(IR+ER), 3))
            else:
                total_reads=spliced_counters[start_number]+sum(difference_counters)+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)])
                if total_reads<=10:
                    PSI="NAN"
                else:
                    IR=spliced_counters[start_number]
                    ER= sum(difference_counters)+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)])
                    PSI=str(round(IR/(IR+ER), 3))
        
        psi.append(PSI)
    #The psi scores could potentially be in a different order than the starts in the events, since we sorted. So lets fix that.
    psi_sorted=[""]*len(starts)
    for s in starts:
        for i in range(0, len(event)):
            if str(s) in event[i]:
                psi_sorted[i]=psi[starts.index(s)]
                break
                
    return psi_sorted
    
def PSI_AD(gene, sample, event):
    chrom = current_gene[1]
    strand=current_gene[2]
    gene_start=int(current_gene[3])
    gene_stop=int(current_gene[4])
    
    #Find .bam file corresponding to sample name.
    for file in bam_file_list:
        if sample in file:
            bam_file=file
            index_file=file[0:-1]+"i"
            break
    #Get all starts
    stops=[]
    for i in event:
        #print(event, i)
        stops.append(int(i[2]))
        
        
    #Sort starts by coordinate.
    stops=sorted(stops)
    
    #We only want to count each read of a read pair once per event. So we have a list of reads that have already been counted.
    counted=[]
    
    #Spliced Counters
    spliced_counters=[0]*len(event)
    #open bam file.
    samfile=pysam.AlignmentFile(bam_file, 'rb', index_filename=index_file)
    spliced_reads=samfile.fetch(chrom, gene_start, gene_stop)
    
    #Go through reads
    for read in spliced_reads:
        #Filter
        if Filter_Reads(read, strand)==True:
            continue
        #only spliced reads
        if not re.search(r'\d+M(?:\d+[I,D,S,H])?\d+N(?:\d+[I,D,S,H])?\d+M',read.cigarstring):
            #This is spliced, but its could have insertions and/or deletions around the intron... 
            if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
                excluded_reads.add(read.query_name)
                continue
        
        #Allow for several splice junctions in one read.
        current_cigar = read.cigarstring
        current_start = int(read.reference_start)
        read_name=read.query_name
        
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
                
            #skip alignments with less than 3 matching bases in an exon.
            if exon1<3 or exon2<3:
                # update cigar string
                current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
                current_start= exon2_start
                continue
            
            #Update counters based on matching starts
            if strand=="+":
                #Then the AD stops have to match the exon1_stops
                for i in range(0, len(stops)):
                    #The coordinates on plus strand seem to be shifted for AD, by 1. Because we calculate it to be exclusive.
                    if str(stops)==str(exon1_end-1):
                        if read_name not in counted:
                            spliced_counters[i]+=1
                            counted.append(read_name)
                
            else:
                #Then the AD starts have to match the exon2_starts
                for i in range(0, len(stops)):
                    if str(stops[i])==str(exon2_start):
                        if read_name not in counted:
                            spliced_counters[i]+=1
                            counted.append(read_name)
            
            # update cigar string
            current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
            current_start= exon2_start
    
    #Now we need to still count the unspliced reads in the difference regions
    difference_counters=[]
    difference_lengths=[] #For normalization
    for i in range(0, len(stops)-1):
        stop1=stops[i]
        stop2=stops[i+1]
        length=stop2-stop1 #since they are sorted.
        difference_lengths.append(length)
        difference_reads=samfile.fetch(chrom, stop1, stop2)
        counter=0
        
        #Filter reads
        for read in difference_reads:
            #Filter
            if Filter_Reads(read, strand)==True:
                continue
            #no spliced reads
            if re.search(r'\d+M(?:\d+[I,D,S,H])?\d+N(?:\d+[I,D,S,H])?\d+M',read.cigarstring):
                #This is spliced, but its could have insertions and/or deletions around the intron... 
                if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
                    excluded_reads.add(read.query_name)
                continue
            
            #Check if already counted
            read_name=read.query_name
            if read_name not in counted:
                #if not counted, add
                counter+=1
                counted.append(read_name)
            
        #Add difference reads into list
        difference_counters.append(counter)
        difference_lengths.append(length)
        
    #The difference reads need to be normalized
    normalized_diff=[]
    for i in range(0, len(difference_counters)):
        normalized_diff.append(difference_counters[i]/difference_lengths[i])
    
    #Now calculate PSI for all starts with the counts.
    psi=[]
    #give event group a number for AD_counts
    n=len(AD_counts[sample].keys())+1
    AD_counts[sample][n]=dict()
    for s in stops:
        stop_number=stops.index(s)
        AD_counts[sample][n][s]=spliced_counters[stop_number]
        if strand=="+":
            if stop_number!=0:
                total_reads=difference_counters[stop_number-1]+spliced_counters[stop_number]+sum([x for i,x in enumerate(difference_counters) if i!=(stop_number-1)])+sum([x for i,x in enumerate(spliced_counters) if i!=(stop_number)]) 
                #only numerical PSI scores for total reads>10
                if total_reads<=10:
                    PSI="NAN"
                else:
                    IR=normalized_diff[stop_number-1]+spliced_counters[stop_number]
                    ER=sum([x for i,x in enumerate(normalized_diff) if i!=(stop_number-1)])+sum([x for i,x in enumerate(spliced_counters) if i!=(stop_number)]) 
                    PSI=str(round(IR/(IR+ER), 3))
            else:
                total_reads=spliced_counters[stop_number]+sum(difference_counters)+sum([x for i,x in enumerate(spliced_counters) if i!=(stop_number)])
                if total_reads<=10:
                    PSI="NAN"
                else:
                    IR=spliced_counters[stop_number]
                    ER= sum(normalized_diff)+sum([x for i,x in enumerate(spliced_counters) if i!=(stop_number)]) 
                    PSI=str(round(IR/(IR+ER), 3))
        else:
            if stop_number!=(len(stops)-1):
                total_reads=difference_counters[stop_number]+spliced_counters[stop_number]+sum([x for i,x in enumerate(difference_counters) if i!=(stop_number)])+sum([x for i,x in enumerate(spliced_counters) if i!=(stop_number)]) 
                if total_reads<=10:
                    PSI="NAN"
                else:
                    IR=difference_counters[stop_number]+spliced_counters[stop_number]
                    ER=sum([x for i,x in enumerate(normalized_diff) if i!=(stop_number)])+sum([x for i,x in enumerate(spliced_counters) if i!=(stop_number)]) 
                    PSI=str(round(IR/(IR+ER), 3))
            else:
                total_reads=spliced_counters[stop_number]+ sum(difference_counters)+sum([x for i,x in enumerate(spliced_counters) if i!=(stop_number)]) 
                if total_reads<=10:
                    PSI="NAN"
                else:
                    IR=spliced_counters[stop_number]
                    ER= sum(normalized_diff)+sum([x for i,x in enumerate(spliced_counters) if i!=(stop_number)]) 
                    PSI=str(round(IR/(IR+ER), 3))
        
        psi.append(PSI)
        
    #The psi scores could potentially be in a different order than the starts in the events, since we sorted. So lets fix that.
    psi_sorted=[""]*len(stops)
    for s in stops:
        for i in range(0, len(event)):
            if str(s) in event[i]:
                psi_sorted[i]=psi[stops.index(s)]
                break

    return psi_sorted

def PSI_IR(sample, entry, gene):
    """

    Parameters
    ----------
    sample : string
        Sample name
    entry: string separated by _
        contains strand, chromosome, exon1stop, exon2start

    Returns
    -------
    str
        PSI score

    """
    chrom, strand, exon1stop, exon2start=entry
    exon1stop=int(exon1stop)
    exon2start=int(exon2start)
    #Find .bam file corresponding to sample name.
    for file in bam_file_list:
        if sample in file:
            bam_file=file
            index_file=file[0:-1]+"i"
            break
    
    #Check: If there is a CE in the IR we need to check its PSI score
    #Extract start and stop
    for entry in CE_PSI[sample]:
        #The dictionary contains only the CE that have PSI 1 
        CE_start=int(entry[2])
        CE_stop=int(entry[3])
        #Check if within intron
        if CE_start>exon1stop and CE_stop<exon2start:
            #If this PSI=1
            PSI="0.0"
            return PSI
    
    #Insert read confidence interval
    insert_confidence=[insert_mean-1.96*(insert_sd/math.sqrt(len(sample_names))),insert_mean+1.96*(insert_sd/math.sqrt(len(sample_names)))]
    
    #We only want to count each read of a read pair once per event. So we have a list of reads that have already been counted.
    counted=dict()
    #Note that we prioritize counting spliced reads, so those are checked first.
    #Initialize opening of file
    samfile=pysam.AlignmentFile(bam_file, 'rb', index_filename=index_file)
    
    #Check if there is AD/AA events involving the IR coordinates
    AA_AD=False
    spliced_counter=0
    a=0
    b=0
    for event in AA_counts[sample]:
        sum_counts=0
        for start in AA_counts[sample][event]:
            sum_counts+=AA_counts[sample][event][start]
            if start ==exon2start or start== exon1stop:
                #We have an AA event overlapping with the intron to be scored. Use splice counts from that.
                AA_AD=True
                spliced_counter+=sum_counts
                #get coordinates for overlap reads.
                if strand=="+":
                    a=min(AA_counts[sample][event].keys())
                else:
                    a=max(AA_counts[sample][event].keys())
        #if event is found, stop iterating.
        if a!=0:
            break
                
    for event in AD_counts[sample]:
        sum_counts=0
        for stop in AD_counts[sample][event]:
            sum_counts+=AD_counts[sample][event][stop]
            if stop==exon2start or stop==exon1stop:
                #We have an AD event overlapping with the intron to be scored. Use splice counts from that.
                spliced_counter+=sum_counts
                AA_AD=True
                #get coordinates for overlap reads
                if strand=="+":
                    b=max(AD_counts[sample][event].keys())
                else:
                    b=min(AD_counts[sample][event].keys())
        #can stop iterating through events if event is found
        if b!=0:
            break
                
    
    #If there is no overlap with an AD or AA event, then we need to manually count spliced reads.
    if AA_AD==False:
        #Get Spliced Reads
        #Open bam file in gene range.
        gene_start=int(current_gene[3])
        gene_stop=int(current_gene[4])
        spliced_reads=samfile.fetch(chrom, gene_start, gene_stop)
        
        for read in spliced_reads:
            #Filter
            if Filter_Reads(read, strand)==True:
                continue
            #only spliced reads
            if not re.search(r'\d+M(?:\d+[I,D,S,H])?\d+N(?:\d+[I,D,S,H])?\d+M',read.cigarstring):
                #This is spliced, but its could have insertions and/or deletions around the intron... 
                if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
                    excluded_reads.add(read.query_name)
                    continue
    
            read_start=int(read.reference_start)
            read_length=int(read.infer_query_length())
            read_name=read.query_name
            #Allow for several splice junctions in one read.
            current_cigar = read.cigarstring
            current_start = int(read_start)
            
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
                exon2_end=exon2_start+exon2
                    
                #skip alignments with less than 3 matching bases in an exon.
                if exon1<3 or exon2<3:
                    # update cigar string
                    current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
                    current_start= exon2_start
                    continue
                
                #Check if the exons match the IR event
                #print(exon1_start)
                if exon1_start< exon2start and exon2_end>exon1stop:
                    if read_name not in counted:
                        spliced_counter+=1
                        counted[read_name]="spliced"
    
                # update cigar string
                current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
                current_start= exon2_start
    #If there is an AA, AD event, we can use their splice counts
    else:
        spliced_counter+=sum_counts
        
        
    #To count overlapping reads, we need to know which points they need to overlap with:
    #If there is an AA/AD event then there might be additional coordinates to check apart from the IR.
    if AA_AD==True:
        points=dict()
        intron_start=exon1stop
        intron_stop=exon2start
        #If a is not 0, then there is an AA event
        if a!=0:
            #If the strand is +, then a is on the right of the IR
            if strand=="+":
                points[a]="right"
                if exon2start != a:
                    intron_stop=a
                    points[exon2start]="right"
            #if the strand is -, then a is on the left of the IR.
            else:
                points[a]="left"
                if exon1stop!=a:
                    intron_start=a
                    points[exon1stop]="left"
        #If b is not 0, then there is an AD event.
        if b!=0:
            #if the strand is +, then b is on the left of the IR
            if strand=="+":
                points[b]="left"
                if exon1stop !=b:
                    intron_start=b
                    points[exon1stop]="left"
            if strand=="-":
                points[b]="right"
                if exon2start != b:
                    intron_stop=b
                    points[exon2start]="right"

    #If there is no AA/AD event, then we have the usual IR start and stop
    else:
        points={exon1stop:"left", exon2start:"right"}
        intron_start=exon1stop
        intron_stop=exon2start
    
    #Only score IR which have overlapping reads on both sides.
    left=False
    right=False
    overlap_counter=0
    #extract reads
    for p in points:
        #only fetch reads that overlap with the point coordinate
        reads_genome=samfile.fetch(chrom,p, p+1)
        for read in reads_genome:
            #Filter
            if Filter_Reads(read, strand)==True:
                continue
            #no spliced reads
            if re.search(r'\d+M(?:\d+[I,D,S,H])?\d+N(?:\d+[I,D,S,H])?\d+M',read.cigarstring):
                #This is spliced, but its could have insertions and/or deletions around the intron... 
                if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
                    excluded_reads.add(read.query_name)
                continue
            
            #if distance between reads too far, then this is a sign that the second read is in an exon and in between is a break. So no IR.
            if read.tlen>insert_confidence[1]:
                continue
            read_start=int(read.reference_start)
            read_length=int(read.infer_query_length())
            #Since we can have AD, AA events, we check if the read pair is within the intron (not an AA, AD)
            if read_start< intron_start or read_start+read_length> intron_stop:
                continue
            
            #Every read that made it to here, is an overlap read that should be counted.
            overlap_counter+=1
            #Check if left and right and adjust flag
            if points[p]=="left":
                left=True
                if read.query_name not in counted:
                    counted[read.query_name]="left"
            else:
                right=True
                if read.query_name not in counted:
                    counted[read.query_name]="right"
            
            
    
    #If either side has no overlapping reads, then the PSI score will be NAN and the spliced reads do not need to be counted.
    #to make sure we do not exclude events supported by no reads 
    if overlap_counter>0:
        if left==False or right==False:
              #if overlap_counter>0 and spliced_counter+overlap_counter>10:
                  #print("Here a numerical psi is possible if we count them as supporting.")
              #else:
                  #print("The problem is the left and right condition", left, right)
              PSI="NAN"
              return PSI
    
    #Calculate PSI
    IR=overlap_counter
    ER=spliced_counter

    if IR+ER<=10:
        #print("The problem is the read counts.", IR+ER)
        PSI="NAN"
    else:
        PSI=str(round(IR/(IR+ER),2))

    return PSI    
    

#%% Find all samples and their alignment files

#All bam files. 
argument_glob=args.samples+"/**/alignment.bam"
bam_file_list=glob.glob(argument_glob, recursive=True)

#Assuming that each bam.file is in its own sample folder:
sample_names=[i.split("/")[-4] for i in bam_file_list]

#%% Score Splicing events

#Open output file
out=open(args.out, "w")
out.write("Location\t"+"\t".join(sample_names)+"\n")

#to keep track of how many "weird" spliced reads we exclude, initiate set
excluded_reads=set()

current_gene=""
#For genes not within the input range, we have a flag
wrong_range=False
#Open input file
with open(args.input,"r") as infile:
    for line in infile:
        if line.startswith("Location"):
            #Thats the file header. Ignore
            continue
        #Find gene headers
        if line.startswith("#"):
            line=line.strip("#")
            #Check if its the first
            if current_gene=="":
                current_gene=[e.strip(" ") for e in line.strip("\n").split(",")]
                #Check if new gene in input range
                if args.coordinates:
                    if current_gene[1]!= coord_chrom:
                        wrong_range=True
                    elif int(current_gene[3])>coord_stop or int(current_gene[4])<coord_start:
                        wrong_range=True
                    else:
                        wrong_range=False
                events={"AA":{}, "AD":{}, "CE":[], "IR":[]}
                infostrings=[]
                #To save CE per sample that have PSI=1
                CE_PSI=dict()
                #To save AA spliced read counts per sample
                AA_counts=dict()
                #To save AD spliced read counts per sample
                AD_counts=dict()
                
                
            
            #We have encountered the next gene! Process all events of previous gene.
            #If we have a previous gene that is outside of range, then we only reset.
            elif wrong_range==True:
                #Reset events and current gene
                #list of gene information: gene name, chrom, strand, min, max
                current_gene=[e.strip(" ") for e in line.strip("\n").split(",")]
                #Check if new gene in input range
                if args.coordinates:
                    if current_gene[1]!= coord_chrom:
                        wrong_range=True
                    elif int(current_gene[3])>coord_stop or int(current_gene[4])<coord_start:
                        wrong_range=True
                    else:
                        wrong_range=False
                events={"AA":{}, "AD":{}, "CE":[], "IR":[]}
                infostrings=[]
                #To save CE per sample that have PSI=1
                CE_PSI=dict()
                #To save AA spliced read counts per sample
                AA_counts=dict()
                #To save AD spliced read counts per sample
                AD_counts=dict()
            
            #Previous gene is not outside of range! Horray! Make Scores
            else:
                #"Go through samples, extract reads, loop through events, do PSI"
                scores=dict()
                for sample in sample_names:
                    #initialize CE list for sample, to append CE with PSI=1
                    CE_PSI[sample]=[]
                    #Initialize AD/AA count dict per sample
                    AA_counts[sample]=dict()
                    AD_counts[sample]=dict()
                    #calculate PSI scores
                    PSI_scores=PSI_for_Gene(current_gene, events, sample)
                    scores[sample]=PSI_scores
                #Write into output file, separately for each AS
                #Can merge them later, but it results in smaller files that way.
                for i in range(0, len(infostrings)):
                    #add AS to infostring so we know what type of event it is.
                    newline=infostrings[i]
                    for sample in scores:
                        newline+="\t"+scores[sample][i]
                    out.write(newline+"\n")     
                    
                #Reset events and current gene
                #list of gene information: gene name, chrom, strand, min, max
                current_gene=[e.strip(" ") for e in line.strip("\n").split(",")]
                #Check if in range
                if args.coordinates:
                    if current_gene[1]!= coord_chrom:
                        wrong_range=True
                    elif int(current_gene[3])>coord_stop or int(current_gene[4])<coord_start:
                        wrong_range=True
                    else:
                        wrong_range=False
                events={"AA":{}, "AD":{}, "CE":[], "IR":[]}
                infostrings=[]
                #To save CE per sample that have PSI=1
                CE_PSI=dict()
                #To save AA spliced read counts per sample
                AA_counts=dict()
                #To save AD spliced read counts per sample
                AD_counts=dict()
            
        #If there is no header, it is an AS event entry
        #These are sorted by alphabet. So it is always AA, AD, CE, IR
        else:
            #If we are in the wrong range, dont save elements
            if wrong_range==True:
                continue
            #Extract all events for this gene
            current_AS=line.strip("\n").split("\t")[1].strip(" ")
            infostrings.append(current_AS+"_"+line.split("\t")[0])
            info=line.split("\t")[0].split("_")
            #Different ways of saving them for different AS. First AA.
            if current_AS=="AA":
                #check if this is an AS event of interest
                if "AA" not in inputs:
                    continue
                
                #1.column has format eventnumber, chrom, strand, start coord exon 2
                #events for AA will be eventnumber: [chrom, strand, startcoord]
                if info[0] not in events["AA"]:
                    events["AA"][info[0]]=[info[1:]]
                else:
                    events["AA"][info[0]].append(info[1:])
                
            #AD
            elif current_AS=="AD":
                #check if this is an AS event of interest
                if "AD" not in inputs:
                    continue
                
                #1.column has format eventnumber, chrom, strand, stopcoord exon 1
                #events for AD will be eventnumber: [chrom, strand, stopcoord]
                if info[0] not in events["AD"]:
                    events["AD"][info[0]]=[info[1:]]
                else:
                    events["AD"][info[0]].append(info[1:])
            
            #CE
            elif current_AS=="CE":
                #check if this is an AS event of interest
                if "CE" not in inputs:
                    continue
                
                #1. column has format: chrom, strand, exonstart, exonend
                #Save to list
                events["CE"].append(info)
            
            #IR
            elif current_AS == "IR":
                #check if this is an AS event of interest
                if "IR" not in inputs:
                    continue
                
                #1. column has format: chrom, strand, intronstart, intronend
                #Save to list
                events["IR"].append(info)

#Last genes scores still need to be written:
if wrong_range==False:
    #Previous gene is not outside of range! Horray! Make Scores
    #"Go through samples, extract reads, loop through events, do PSI"
    scores=dict()
    for sample in sample_names:
        #initialize CE list for sample, to append CE with PSI=1
        CE_PSI[sample]=[]
        #Initialize AD/AA count dict per sample
        AA_counts[sample]=dict()
        AD_counts[sample]=dict()
        #Calculate PSI scores.
        PSI_scores=PSI_for_Gene(current_gene, events, sample)
        
        scores[sample]=PSI_scores
    #Write into output file, separately for each AS
    #Can merge them later, but it results in smaller files that way.
    for i in range(0, len(infostrings)):
        #add AS to infostring so we know what type of event it is.
        newline=infostrings[i]
        for sample in scores:
            newline+="\t"+scores[sample][i]
        out.write(newline+"\n")     
        
    #Reset events and current gene not necessary.

out.close()



#%% End time

print("Run time: {:.2f} seconds.".format(time.time()-start_time))  

