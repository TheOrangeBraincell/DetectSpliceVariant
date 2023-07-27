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
                                         "IR":{"count":0, "reads":{"spliced":[], "not":[]}, "count_n":0, 
                                               "normalize":int(event.split("_")[4])-int(event.split("_")[3])}}
        elif event.startswith("AA"):
            #dict: {AD_#: {start1: dict, start2: dict}}
            starts=sorted(list(AS_events[event].keys()))
            #The order of the starts is reversed if we are on the minus strand! 
            if gene_strand=="+":
                starts=starts[::-1]
            for start in AS_events[event]:
                #Calculate normalize as difference between previous start and current.
                current_index=starts.index(start)
                if current_index==0:
                    #Its the start furtherst of the intron, it wont have any difference reads. hence normalization=1
                    normalize=1
                else:
                    previous_index=current_index-1
                    current_start=int(start.split("_")[2])
                    previous_start=int(starts[previous_index].split("_")[2])
                    normalize=max([current_start, previous_start])-min([current_start, previous_start])
                
                AS_events[event][start]={"IR": {"count":0,"count_n":0, "normalize":normalize}}
            #Reads should not be counted several for different starts, so we keep track of it outside the counting dict.
            counted[event]={"spliced":[], "not":[]}
            
        elif event.startswith("AD"):
            #dict: {AD_#: {stop1: dict, stop2: dict}}
            stops=sorted(list(AS_events[event].keys()))
            #on the plus strand the innermost acceptor has the lowest coordinate, so we reverse the order
            if gene_strand=="-":
                stops=stops[::-1]
            for stop in AS_events[event]:
                #Calculate normalize as difference between previous stop and current.
                current_index=stops.index(stop)
                if current_index==0:
                    #This is the outermost stop, so there is no difference to calculate.
                    normalize=1
                else:
                    previous_index=current_index-1
                    current_stop=int(stop.split("_")[2])
                    previous_stop=int(stops[previous_index].split("_")[2])
                    normalize=max([current_stop, previous_stop])-min([current_stop, previous_stop])
                
                AS_events[event][stop]={"IR": {"count":0,"count_n":0, "normalize":normalize}}
            #Reads should not be counted several for different stops, so we keep track of it outside the counting dict.
            counted[event]={"spliced":[], "not":[]}
        
        else:
            #IR event.
            #dict: {}
            AS_events[event]={"ER":0,  #spliced across
                              "IR": {}, #overlap reads
                              "reads":{"spliced":[], "left":[], "right":[]}}
            #Note that to add all points to IR, ill have to check them for overlaps with AA/AD events
            #But the intron coordinates of this we can definitely put in.
            #we assign left and right based on bigger and smaller, not based on biology.
            AS_events[event]["IR"][event.split("_")[3]]={"count":0, "side":"left"}
            AS_events[event]["IR"][event.split("_")[4]]={"count":0, "side":"right"}
            
            #Check for overlapping AA, AD events.
            for e in AS_events:
                if e.startswith("A"):
                    #check if theres a match:
                    if e.split("_")[4]== event.split("_")[3]:
                        #then the point is on the left. there will be more left points.
                        extra_points=[i for i in list(AS_events.keys()) if i.split("_")[1].startswith(e.split("_")[1])]
                        #These go into dictionary with left side
                        for point in extra_points:
                            AS_events[event]["IR"][point]={"count":0, "side":"left"}
                        break
                    elif e.split("_")[4]==event.split("_")[4]:
                        #then the point is on the right. there will be more right points.
                        extra_points=[i for i in list(AS_events.keys()) if i.split("_")[1].startswith(e.split("_")[1])]
                        #These go into dictionary with right side
                        for point in extra_points:
                            AS_events[event]["IR"][point]={"count":0, "side":"right"}
                        break
            
    #open the bam file, and go through the reads.
    samfile=pysam.AlignmentFile(files[sample], 'rb', index_filename=files[sample][0:-1]+"i")
    reads=samfile.fetch(chrom, gene_start, gene_stop)
    
    #Iterate through reads, count them for events.
    for read in reads:
        #Filter
        if Filter_Reads(read, gene_strand)==True:
            continue
        #Now we go through all the events and try to figure out what to count them for.
        for event in AS_events:
            #event=event id later used in first column of output.
            #Extract CE info from event key.
            if event.startswith("CE"):
                skip=CE(sample, event, read)
                if skip==True:
                    #The read is invalid all together, so we want to break the inner loop.
                    break
                #Else the read just doesnt match the current event and we let the program continue.
            elif event.startswith("AA"):
                skip=AA(sample, event, read)
                if skip==True:
                    #The read is invalid all together and does not need to be investigated further.
                    break
            elif event.startswith("AD"):
                skip=AD(sample, event, read)
                if skip==True:
                    #The read is invalid all together and does not need to be investigated further.
                    break
                
            elif event.startswith("IR"):
                continue  
    #Now we have counted all the reads, we can calculate the PSI scores.
    sample_PSI={"sample":sample}
    for event in AS_events:
        if event.startswith("CE"):
            #PSI is NAN if there are 10 or less reads.
            if AS_events[event]["IR"]["count"]+AS_events[event]["IR"]["count_n"]+AS_events[event]["ER"]["count"]<11:
                PSI="NAN"
            else:
                #sum IR, normalize count_n
                IR=AS_events[event]["IR"]["count"]+(AS_events[event]["IR"]["count_n"]/AS_events[event]["IR"]["normalize"])
                #assign ER
                ER=AS_events[event]["ER"]["count"]
                PSI=str(round(IR/(IR+ER), 3))
            sample_PSI[event]=PSI
            
        elif event.startswith("A"):
            #several PSI scores for one AA /AD event
            #PSI is NAN if there are 10 or less reads.
            total_reads=0
            for coord in AS_events[event]:
                total_reads+= AS_events[event][coord]["IR"]["count"]+AS_events[event][coord]["IR"]["count_n"]
            
            if total_reads <11:
                PSI="NAN"
                for coord in AS_events[event]:
                    sample_PSI[event+"_"+coord]=PSI 
            else:
                #The denominator of the PSI fraction is the same for all of the coords in the same ad event. so lets do it first.
                sum_ER_IR=0
                for coord in AS_events[event]:
                    sum_ER_IR+= AS_events[event][coord]["IR"]["count"]+AS_events[event][coord]["IR"]["count_n"]/AS_events[event][coord]["IR"]["normalize"]
                
                
                #Now that we have the sum of IRs, we can calculate the PSI for each IR
                for coord in AS_events[event]:
                    IR=AS_events[event][coord]["IR"]["count"]+AS_events[event][coord]["IR"]["count_n"]/AS_events[event][coord]["IR"]["normalize"]
                    PSI=str(round(IR/(sum_ER_IR), 3))
                    sample_PSI[event+"_"+coord]=PSI 
                
        elif event.startswith("IR"):
            
            PSI="NAN"
             

    return sample_PSI

def CE(sample, event, read):
    CE_start=int(event.split("_")[3])
    CE_stop=int(event.split("_")[4])
    
    #read information
    read_start=int(read.reference_start)
    read_range=sum([int(i) for i in re.findall(r'\d+', read.cigarstring)])
    read_stop=read_start+read_range
    
    #If both start and stop of read are before or after CE, then we can go on.
    if read_start < CE_start and read_stop < CE_start:
        return False
    if read_start > CE_stop and read_stop > CE_stop:
        return False
    
    #We are in the right range, lets process the read.
    current_cigar = read.cigarstring
    read_name=read.query_name
    #Spliced or not?
    if re.search(r'\d+M(?:\d+[I,D,S,H])?\d+N(?:\d+[I,D,S,H])?\d+M',read.cigarstring):
        #This is spliced, but its could have insertions and/or deletions around the intron... 
        if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
            #we save excluded weird reads like this one for the logfile.
            excluded_reads.add(read_name)
            #Then this read doesnt need to be processed at all.
            return True
        
        current_start = read_start
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
                return False
            
            #Number of spliced reads from previous to CE.
            if exon1_end< CE_start and CE_stop >= exon2_end>=CE_start:
                #count!
                if read_name not in AS_events[event]["IR"]["reads"]["spliced"]:
                    AS_events[event]["IR"]["count"]+=1
                    AS_events[event]["IR"]["reads"]["spliced"].append(read_name)
                    if read_name in AS_events[event]["IR"]["reads"]["not"]:
                        #That means we have counted this reads pair as exon read. 
                        #But we prioriize spliced reads in the counts.
                        AS_events[event]["IR"]["count_n"]-=1
                        AS_events[event]["IR"]["reads"]["not"].remove(read_name)
                
            #Number of spliced reads from CE to next.   
            elif CE_start <= exon1_start<=CE_stop and exon2_start >CE_stop:
                #count!
                if read_name not in AS_events[event]["IR"]["reads"]["spliced"]:
                    AS_events[event]["IR"]["count"]+=1
                    AS_events[event]["IR"]["reads"]["spliced"].append(read_name)
                    if read_name in AS_events[event]["IR"]["reads"]["not"]:
                        #That means we have counted this reads pair as exon read. 
                        #But we prioriize spliced reads in the counts.
                        AS_events[event]["IR"]["count_n"]-=1
                        AS_events[event]["IR"]["reads"]["not"].remove(read_name)
                    
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
        if read_name not in AS_events[event]["IR"]["reads"]["spliced"] and read_name not in AS_events[event]["IR"]["reads"]["not"]:
            AS_events[event]["IR"]["count_n"]+=1
            AS_events[event]["IR"]["reads"]["not"].append(read_name)
    return False
    
def AD(sample, event, read):
    #event= AD_#
    #AS_events[event]={stop1:{IR:{count:0, count_n:0, normalize=##}},
    #                 stop2: etc }
    
    stops=sorted(list(AS_events[event].keys()))
    gene_strand=gene[2]
    #read information
    read_start=int(read.reference_start)
    read_range=sum([int(i) for i in re.findall(r'\d+', read.cigarstring)])
    read_stop=read_start+read_range

    #If both start and stop of the read are before the smallest stop in AD, no counting
    if read_start < int(stops[0].split("_")[2]) and read_stop < int(stops[0].split("_")[2]):
        #Read is not invalid, it just doesnt match event, hence we return false
        return False
    #if both start and stop of read are after the biggest stop, no counting
    if read_start > int(stops[-1].split("_")[2]) and read_stop > int(stops[-1].split("_")[2]):
        return False

    #We are in event range, lets process the read.
    #Is it spliced or not
    read_name=read.query_name
    if re.search(r'\d+M(?:\d+[I,D,S,H])?\d+N(?:\d+[I,D,S,H])?\d+M',read.cigarstring):
        #spliced
        #This is spliced, but its could have insertions and/or deletions around the intron... 
        if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
            #This is invalid. we dont process it. So the function returns True.
            excluded_reads.add(read.query_name)
            return True
        #Now process spliced read.
        #Allow for several splice junctions in one read.
        current_cigar = read.cigarstring
        current_start = int(read.reference_start)
        
        
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
            
            #Update counters based on matching stops.
            match=False
            for stop in stops:
                #Then the AD stop have to match the exon1_stop
                if gene_strand=="+":
                    #the coordinates on plus strand are shifted by one, as the exon end is exclusive.
                    if str(stop.split("_")[2])==str(exon1_end-1):
                        match=True
                        break
                    
                else:
                    #Then AD stops have to match exon2_starts
                    if str(stop.split("_")[2])==str(exon2_start):
                        match=True
                        break
            #Update counters for matching stop.
            #The read counts as ER for every stop except the matched.
            if match==True and read_name not in counted[event]["spliced"]:
                counted[event]["spliced"].append(read_name)
                #This counts as an ER count for every stop except the one it matches.
                #for that one its an IR count.
                AS_events[event][stop]["IR"]["count"]+=1

                #The read (or its mate) could already be counted as a difference read. 
                #But we prioritize spliced. So then we need to remove the difference counts.
                if read_name in counted[event]["not"]:
                    #remove
                    counted[event]["not"].remove(read_name)
                    #Deduct counts
                    #IR for this stop
                    AS_events[event][stop]["IR"]["count_n"]-=1

            # update cigar string
            current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
            current_start= exon2_start
    else:
        #not spliced reads.
        #We already excluded all reads outside the range of the event
        #So we just need to figure out which stop they belong to.
        
        #The read stop defined above includes non matching regions, making the end coordinate not the
        #coordinate in the reference genome. So we only want M parts.
        match_length=sum([int(i.strip("M")) for i in re.findall(r'\d+M', read.cigarstring)])
        read_stop=read_start+match_length
        for i in range(1, len(stops)):
            stop1=int(stops[i-1].split("_")[2])
            stop2=int(stops[i].split("_")[2])
            #Does it overlap with this region, but not get into the next
            if gene_strand=="+":
                if stop2 >= read_stop > stop1:
                    if read_name not in counted[event]["not"]:
                        #match
                        match=stops[i]
                        AS_events[event][match]["IR"]["count_n"]+=1
                        counted[event]["not"].append(read_name)
                    break
            else:
                if stop1 >= read_start > stop2:
                    if read_name not in counted[event]["not"]:
                        #match
                        match=stops[i-1]
                        AS_events[event][match]["IR"]["count_n"]+=1
                        counted[event]["not"].append(read_name)
                    break

    return False

def AA(sample, event, read):
    #event=AA_#
    #AS_events[event]={start1:{IR:{count:0, count_n:0, normalize=##}},
    #                 start2: etc }
    
    starts=sorted(list(AS_events[event].keys()))
    gene_strand=gene[2]
    #read information
    read_start=int(read.reference_start)
    read_range=sum([int(i) for i in re.findall(r'\d+', read.cigarstring)])
    read_stop=read_start+read_range

    #If both start and stop of the read are before the smallest stop in AD, no counting
    if read_start < int(starts[0].split("_")[2]) and read_stop < int(starts[0].split("_")[2]):
        #Read is not invalid, it just doesnt match event, hence we return false
        return False
    #if both start and stop of read are after the biggest stop, no counting
    if read_start > int(starts[-1].split("_")[2]) and read_stop > int(starts[-1].split("_")[2]):
        return False
    
    #We are in event range, lets process the read.
    #Is it spliced or not
    read_name=read.query_name
    if re.search(r'\d+M(?:\d+[I,D,S,H])?\d+N(?:\d+[I,D,S,H])?\d+M',read.cigarstring):
        #spliced
        #This is spliced, but its could have insertions and/or deletions around the intron... 
        if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
            #This is invalid. we dont process it. So the function returns True.
            excluded_reads.add(read.query_name)
            return True
        #Now process spliced read.
        #Allow for several splice junctions in one read.
        current_cigar = read.cigarstring
        current_start = int(read.reference_start)
        
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
            
            #Update counters based on matching stops.
            match=False
            for start in starts:
                #Then the AA start has to match the exon2 start.
                if gene_strand=="+":
                    if str(start.split("_")[2])==str(exon2_start):
                        match=True
                        break
                    
                else:
                    #the coordinates on plus strand are shifted by one, as the exon end is exclusive.
                    #Then AD start has to match the exon 1 end.
                    if str(start.split("_")[2])==str(exon1_end-1):
                        match=True
                        break
            #Update counters for matching stop.
            #The read counts as ER for every stop except the matched.
            if match==True and read_name not in counted[event]["spliced"]:
                counted[event]["spliced"].append(read_name)
                #This counts as an ER count for every stop except the one it matches.
                #for that one its an IR count.
                AS_events[event][start]["IR"]["count"]+=1

                #The read (or its mate) could already be counted as a difference read. 
                #But we prioritize spliced. So then we need to remove the difference counts.
                if read_name in counted[event]["not"]:
                    #remove
                    counted[event]["not"].remove(read_name)
                    #Deduct counts
                    #IR for this stop
                    AS_events[event][start]["IR"]["count_n"]-=1

            # update cigar string
            current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
            current_start= exon2_start
    
    else:
        #Unspliced reads. We already excluded all reads outside the reange of the event,
        #so we just need to figure out in which difference they are in.
        
        #The read stop defined further up includes non matching regions (S, D, I), making
        #the end coordinate not the right coordinate in the reference genome. So we only want M parts.
        match_length=sum([int(i.strip("M")) for i in re.findall(r'\d+M', read.cigarstring)])
        read_stop=read_start+match_length
        for i in range(1, len(starts)):
            start1=int(starts[i-1].split("_")[2])
            start2=int(starts[i].split("_")[2])
            #Does it overlap with this region, but not get into the next
            if gene_strand=="-":
                if start2 >= read_stop > start1:
                    if read_name not in counted[event]["not"]:
                        #match
                        match=starts[i-1]
                        AS_events[event][match]["IR"]["count_n"]+=1
                        counted[event]["not"].append(read_name)
                    break
            else:
                if start1 >= read_start > start2:
                    if read_name not in counted[event]["not"]:
                        #match
                        match=starts[i]
                        AS_events[event][match]["IR"]["count_n"]+=1
                        counted[event]["not"].append(read_name)
                    break

    return False

def IR(sample, event, read):
    IR_start=int(event.split("_")[3])
    IR_stop=int(event.split("_")[4])
    
    #read_information
    read_start=int(read.reference_start)
    read_range=sum([int(i) for i in re.findall(r'\d+', read.cigarstring)])
    read_stop=read_start+read_range
    
    #If both start and stop of read are before or after IR, then we can go on.
    if read_start < IR_start and read_stop < IR_start:
        return False
    if read_start > IR_stop and read_stop > IR_stop:
        return False
    
    #we are in range. Lets process.
    current_cigar = read.cigarstring
    read_name=read.query_name
    #Spliced or not?
    if re.search(r'\d+M(?:\d+[I,D,S,H])?\d+N(?:\d+[I,D,S,H])?\d+M',read.cigarstring):
        #This is spliced, but its could have insertions and/or deletions around the intron... 
        if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
            #we save excluded weird reads like this one for the logfile.
            excluded_reads.add(read_name)
            #Then this read doesnt need to be processed at all.
            return True
        
        #spliced. then it can only be read across the intron for it to be interesting for us.
        current_start = read_start
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
                return False
            
            #Number of spliced reads from previous to CE.
            if exon1_start< IR_start and exon2_end > IR_stop:
                #count!
                if read_name not in AS_events[event]["reads"]["spliced"]:
                    AS_events[event]["ER"]+=1
                    AS_events[event]["reads"]["spliced"].append(read_name)
                    #Reads cannot be spliced across and havea pair within the intron,
                    #so we dont need to check for that.
                    
            # update cigar string
            current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
            current_start= exon2_start
    else:
        #Not spliced reads
        #We need to check for every point in IR if the matching part of the 
        #read overlaps with it.
        current_start=read_start
        #Find out coordinates of M part of read.
        while re.search(r'(\d+)M', current_cigar):
            #could have several matching regions.
            match=re.search(r'(\d+)M', current_cigar)
            match_length=int(match.group(1))
            #assign variables.
            before=current_cigar.split("M")[0]
            #D are positions on ref that are absent in query
            #I are positions on query that are absent in ref
            #M are match
            #S are soft clipped, i.e. not in reference
            #N should not be here.
            #Since we are wanting to match to reference coordinates, we add length of part before M together accordingly
            possible=["I", "D", "S"]
            current=""
            before_length=0
            for i in before:
                if i in possible:
                    if i =="D":
                        before_length+=int(current)
                    current=""
                current+=i
            
            #match coordinate can be calculated
            match_start=current_start+before_length
            match_stop=match_start+match_length
            #check overlap points.
            
            """YOU ARE HERE"""
                
            
        
        
    return False
#%% 0.3 Start Timer

start_time=time.time()

#%% 1. Read in AS events

#We will save the identifiers in this dict and then add the regions that should be counted.
AS_events=dict()
#the counts for AA and AD are saved in a dictionary counted as per group of stops/starts
#there should no read be counted double.
counted=dict()

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
            print(gene)
            continue
        
        #All thats left now is AS entries.
        Location, AS_Type, Gene_name, Database, Gencode_ID, Refseq_ID= line.strip().split("\t")
        #Initialize dictionary, different for AA/AD then for IR/CE
        if AS_Type.startswith("A"):
            if AS_Type+"_"+Location.split("_")[0] not in AS_events:
                AS_events[AS_Type+"_"+Location.split("_")[0]]={"_".join(Location.split("_")[1:]): dict()}
            else:
                AS_events[AS_Type+"_"+Location.split("_")[0]]["_".join(Location.split("_")[1:])]=dict()
            counted[AS_Type+"_"+Location.split("_")[0]]={"spliced":[], "not":[]}
        else:
            AS_events[AS_Type+"_"+Location]=dict()

#%% 2. Read in list of bam files. 

print("Reading in BAM files...", end="\r")

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

#%% 3. Read counts, per sample, allow parallelization, write output.

excluded_reads=set()

with mp.Pool(3) as pool:
    #result is a list. i.e. two gene_dicts.
    result=pool.map(PSI_for_Sample, sample_names)
# result=[]
# for sample in sample_names:
#     result.append(PSI_for_Sample(sample))

print("results are finished")

with open(args.out, "w") as outfile:
    #write header
    outfile.write("Event\t"+"\t".join(sample_names)+"\n")
    all_events=list(result[0].keys())
    #Theres an identifier in there to signal which sample, we dont have a psi score for that one.
    all_events.remove("sample")
    for event in all_events:
        new_line=event
        for sample in sample_names:
            #find right scores in PSI_dict
            for dictionary in result:
                if dictionary["sample"]==sample:
                    new_line+="\t"+ dictionary[event]
        outfile.write(new_line+"\n")
    

#%% End time

print("Run time: {:.2f} seconds.".format(time.time()-start_time))  


