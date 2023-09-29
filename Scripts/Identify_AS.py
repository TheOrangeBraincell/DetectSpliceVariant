# -*- coding: utf-8 -*-
"""
Date: Wed Jan 18 09:22:39 2023
File Name: Identify_AS.py
Author: Mirjam Karlsson-Müller

Description:
    Finds AS events in annotated genes from GENCODE and RefSeq.
    Returns Table for a gene, containing all AS event locations and types.
    Also returns .bed files if wished for, one per type event and gene, with event location entries, possible to view in f.e. igv.

Abbreviations:
    AS=     Alternative Splicing
    CE=     Casette Exon
    AA=     Alternative Acceptor
    AD=     Alternative Donor
    IR=     Intron Retention
    
List of Functions:
    add_to(): 
        Adds exon to AD and/or AA event. If the event this exon belongs to does not exist yet, creates new one.
        Specific cases for plus and negative strand.
    
    database_read():
        Processes NCBI and GENCODE input, in a function so it can be run parallel.
    
Procedure: 
    1. Creates a dictionary with transcript IDs belonging to the gene of interest and their annotated exons from databases GENCODE and RefSeq.
    2. Going through this dictionary, finds potential AS events:
        AD/AA: Find exons with overlap, put them into AA, AD dictionary according to type of overlap and strand.
        CE: Every exon thats not first or last is a potential CE.
        IR: Every space between two exons in the same transcript can be a potential IR event.
    
Useage:
    Inputs:
        - Database files for GENCODE and RefSeq respectively -g/-r
        - output file name for table of eventlocations/types -o
        - gene name -n
        - Whether bed file is wished for -b (True or False)
    
    Outputs:
        - table of event locations and types, as well as gene information and transcript information.
        - (optional) a .bed file containing the location of AS events.


    Instructions:
        Run in command line. For example.
        
        #Estrogen Receptor alpha (ESR1)
        python gitrepo/Identify_AS.py -o AS_events_ESR1.tsv -n "ESR1" -g ~/MasterProject/Database/hg38_GENCODE39_all.tsv -r ~/MasterProject/Database/hg38_NCBI_all.tsv 
        
        #BRCA1 (neg strand)
        python ../../gitrepo/Identify_AS.py -o AS_events_BRCA1.tsv -n "BRCA1" -g ~/MasterProject/Database/hg38_GENCODE39_all.tsv -r ~/MasterProject/Database/hg38_NCBI_all.tsv        
       
    
Possible Bugs:
    -not compatible with other database tables, as specifically tailored to the two used.
    Instructions on how to retrieve the right ones from the ucsc table browser can be found
    on the github.
"""


#%% 0.0 Imports

import argparse
import re
import time
from multiprocessing import Pool

#%% 0.1 argparse

"Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='Identify Alternative Splicing',
                                 usage='%(prog)s -o OUTPUT-FILE \
                                     -g GENCODE-FILE -r REFSEQ-FILE \
                                         [-c] "chrX:XXXXXX-XXXXXX" -as AS-TYPE',
                                 description="""Per AS event of interest, creates
                                 a table with the PSI scores supporting said event
                                 per sample in sample folder.""")

parser.add_argument('--out', '-o', required=True,
                    help="Output file, containing AS events for the input gene.")
parser.add_argument('--gene', '-n', required=True,
                    help="Name of the gene we are finding AS for.")
parser.add_argument('--gencode', '-g', required=True,
                    help="""tsv file containing bed file information on 
                    annotated exons from GENCODE39 as well as gene names.""")
parser.add_argument('--refseq', '-r', required=True,
                    help="""tsv file containing bed file information on 
                    annotated exons from RefSeq as well as gene names.""")
parser.add_argument('--bed', '-b', type=bool,
                    help="If an output file bed file for each type of event is wished for,\
                        set to True")


args = parser.parse_args()

#%% 0.2 Functions

def database_read(file):
    gene_chrom="NA"
    with open(file, "r") as infile:
        for line in infile:
            # To exclude potential title lines/empty lines, formatting mistakes
            # Only takes chr[] and chr[]_random lines, in accordance with bam.
            if re.search(r"(.+)\t(?:([a-z]{3}[X,M,Y]?\d*)|([a-z]{3}[X,M,Y]?\d*)"
                         r".+_random)\t(\-?\+?)\t(\d+)\t(\d+)\t(\d+)\t([\d,]+)"
                         r"\t([\d,]+)\t(.+)", line):
                # specify groups.
                entry = re.search(r"(.+)\t([a-z]{3}[X,M,Y]?\d*).*\t"
                             r"(\-?\+?)\t(\d+)\t(\d+)\t(\d+)\t([\d,]+)"
                             r"\t([\d,]+)\t(.+)", line)
                #Assign variables to groups
                trans_ID=entry.group(1)
                chrom=entry.group(2)
                strand=entry.group(3)
                gene_name=entry.group(9)
                
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
                
                """To not get caught in exon start/stop, -/+ strand complications,
                coordinates are referred to as bigger and smaller instead."""
                
                number_exons=int(entry.group(6))
                exon_smaller=entry.group(7).split(",")[0:-1]
                exon_bigger=entry.group(8).split(",")[0:-1]
                gene_name=entry.group(9)
                if file==args.gencode:
                    db="G"
                elif file==args.refseq:
                    db="R"
                else:
                    db="You got a bug concerning the db check in database_read function."
                
                #Add transcript ID to gene dictionary.
                gene_dict[trans_ID]=[]
                
                #make entries for each exon.
                for i in range(0, number_exons):
                    if i==0:
                        if strand=="+":
                            position= "first"
                        else:
                            position="last"
                    elif i==number_exons-1:
                        if strand=="+":
                            position="last"
                        else:
                            position="first"
                    else:
                        position="middle"
                    gene_dict[trans_ID].append([chrom, exon_smaller[i], 
                                                exon_bigger[i], strand, position,db])
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

def add_to(events):
    for e in events:
        #e can be AA or AD or a list of both.
        if e=="AA":
            #print(entry, coord_exons[key_string])
            if entry[4]=="first" or coord_exons[key_string][4]=="first":
                continue
            #We want to save the actual exon starts (not just smaller coordinate.)
            strand=entry[3]
            if strand=="+":
                start1=int(entry[1])
                start2=int(coord_exons[key_string][1])
            
            else:
                start1=int(entry[2])
                start2=int(coord_exons[key_string][2])
            #Go through previous entries in potential_AA to make sure theres no duplicates.
            event_exists=False
            #potential_AA: {1:[[start1, start2], {start1:{db, trans_id_dict(trans ids in set)}, start2:{db, trans_id_dict(trans ids in set)}}]}
            for events in potential_AA:
                if start1 in potential_AA[events][0] and start2 not in potential_AA[events][0]:
                    potential_AA[events][0].append(start2)
                    potential_AA[events][1][start2]=[coord_exons[key_string][5],coord_exons[key_string][6]]
                    event_exists=True
                elif start2 in potential_AA[events][0] and start1 not in potential_AA[events][0]:
                    potential_AA[events][0].append(start1)
                    potential_AA[events][1][start1]=[entry[5], entry[6]]
                    event_exists=True
                elif start1 in potential_AA[events][0] and start2 in potential_AA[events][0]:
                    #Both starts are already found in there, but that means that theres two exons with one of these starts and we need the db and trans_ID info of the new one
                    event_exists=True
                    #first update db if necessary
                    if potential_AA[events][1][start1][0]!=entry[5]:
                        #either the potential already has both db, or it has one and the new one has the other. so either way, it will be both after this, since they arent the same.
                        potential_AA[events][1][start1][0]="GR"
                    if  potential_AA[events][1][start2][0]!=coord_exons[key_string][5]:
                        #either the potential already has both db, or it has one and the new one has the other. so either way, it will be both after this, since they arent the same.
                        potential_AA[events][1][start2][0]="GR"
                    #if the databases are the same, then theres no need to update.
                    
                    #Add all transcript ids to set, duplicates are not saved anyways. start1 first
                    for letter in entry[5]:
                        if potential_AA[events][1][start1][1][letter]=="-":
                            potential_AA[events][1][start1][1][letter]=entry[6][letter]
                        else:
                            #Add transcript ids of this entry corresponding to this db into the set of the potential_AA
                            potential_AA[events][1][start1][1][letter].update(entry[6][letter])
                    #and for start2
                    for letter in coord_exons[key_string][5]:
                        if potential_AA[events][1][start2][1][letter]=="-":
                            potential_AA[events][1][start2][1][letter]=coord_exons[key_string][6][letter]
                        else:
                            potential_AA[events][1][start2][1][letter].update(coord_exons[key_string][6][letter])
            
            #If the event is not found, make a new one.
            if event_exists==False:
                new_index=len(list(potential_AA.keys()))
                potential_AA[new_index]=[[int(start1), int(start2)], {start1: [entry[5], entry[6]], start2:[coord_exons[key_string][5], coord_exons[key_string][6]]}]
                
            
        if e=="AD":
            if entry[4]=="last" or coord_exons[key_string][4]=="last":
                continue
            
            #save exon stops (actual, not just bigger coordinate)
            strand=entry[3]
            if strand=="+":
                stop1=int(entry[2])
                stop2=int(coord_exons[key_string][2])
            
            else:
                stop1=int(entry[1])
                stop2=int(coord_exons[key_string][1])
            #Go through previous entries in potential_AA to make sure theres no duplicates.
            event_exists=False
            for events in potential_AD:
                if stop1 in potential_AD[events][0] and stop2 not in potential_AD[events][0]:
                    potential_AD[events][0].append(stop2)
                    potential_AD[events][1][stop2]=[coord_exons[key_string][5],coord_exons[key_string][6]]
                    event_exists=True
                elif stop2 in potential_AD[events][0] and stop1 not in potential_AD[events][0]:
                    potential_AD[events][0].append(stop1)
                    potential_AD[events][1][stop1]=[entry[5], entry[6]]
                    event_exists=True
                elif stop1 in potential_AD[events][0] and stop2 in potential_AD[events][0]:
                    event_exists=True
                    #Both starts are already found in there, but that means that theres two exons with one of these starts and we need the db and trans_ID info of the new one
                    #first update db if necessary
                    if potential_AD[events][1][stop1][0]!=entry[5]:
                        #either the potential already has both db, or it has one and the new one has the other. so either way, it will be both after this, since they arent the same.
                        potential_AD[events][1][stop1][0]="GR"
                    if  potential_AD[events][1][stop2][0]!=coord_exons[key_string][5]:
                        #either the potential already has both db, or it has one and the new one has the other. so either way, it will be both after this, since they arent the same.
                        potential_AD[events][1][stop2][0]="GR"
                    #if the databases are the same, then theres no need to update.
                    
                    #Add all transcript ids to set, duplicates are not saved anyways. start1 first
                    for letter in entry[5]:
                        if potential_AD[events][1][stop1][1][letter]=="-":
                            #No transcript ids yet for this dataset.
                            potential_AD[events][1][stop1][1][letter]=entry[6][letter]
                        else:
                            #Add transcript ids of this entry corresponding to this db into the set of the potential_AA
                            potential_AD[events][1][stop1][1][letter].update(entry[6][letter])
                    #and for start2
                    for letter in coord_exons[key_string][5]:
                        if potential_AD[events][1][stop2][1][letter]=="-":
                            potential_AD[events][1][stop2][1][letter]=coord_exons[key_string][6][letter]
                        else:
                            potential_AD[events][1][stop2][1][letter].update(coord_exons[key_string][6][letter])
                    
            #If the event is not found, make a new one.
            if event_exists==False:
                new_index=len(list(potential_AD.keys()))
                potential_AD[new_index]=[[int(stop1), int(stop2)], {stop1: [entry[5], entry[6]], stop2:[coord_exons[key_string][5], coord_exons[key_string][6]]}]

#%% 0.3 Timer

start_time=time.time()

print("Starting Identify AS script! ")
#%% 1. Process databases in parallel.

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

#Save gene range and other information for later.
starts=[]
stops=[]
for trans_ID in gene_dict:
    starts.append(int(gene_dict[trans_ID][0][1]))
    stops.append(int(gene_dict[trans_ID][-1][2]))
    chrom=gene_dict[trans_ID][0][0]
    strand=gene_dict[trans_ID][0][3]

#Take smallest start and biggest stop as range of gene.
gene_ranges=[chrom, strand, min(starts), max(stops)]

#%% 2. Prepare output file for identified AS events
#We will write them as we identify them.

out=open(args.out,"w")
out.write("Location\tAS_Type\tGene_Name\tDatabase\tGencode_Transcript_IDs\tRefseq_Transcript_IDs\n")

#Inputs previously, now is just all types of AS we are investigating.
inputs=["AA", "AD", "CE", "IR"]

#Write gene header for output file
out.write("#"+args.gene+" ,"+str(gene_ranges[0])+" ,"+str(gene_ranges[1])+
          " ,"+str(gene_ranges[2])+" ,"+str(gene_ranges[3])+"\n")

#%% 3. Alternative Acceptors and Donors

print("Identifying Alternative Acceptors and Donors...", end="\r")

#For both AA and AD, we compare the exon coordinates in different transcripts for any overlaps.
#We first collect all exons for this gene, and then identify the potential events.
#Remove transcript ids and duplicate exons. Collect all exons in nested list.
gene_exons=[]
for trans_ID in gene_dict:
    for exon in gene_dict[trans_ID]:
        #exon is now: [chrom, small, big, strand, position, db]
        #Because db and transid info wont match necessarily, we only check for up to -2
        no_match=True
        for index,entry in enumerate(gene_exons):
            if exon[0:-1] == entry[0:-2]:
                no_match=False
                #if its already there, update db and Transcript ID

                #update db according to whether this one is included already.
                #The if statement will catch it only if its the opposite database by itself at this point in the entry.
                if exon[-1] not in entry[-2]:
                    entry[-2]=entry[-2]+exon[-1]
                    #Then the dictionary at the end of entry will also have - for the new database.
                    entry[-1][exon[-1]]={trans_ID}
                else:
                    #there should already be a transcript id from that database in the dictionary at the last spot. append.
                    entry[-1][exon[-1]].add(trans_ID)
                    
                #Update gene_exons entry
                gene_exons[index]=entry
        #If there was no match found, then we have a completely new exon and append it to gene_exons
        if no_match==True:
            if exon[-1]=="G":
                temp_list=[i for i in exon]
                temp_list.append({exon[-1]:{trans_ID}, "R":"-"})
                gene_exons.append(temp_list)
            elif exon[-1]=="R":
                temp_list=[i for i in exon]
                temp_list.append({exon[-1]:{trans_ID}, "G":"-"})
                gene_exons.append(temp_list)
            else:
                print("WTF? ",  exon)
                
#Go through all exons per gene to find alternative donors/acceptors.
#Initiate potential AA/AD for each gene.
potential_AA=dict()
potential_AD=dict()
#nested with start, stop
coordinates=[]
#key=start_stop, value=entry.
coord_exons=dict()
#Go through exons to find potentials
for entry in gene_exons: #list of [chrom, small, big, strand, position, db, {G:[], R:[]}]
    start=int(entry[1])
    stop=int(entry[2])
    strand=entry[3]
    
    """If the start or stop is shared with another exon, then we found potential AS.
    But one of the exons can also be completely within another exon, or completely
    surrounding the other one. Or having one coordinate within the other. That also counts. So we find any type of overlap"""
    
    for coordinate in coordinates:
        key_string=str(coordinate[0])+"_"+str(coordinate[1])
        #if the start coordinate lies within a previous exon.
        if coordinate[0]<=start<coordinate[1] and coordinate[1]!=stop:
            #Then either they share start, and we only have one event
            if start==coordinate[0]:
                # if + AD, if - AA.
                if strand=="+":
                    add_to(["AD"])
                else:
                    add_to(["AA"])
            #or we have both events
            else:
                #both
                add_to(["AA","AD"])
                
        #if the stop coordinate lies within a previous exon.
        elif coordinate[0]+1< stop <= coordinate[1] and start!= coordinate[0]:
            #either they share stop, and we have one event
            if stop==coordinate[1]:
                #if - AD, if + AA.
                if strand=="-":
                    add_to(["AD"])
                else:
                    add_to(["AA"])
            #Or we have both events.
            else:
                #both
                add_to(["AA","AD"])
        # if this exon contains a previous exon completely then we have both
        elif coordinate[0] > start and coordinate[1] < stop:
            #both
            add_to(["AA","AD"])
    
    #Add new coordinates to compare new entries to.
    coord_exons[str(start)+"_"+str(stop)]=entry
    coordinates.append([start, stop])
#print(potential_AD, "\n")
print("Identifying Alternative Acceptors and Donors: Done! \n", end="\r")

#%% 3.1 Alternative Acceptors

print("Writing Alternative Acceptors to Output ...", end="\r")

to_sort_by=dict()
#Sort starts for each aa_event.
for aa_event in potential_AA:
    potential_AA[aa_event][0].sort()
    to_sort_by[potential_AA[aa_event][0][0]]=aa_event
#we want the AA events sorted too, which might mess with the indeces given to the events prior.
new_order=sorted(list(to_sort_by.keys()))
sorted_AA=dict()
for index, item in enumerate(new_order):
    sorted_AA[index+1]=potential_AA[to_sort_by[item]]

#write sorted events
for aa_event in sorted_AA:
    for starts in sorted_AA[aa_event][0]:
        #Use entry to write output file.
        out.write("{}_{}_{}_{}\t{}\t{}\t{}\t{}\t{}\n".format(aa_event, gene_ranges[0], gene_ranges[1], starts, "AA", args.gene, 
                                                             sorted_AA[aa_event][1][starts][0], ",".join(sorted_AA[aa_event][1][starts][1]["G"]), 
                                                             ",".join(sorted_AA[aa_event][1][starts][1]["R"])))
#If a bed output file is wished for, it is created here.
if args.bed==True:
    bed=open(args.out+"AA.bed", "w")
    if len(potential_AA)!=0:
        #Extract strand information for header (for each gene)
        strand=gene_ranges[1]
        chrom=gene_ranges[0]
    
    for entry in potential_AA:
        smallest_start=min(entry)
        biggest_start=max(entry)
        name="AA_"+chrom+"_"+str(smallest_start)+"_"+str(biggest_start)
        bed.write("{}\t{}\t{}\t{}\t.\t{}\n".format(chrom, str(smallest_start), str(biggest_start), name, strand))
    bed.close()

print("Writing Alternative Acceptors to Output: Done! \n", end="\r")

#%% 3.2 Alternative Donors.

print("Writing Alternative Donors to Output ...", end="\r")

to_sort_by=dict()
#Sort entries after location and coordinates
for ad_event in potential_AD:
    potential_AD[ad_event][0].sort()
    to_sort_by[potential_AD[ad_event][0][0]]=ad_event

#we want the AA events sorted too, which might mess with the indeces given to the events prior.
new_order=sorted(list(to_sort_by.keys()))
sorted_AD=dict()
for index, item in enumerate(new_order):
    sorted_AD[index+1]=potential_AD[to_sort_by[item]]

#write sorted events
for ad_event in sorted_AD:
    for stops in sorted_AD[ad_event][0]:
        #Use entry to write output file.
        out.write("{}_{}_{}_{}\t{}\t{}\t{}\t{}\t{}\n".format(ad_event, gene_ranges[0], gene_ranges[1], stops, "AD", args.gene, 
                                                             sorted_AD[ad_event][1][stops][0], ",".join(sorted_AD[ad_event][1][stops][1]["G"]), 
                                                             ",".join(sorted_AD[ad_event][1][stops][1]["R"])))
#If a bed file is wished for, it is created here.
if args.bed==True:
    bed=open(args.out+"AD.bed", "w")
    if len(potential_AD)!=0:
        #Extract strand information for header (for each gene)
        strand=gene_ranges[1]
        chrom=gene_ranges[0]
    for entry in potential_AD:
        smallest_stop=min(entry)
        biggest_stop=max(entry)
        name="AA_"+chrom+"_"+str(smallest_stop)+"_"+str(biggest_stop)
        bed.write("{}\t{}\t{}\t{}\t.\t{}\n".format(chrom, str(smallest_stop), str(biggest_stop), name, strand))
    bed.close()  
    
print("Writing Alternative Donors to Output: Done! \n", end="\r")

#%% 4. Casette Exons.

print("Identifying Casette Exons ...", end="\r")
#Remove first and last exon for each transcript, as they cannot be CE
CE_exons=[]
for trans_ID in gene_dict:
    #If a transcript has 2 or less exons, there is no CE
    if len(gene_dict[trans_ID])<=2:
        continue
    for exon in gene_dict[trans_ID]:
        if exon[4]=="middle":
            #Check if they are already in CE exons
            no_match=True
            for index,entry in enumerate(CE_exons):
                if exon[0:-1] == entry[0:-2]:
                    no_match=False
                    #if its already there, update db and Transcript ID

                    #update db according to whether this one is included already.
                    #The if statement will catch it only if its the opposite database by itself at this point in the entry.
                    if exon[-1] not in entry[-2]:
                        entry[-2]=entry[-2]+exon[-1]
                        #Then the dictionary at the end of entry will also have - for the new database.
                        entry[-1][exon[-1]]=[trans_ID]
                    else:
                        #there should already be a transcript id from that database in the dictionary at the last spot. append.
                        entry[-1][exon[-1]].append(trans_ID)
                    #Update gene_exons entry
                    CE_exons[index]=entry
            #If there was no match found, then we have a completely new exon and append it to gene_exons
            if no_match==True:
                if exon[-1]=="G":
                    temp_list=[i for i in exon]
                    temp_list.append({exon[-1]:[trans_ID], "R":"-"})
                    CE_exons.append(temp_list)
                elif exon[-1]=="R":
                    temp_list=[i for i in exon]
                    temp_list.append({exon[-1]:[trans_ID], "G":"-"})
                    CE_exons.append(temp_list)
                else:
                    print("WTF? ",  exon)
            
#Sort exons after start coordinate
CE_exons.sort(key=lambda x: int(x[1]))
#Write potential CE into output for this gene.
for exon in CE_exons:
    out.write("{}_{}_{}_{}\t{}\t{}\t{}\t{}\t{}\n".format(exon[0], exon[3],
                                  exon[1], exon[2], "CE", args.gene, exon[-2], ",".join(exon[-1]["G"]), ",".join(exon[-1]["R"])))
#If bed file wanted
if args.bed==True:
    bed=open(args.out+"CE.bed", "w")
    for exon in CE_exons:
        chrom, start, stop, strand = entry[0:4]
        name="CE_"+chrom+"_"+start+"_"+stop
        bed.write("{}\t{}\t{}\t{}\t.\t{}\n".format(chrom, start, stop, name, strand))
    bed.close()
    
print("Identifying Casette Exons: Done! \n", end="\r")

#%% 5. Intron Retention
#Retention of introns is usually not annotated. So we do the same as for CE.
#And check for every gap between exons, if there is read showing intron retention.
#That means at this point we save all potential gaps.

print("Identifying Intron Retention... ", end="\r")

#Remove transcript ids and duplicate exons.
IR_coord=[]
for trans_ID in gene_dict:
    for i in range(0, len(gene_dict[trans_ID])-1):
        chrom=gene_dict[trans_ID][i][0]
        strand=gene_dict[trans_ID][i][3]
        db=gene_dict[trans_ID][i][-1]
        IR=[chrom,strand,gene_dict[trans_ID][i][2],gene_dict[trans_ID][i+1][1], db, trans_ID]
        #Check if already in list
        no_match=True
        for entry in IR_coord:
            if IR[0:-2]== entry[0:-2]:
                #theres a match
                no_match=False
                #if its already there, update db and Transcript ID
                #if the db is not in the entries db, then it needs to be added to db and dict.
                if IR[-2] not in entry[-2]:
                    #add the db strings
                    entry[-2]=entry[-2]+IR[-2]
                    #replace "-" in dictionary, initialize transcript list for new db
                    entry[-1][IR[-2]]=[IR[-1]]
                else:
                    #just add the trans ID to the dict, as db is already ok.
                    entry[-1][IR[-2]].append(IR[-1])

        #If there was no match found, then we have a completely new exon and append it to gene_exons
        if no_match==True:
            temp_list=IR[0:-1]
            if IR[-2]=="G":
                temp_list.append({"G":[IR[-1]], "R":"-"})
                IR_coord.append(temp_list)
            elif IR[-2]=="R":
                temp_list.append({"R":[IR[-1]], "G":"-"})
                IR_coord.append(temp_list)
            else:
                print("WTF ", IR)
            
#Sort IR entries by start of intron
IR_coord.sort(key=lambda x: x[2])
#Write IR entries into output file
for IR in IR_coord:
    out.write("{}_{}_{}_{}\t{}\t{}\t{}\t{}\t{}\n".format(IR[0], IR[1], IR[2], IR[3], "IR", args.gene, IR[-2], ",".join(IR[-1]["G"]), ",".join(IR[-1]["R"])))

#If bed file wished for
if args.bed==True:
    bed=open(args.out+"IR.bed", "w")
    for gene in IR_coord:
        strand=IR_coord[0].split("_")[0]
        
        for entry in IR_coord:
            chrom = entry.split("_")[1]
            start=entry.split("_")[2]
            stop=entry.split("_")[3]
            name="IR_"+chrom+"_"+start+"_"+stop
            bed.write("{}\t{}\t{}\t{}\t.\t{}\n".format(chrom, start, stop, name, strand))
    bed.close()
print("Identifying Intron Retention: Done! \n", end="\r")

#%% 6. Close output file
out.close()

#%% Timer

print("Run time: {:.2f} seconds.".format(time.time()-start_time)) 
