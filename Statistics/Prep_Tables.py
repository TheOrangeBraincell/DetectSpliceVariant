# 29.10.2024
# Mirjam Müller
# Prep_Tables.py

# Prepares the genotype and psi tables by filtering out the rows that have less than 100 valid entries.
# I also do not need the entries that are not in any pairs of the pair table.
# And lastly I want to remove all entries in the genotype file that have only one genotype i.e. one unique value across the row after the first column, post removing NE and ND.
# write out one file with the remaining pairs and one each for the genotype and psi table.

# Runs for one gene with:
# python Prep_Tables.py Test_Set/ESR1_exgeno.tsv Test_Set/ESR1_PSI.tsv Test_Set/ESR1_pairs.tsv ESR1_filt.tsv
#
# Runs for all genes with:
# cat test_genes.txt | while read gene; do python Prep_Tables.py Test_Set/${gene}_exgeno.tsv Test_Set/${gene}_PSI.tsv Test_Set/${gene}_pairs.tsv ${gene}_filt.tsv; done

"""
import pandas as pd
import sys

# Read in files
genotype_table = pd.read_csv(sys.argv[1], sep="\t")
psi_table = pd.read_csv(sys.argv[2], sep="\t")

# Filter out rows that are not in the pair table
# Read in pair table
pair_table = pd.read_csv(sys.argv[3], sep="\t")

#The variants column has several variants in a comma separated list. I want to split them into separate rows
pair_table2 = pair_table.assign(variants=pair_table["variants"].str.split(",")).explode("variants")
pair_table=pair_table2

# remove all entries from genotype table and psi table that are not part of a pair
genotype_table3= genotype_table[genotype_table["Location"].isin(pair_table["variants"])]

psi_table1 = psi_table[psi_table["Event"].isin(pair_table["AS_event"])]

# Filter out rows with less than 100 valid entries
#Valid entries for PSI scores are all values that are not NAN
#Valid entries for Genotype scores are all values that are not ND or NE
genotype_table= genotype_table3[genotype_table3.apply(lambda x: x[1:].isin(["ND", "NE"]).sum() < 100, axis=1)]
#genotype_table = genotype_table3.loc[genotype_table3.iloc[:, 1:].isin(["ND", "NE"]).sum(axis=1) < 100] #start here

#psi_table2 = psi_table1[psi_table1.apply(lambda x: x[1:].isin(["NAN"]).sum() < 100, axis=1)]
psi_table2 = psi_table1.loc[psi_table1.iloc[:, 1:].isna().sum(axis=1) < 100]
	
# If the scores in the psi table are all the same across the row, remove the row
psi_table3 = psi_table2[psi_table2.apply(lambda x: x[1:].nunique() > 1, axis=1)]

genotype_table3 = genotype_table[genotype_table.apply(lambda x: x[1:].nunique() > 1, axis=1)]




# Update pairs table to only include Locations that are still in th genotype table and Events that are still in the psi table
pair_table2= pair_table[pair_table["variants"].isin(genotype_table["Location"])]
pair_table = pair_table2
pair_table2 = pair_table[pair_table["AS_event"].isin(psi_table["Event"])]
pair_table = pair_table2

# Write out files
genotype_table3.to_csv(sys.argv[4], sep="\t", index=False)
psi_table3.to_csv(sys.argv[5], sep="\t", index=False)
pair_table.to_csv(sys.argv[6], sep="\t", index=False)

"""

#Lets do an alternative way without pandas cause shit aint working.

import sys


paired_events=[]
paired_variants=[]
all_pairs=[]

with open(sys.argv[3], "r") as pairs:
    for line in pairs:
        if line.startswith("AS_event"):
            pairs_header=line
            continue
        
        if len(line.strip().split("\t"))==2:
            event = line.split("\t")[0]
            variants=line.strip().split("\t")[1].split(",")

        else:
            continue
            
        paired_events.append(event)
        for v in variants:
            paired_variants.append(v)
            all_pairs.append([event, v])
        

genotype_dict=dict()

with open(sys.argv[1], "r") as genotype:
    for line in genotype:
        if line.startswith("Location"):
            genotype_header=line
            samples=line.strip().split("\t")[1:]
            continue
        
        id = line.strip().split("\t")[0]
        entries=line.strip().split("\t")[1:]
        """
        #count how many are not ND or NE
        valid_count=0
        valid_list=[]
        for e in entries:
            if e not in ["NE", "ND"]:
                valid_count+=1
                if e not in valid_list:
                    valid_list.append(e)
        
        #see if the line has more than 100 valid entries
        if valid_count<100:
            continue

        #If there is only one genotype for this location, no test possible, skip.
        if len(valid_list)==1:
            continue
        
        if id not in paired_variants:
            continue
        """
        genotype_dict[id]=dict()
        for i in range(len(entries)):
            genotype_dict[id][samples[i]]=entries[i]


psi_dict=dict()
with open(sys.argv[2], "r") as psi:
    for line in psi:
        if line.startswith("Event"):
            psi_header=line
            samples=line.strip().split("\t")[1:]
            continue

        id = line.split("\t")[0]
        entries=line.split("\t")[1:]
        """
        #count how many are not ND or NE
        valid_count=0
        valid_list=[]
        for e in entries:
            if e!="NAN":
                valid_count+=1
                if e not in valid_list:
                    valid_list.append(e)
        
        #see if the line has more than 100 valid entries
        if valid_count<100:
            #print("removed line " + id + " with only " + str(valid_count) + " valid entries")
            continue

        #If there is only one genotype for this location, no test possible, skip.
        if len(valid_list)==1:
            #print("removed line " + id + " with only one unique value")
            continue
        
        if id not in paired_events:
            #print("removed line " + id + " not in pairs")
            continue
        """
        psi_dict[id]=dict()
        for i in range(len(entries)):
            psi_dict[id][samples[i]]=entries[i]

#Check if both genotype and psi dictionaries have entries for the same sample
result_dict=dict() 
for i in all_pairs:
    result_dict["-".join(i)] = dict()
    for sample in genotype_dict[i[1]]:
        if genotype_dict[i[1]][sample]!="ND" and genotype_dict[i[1]][sample]!="NE":
            if sample in psi_dict[i[0]] and psi_dict[i[0]][sample]!="NAN":
                result_dict["-".join(i)][sample] = [genotype_dict[i[1]][sample], psi_dict[i[0]][sample]]

#Now make sure that every pair has 100 samples left.
final_dict=dict()
for i in result_dict:
    if len(list(result_dict[i].keys()))>100:
        unique_genotypes=[]
        unique_psi=[]
        for sample in result_dict[i]:
            if result_dict[i][sample][0] not in unique_genotypes:
                unique_genotypes.append(result_dict[i][sample][0])
            if result_dict[i][sample][1] not in unique_psi:
                unique_psi.append(result_dict[i][sample][1])
        
        if len(unique_genotypes)>1 and len(unique_psi)>1:
            final_dict[i]=result_dict[i]

with open(sys.argv[4], "w") as out:
    out.write("AS_event\tVariant\tSample\tGenotype\tPSI\n")
    for i in final_dict:
        event=i.split("-")[0]
        variant=i.split("-")[1]
        for sample in final_dict[i]:
            out.write(event+"\t"+variant+"\t"+sample+"\t"+final_dict[i][sample][0]+"\t"+final_dict[i][sample][1]+"\n")

"""
updated_events=[]
updated_variants=[]
with open(sys.argv[6], "w") as out_pairs:
    out_pairs.write(pairs_header)
    for i in all_pairs:
        if i[0] in psi_dict and i[1] in genotype_dict:
            out_pairs.write(i[0]+"\t"+i[1]+"\n")
            updated_events.append(i[0])
            updated_variants.append(i[1])

with open(sys.argv[4], "w") as out_genotype:
    out_genotype.write(genotype_header)
    for k in genotype_dict:
        if k in updated_variants:
            out_genotype.write(genotype_dict[k])

with open(sys.argv[5], "w") as out_psi:
    out_psi.write(psi_header)
    for k in psi_dict:
        if k in updated_events:
            out_psi.write(psi_dict[k])


#I forgot a thing! They also need to have the values in the same samples....
"""