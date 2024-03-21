#!/usr/bin/env Rscript 
#27.10.22
# Sequence Depth Script

library(tidyverse)

start_time <-Sys.time()
# setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs\\TNBC_comparison")
# 
# dna_variants<- read_tsv("tnbc_variants_fpkm.tsv")
# rna_variants_filtered<-read_tsv("rna_variants_filtered.txt")
# rna_variants_unfiltered<-read_tsv("rna_variants.txt")

#Bianca
dna_variants<- read_tsv("/proj/sens2019019/mirjam/TNBC_Comparison/tnbc_variants_fpkm.tsv")
rna_variants_filtered<-read_tsv("/proj/sens2019019/mirjam/TNBC_Comparison/rna_variants_filtered.txt")
rna_variants_unfiltered<-read_tsv("/proj/sens2019019/mirjam/TNBC_Comparison/rna_variants.txt")
alignment_paths<-read_tsv("/proj/sens2019019/mirjam/TNBC_Comparison/alignment_file_list.txt", col_names=FALSE)

#Adjust heterozygous genotype format
rna_variants_filtered$genotype<-str_replace(rna_variants_filtered$genotype, "1/0", "0/1")
rna_variants_unfiltered$genotype<-str_replace(rna_variants_unfiltered$genotype, "1/0", "0/1")


#Remove samples that are only in one cohort.
rna_variants_filtered %>% 
  select(Sample) %>% 
  distinct() %>% 
  pull(Sample)-> rna_samples

dna_variants %>% 
  select(Sample) %>% 
  distinct() %>% 
  pull(Sample)->dna_samples


rna_variants_filtered %>% 
  filter(Sample %in% dna_samples) %>% 
  filter(nchar(alt)==1 & nchar(ref)==1) %>% 
  filter(Location=="Exon") %>% 
  select(!Location)->rna_subs_f #rna substitutions filtered

rm(rna_variants_filtered)

rna_variants_unfiltered %>% 
  filter(Sample %in% dna_samples) %>% 
  filter(nchar(alt)==1 & nchar(ref)==1) %>% 
  filter(Location=="Exon") %>% 
  select(!c(Location, FPKM))->rna_subs

rm(rna_variants_unfiltered)

dna_variants %>% 
  filter(Sample %in% rna_samples) %>% 
  filter(nchar(alt)==1 & nchar(ref)==1) %>% 
  filter(Location=="Exon") %>% 
  select(!Location)->dna_subs #dna substitutions

rm(dna_variants)

#Small version: take subsample of number of samples. run for those!
#get all samples, from any of the three
dna_subs %>% 
  select(Sample) %>% 
  distinct() %>% 
  pull(Sample)->samples

#subset
subset_sample <- sample(samples, size=1, replace=F)


#Make all three sets smaller.
rna_subs %>% 
  filter(Sample %in% subset_sample) -> rna_subs

dna_subs %>% 
  filter(Sample %in% subset_sample)-> dna_subs

rna_subs_f %>% 
  filter(Sample %in% subset_sample)-> rna_subs_f

alignment_paths%>%
  mutate(Sample=substring(full_Path,72,78))-> paths

#print(paths)

rna_subs %>%
  inner_join(paths, by="Sample") %>%
  select(Sample, full_Path) %>%
  distinct()-> test

rna_subs_f%>% 
  mutate(Match="Filtered") %>% 
  right_join(rna_subs) %>% 
  mutate_at(vars(Match), ~replace_na(., "unfiltered")) %>% 
  right_join(dna_subs) %>% 
  mutate_at(vars(Match), ~replace_na(.,"absent")) %>% 
  #filter(FPKM>=10) %>% 
  inner_join(paths, by="Sample") %>% 
  mutate(Start= position-1) %>%  #Because same format as bed file, i.e. 0 based.
  mutate(Stop= position) %>% 
  mutate(cmd_args = paste0("depth -r ", chrom, ":", Start, "-", Stop, " ", full_Path )) %>% 
  mutate(depth = map(cmd_args, ~system2(command = "samtools", args = .x, stdout = T))) %>%
  unnest(depth) %>%
  separate(depth, into = c("chr", "nt", "depth"), extra = "merge", sep = "\t") %>% 
  select(-cmd_args, -chr, -Start, -Stop, -full_Path, -Location, -nt, -FPKM) %>% 
  group_by(Sample) %>% 
  group_walk(~write_tsv(.x, paste0(.y$Sample, "_seq_deph.tsv")))

stop_time <-Sys.time()   

print(stop_time-start_time)