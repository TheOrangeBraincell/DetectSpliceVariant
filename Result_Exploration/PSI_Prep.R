# 20.08.2024
# Mirjam Müller
# PSI_Prep.R
# Create median tables for each event type separately

if (!require("tidyverse")) install.packages("tidyverse")

library(tidyverse)

args <- commandArgs(trailingOnly=T)

as_type <- args[1]

#Create Median tables

##
psi_scores<- read_tsv(paste0(as_type, "_summary.tsv"))

psi_scores %>% 
  filter(n_distinct(PSI)>1) %>%
  group_by(Event) %>% 
  summarize(average_psi = median(PSI)) %>% 
  write_tsv(paste0(as_type,"_medians.tsv"))

