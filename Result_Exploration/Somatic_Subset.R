# 09.08.2024
# Mirjam Müller
# Somatic_Subset.R


if (!require("tidyverse")) install.packages("tidyverse")

library(tidyverse)


#Set up input file, germline variants
args <- commandArgs(trailingOnly=T)

frequent<-read_tsv(args[1], skip=1)
lines <- args[2]

frequent %>% 
  sample_n(lines) %>% 
  write_tsv("somatic_subset_genotypes.tsv")

