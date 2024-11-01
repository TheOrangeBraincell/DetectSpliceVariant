# 29.10.2024
# Mirjam Müller
# Kruskal Wallis tests
# KW_tests.R
# 
# Reads in a triple of files: PSI table, genotype table and the pairs file matching their events.
# These have been prefiltered with python to lessen the load.
# 
# cat test_genes.txt | while read gene; do Rscript KW_tests.R ${gene}_filt.tsv ${gene}_kw.tsv; done
#
# Does a KW test and then puts out a p-value file.

library(tidyverse)
# Input from command line
args <- commandArgs(trailingOnly = TRUE)

# 1 for input file
test_frame <- read_tsv(args[1])


# Step 1: Kruskal-Wallis Test per pair
test_frame %>% 
  group_by(Variant, AS_event) %>% 
  summarise(
    kruskal_p_value = kruskal.test(PSI ~ Genotype)$p.value,
    kruskal_statistic = kruskal.test(PSI ~ Genotype)$statistic
  )-> result

result %>% 
  arrange(kruskal_p_value) %>% 
  write.table(file=args[2], quote=F, sep="\t", row.names=F)


#I still have an error for one of the genes, that all observations are in the same group. check why.








