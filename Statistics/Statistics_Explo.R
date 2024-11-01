#Statistics Exploration
# 03.10.2024
# Mirjam Müller

setwd("~/Documents/PhD/DetectSpliceVariant/Results/StatisticalComparison/")

library(tidyverse)
library(boot)

psi<-read_tsv("Test_Set/A1BG_PSI.tsv")
genotypes <- read_tsv("Test_Set/A1BG_exgeno.tsv")

pairs <- read_tsv("Test_Set/A1BG_pairs.tsv")

#For each pair I need to conduct a test. This means i need to have the data set up the right way.

#Maybe preprocess psi and genotypes and then iterate through pairs?
#Requires pairs to be in one pair per line format...

long_pairs <- pairs %>%
  separate_rows(variants, sep = ",") %>% 
  drop_na()

#Prepare psi and genotypes

psi %>% 
  pivot_longer(cols = -"Event", names_to="Sample", values_to="PSI", values_drop_na = T) -> psi_longer

genotypes %>% 
  pivot_longer(cols = -"Location", names_to="Sample", values_to="genotype") %>% 
  filter(genotype != "NE") %>% 
  filter(genotype != "ND") -> genotype_longer

#Okay the tables are prepared. But the question does remain how i will upscale and implement this...

#Lets first test it tho to know exactly in what shape i need it.

long_pairs %>% 
  head(1) %>% 
  left_join(psi_longer, by=c("AS_event"="Event")) %>% 
  inner_join(genotype_longer, by= c("variants"="Location", "Sample"))-> test_frame


# Step 1: Kruskal-Wallis Test
kruskal_test <- kruskal.test(PSI ~ genotype, data = test_frame)
print(kruskal_test)

# Step 2: Calculate Group Medians
group_medians <- test_frame %>%
  group_by(genotype) %>%
  summarise(median_value = median(PSI))
print(group_medians)

# Step 3: Define a function to calculate the difference in medians
median_diff <- function(data, indices, group1, group2) {
  boot_data <- data[indices, ]
  group1_median <- median(boot_data$PSI[boot_data$genotype == group1])
  group2_median <- median(boot_data$PSI[boot_data$genotype == group2])
  return(group1_median - group2_median)
}

# Step 4: Perform bootstrapping for confidence intervals of median differences
set.seed(123) # For reproducibility

# Example for Group A vs Group B
boot_AB <- boot(data = test_frame, statistic = median_diff, R = 1000, group1 = "0/0", group2 = "0/1")
boot.ci(boot_AB, type = "perc")

#Okay it works! Note that I need to do more pairwise comparisons if i have all 3 genotypes. 
#I also need to have a regulator on how many pairs and for which genotypes depending on the event variant pair.

#But its working! Now we need to implement it on a bigger scale and with more flexibility.

counts <- read_tsv("sample_counts.txt")

# Create bins for the count variable
counts <- counts %>%
  mutate(count_bin = cut(count, breaks = seq(0, 3500, by = 350), include.lowest = TRUE, right = FALSE))

# Create histograms for each filetype with binned counts
ggplot(counts, aes(x = count_bin, fill = filetype)) +
  geom_bar(position = "dodge") +
  facet_wrap(~ filetype) +
  labs(title = "Counts per Event by Score Type (Binned)",
       x = "Count Bins",
       y = "Frequency") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Clearly Helena is onto something. There are more events with less than 345 values.
#Lets zoom in a bit.

counts <- counts %>%
  filter(count<=500) %>% 
  mutate(count_bin = cut(count, breaks = seq(0, 500, by = 50), include.lowest = TRUE, right = FALSE))

# Create histograms for each filetype with binned counts
ggplot(counts, aes(x = count_bin, fill = filetype)) +
  geom_bar(position = "dodge") +
  facet_wrap(~ filetype) +
  labs(title = "Counts per Event by Score Type (Binned)",
       x = "Count Bins",
       y = "Frequency") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Seems like most of these under 350 samples are rare ones, because they are also under 50. 
#We could do a cutoff of 100 perhaps? This would make sense genotype wise. We could argue for the same 
# as we did for the result exploration. 100 samples are needed for an allele frequency of 0.1.
# It is the standard used. so it would make sense. I will suggest it to helena, including the plots.