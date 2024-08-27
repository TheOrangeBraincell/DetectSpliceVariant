# 08.08.2024
# HW_Plot_PT2.R
# Mirjam Müller
#
# Takes the parsed HW test results and makes a plot showing the number of
# significant tests.
#
# IMPORTS
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

#Set up input file, germline variants
args <- commandArgs(trailingOnly=T)

input_file <- args[1]
number_of_tests <- args[2]

hw_results <- read_tsv(input_file)

#Plot the p-value distribution
hw_results %>%
  drop_na() %>% 
  mutate(p_adjusted= p.adjust(`P_Value`, n=as.integer(number_of_tests))) -> hw_results

num_below_0_05 <- sum(hw_results$p_adjusted < 0.05, na.rm = TRUE)
total <- nrow(hw_results)

p <- hw_results %>% 
  ggplot( aes(x=p_adjusted)) +
  geom_histogram(binwidth= 0.025, fill="#69b3a2", color="black", alpha=0.9) +
  ggtitle("Corrected p-values Hardy-Weinberg tests for variants found in more than 1% of the samples.") +
  theme_classic() +
  theme(
    plot.title = element_text(size=9)
  )+
  geom_vline(xintercept =0.05, color="red")+
  annotate("text", x = 0.75, y = Inf, label = paste0("p_adjusted < 0.05: ", num_below_0_05, 
                                                     "\nTotal p_adjusted: ", total), 
           hjust = 0.5, vjust = 2, size = 3, color = "black", fontface = "bold")

ggsave(paste0("hw_results_plot_", Sys.Date(), ".png"), plot = p, width = 7, height = 4)