# 08.08.2024
# Summary_Table.R
# Mirjam Müller
#
# Makes a summary table over AS events and variants
#
#
# IMPORTS
if (!require("grid")) install.packages("grid")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("patchwork")) install.packages("patchwork")
if (!require("tidyverse")) install.packages("tidyverse")

library(gridExtra)
library(grid)
library(patchwork)
library(tidyverse)


args <- commandArgs(trailingOnly=T)

psi <- args[1]
genotype <- args[2]


psi_counts<-read_delim(psi, col_names = FALSE)
genotype_counts<-read_delim(genotype, col_names = FALSE)

psi_counts %>% 
  filter(X1>1) %>% #to remove empty files
  mutate(counts=X1-1) %>% #cause theres a file header.
  filter(X2!="total") %>%  #theres a sum line at the end
  mutate(gene=X2) %>% 
  select(!c(X1, X2)) %>% 
  summarize(mean_counts = round(mean(counts), 2),
            SD_counts = round(sd(counts), 2),
            median_counts = median(counts),
            IQR_counts = IQR(counts),
            min_counts = min(counts),
            max_counts = max(counts),
            total_counts=sum(counts)) -> psi_summary_table

table_grob <- tableGrob(psi_summary_table)

# Save the table as PNG
png("psi_summary_table.png", width = 800, height = 100)
grid.draw(table_grob)
dev.off()


#use only exonic variants
genotype_counts %>% 
  filter(X1>1) %>% #to remove empty files
  mutate(counts=X1-1) %>% #cause theres a file header.
  mutate(gene=X2) %>% 
  drop_na() %>%  # To remove total line at the end, doesnt have a gene, thus NA 
  select(!c(X1, X2)) %>% 
  summarize(mean_counts = round(mean(counts), 2),
            SD_counts = round(sd(counts), 2),
            median_counts = median(counts),
            IQR_counts = IQR(counts),
            min_counts = min(counts),
            max_counts = max(counts),
            total_counts=sum(counts)) -> genotype_summary_table


table_grob <- tableGrob(genotype_summary_table)

# Save the table as PNG
png("genotype_summary_table.png", width = 800, height = 100)
grid.draw(table_grob)
dev.off()

#Some overview plots too! For the distribution of variables.
plot1<-genotype_counts %>%
  filter(X1>7) %>% #to remove empty files
  mutate(counts=X1-7) %>% #cause theres a file header.
  mutate(gene=X2) %>% 
  filter(gene!="total")%>% 
  select(!c(X1, X2)) %>% 
  ggplot(aes(x=counts))+
  geom_histogram( fill="turquoise")+
  theme_classic()+
  xlim(0, 400)+
  ylim(0, 3000)+
  labs(title="Exonic Variants", x="Events per gene", y="Frequency")
plot2<-psi_counts %>% 
  filter(X1>1) %>% #to remove empty files
  mutate(counts=X1-1) %>% #cause theres a file header.
  mutate(gene=X2) %>% 
  filter(gene!="total")%>% 
  select(!c(X1, X2)) %>% 
  ggplot(aes(x=counts))+
  geom_histogram(fill="skyblue")+
  theme_classic()+
  xlim(0, 200)+
  ylim(0,3000)+
  labs(title="Alternative Splicing", x="Events per gene", y="Frequency")

combined_plots<-plot1+plot2
ggsave(paste0("Summary_plot_", Sys.Date(), ".png"), plot = combined_plots, width = 7, height = 4)