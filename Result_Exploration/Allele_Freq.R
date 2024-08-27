# 13.08.2024
# Mirjam Müller
# Allele_Freq2.R


if (!require("tidyverse")) install.packages("tidyverse")
if (!require("viridis")) install.packages("viridis")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("hrbrthemes")) install.packages("hrbrthemes")

library(viridis)
library(hrbrthemes)
library(gridExtra)
library(tidyverse)


#Set up input file, germline variants
args <- commandArgs(trailingOnly=T)

#Frequent germline variants:
germline <- read_tsv(args[1])
somatic<-read_tsv(args[2])

germline %>% 
  mutate(Variant_Type="germline") %>% 
  mutate(Allele_Frequency_Observed= as.double(Allele_Frequency_Observed)) %>% 
  mutate(Number_Samples=as.double(Number_Samples))-> germ_temp

somatic %>% 
  mutate(Variant_Type="somatic")-> som_temp

plot_data<- rbind(germ_temp, som_temp)

p1= ggplot(data=plot_data, aes(x=Allele_Frequency_Observed, fill=Variant_Type)) +
  geom_histogram(color="#e9ecef", binwidth=0.05, alpha=.4, position="identity") +
  theme_classic()+
  scale_fill_manual(values=c("#DAA520", "#404080")) +
  labs(fill="")+
  scale_y_log10()+
  theme(legend.position = "none")

p2 =ggplot(data = plot_data, aes(x = Allele_Frequency_Observed, fill = Variant_Type)) +
  geom_density(color = "#e9ecef", alpha = 0.4, position = "identity") +
  theme_classic() +
  scale_fill_manual(values=c("#DAA520", "#404080")) +
  #labs(fill = "")+
  ylim(0, 200)

final_plot<-grid.arrange(p1, p2, ncol=2)
ggsave(paste0("AF_plot_", Sys.Date(), ".png"), plot = final_plot, width = 7, height = 4)