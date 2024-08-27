# 09.08.2024
# Mirjam Müller
# PSI_Plot.R
# make psi score distribution plots.

if (!require("tidyverse")) install.packages("tidyverse")
if (!require("gridExtra")) install.packages("gridExtra")

library(gridExtra)
library(tidyverse)

#make figure
IR_medians<- read_tsv("IR_medians.tsv")
CE_medians<- read_tsv("CE_medians.tsv")
AA_medians<- read_tsv("AA_medians.tsv")
AD_medians<- read_tsv("AD_medians.tsv")


p1<-ggplot(IR_medians, aes(x=average_psi))+
  geom_histogram()+
  theme_classic()+
  ggtitle("Intron Retention")+
  labs(x="Median PSI per event", y="Counts")

p2<-ggplot(CE_medians, aes(x=average_psi))+
  geom_histogram()+
  theme_classic()+
  ggtitle("Cassette Exons")+
  labs(x="Median PSI per event", y="Counts")

p3<-ggplot(AA_medians,aes(x=average_psi))+
  geom_histogram()+
  theme_classic()+
  ggtitle("Alternative Acceptors")+
  labs(x="Median PSI per event", y="Counts")

p4<-ggplot(AD_medians, aes(x=average_psi))+
  geom_histogram()+
  theme_classic()+
  ggtitle("Alternative Donors")+
  labs(x="Median PSI per event", y="Counts")

final_plot<- grid.arrange(p1, p2, p3, p4, ncol=2)
ggsave(paste0("PSI_plot_", Sys.Date(), ".png"), plot = final_plot, width = 7, height = 4)