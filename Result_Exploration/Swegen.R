# 09.08.2024
# Swegen comparison full genome
# Mirjam Müller
#
#
#
# IMPORT 
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)

#Note that this one has the input files hardcoded cause there is just too many.
#I will hard code them in the pre scripts as well, just to make sure it all runs smoothly.
data <- read_tsv("swegen_scanb_variants.tsv", col_names = T)

data %>% 
  #The other 5000 have several AF for several nucleotides. but we have enough data points anyways, so i skip them. could be extracted tho.
  filter(alt_SW==alt_SC) %>% 
  filter(!grepl(",", AF_SW)) %>%
  mutate(AF_SW=as.numeric(AF_SW))-> data_alt 

p1 <- data_alt %>% 
  ggplot(aes(x=AF_SC, y=AF_SW)) +
  geom_jitter()+
  theme_classic()

ggsave(paste0("swegen_", Sys.Date(), ".png"), plot = p1, width = 7, height = 4)


#Do this with a random subset
set.seed(12345)
selected_rows<-sample(nrow(data_alt), 8000)
random_subset<-data_alt[selected_rows,]

p2<- random_subset %>% 
  ggplot(aes(x=AF_SC, y=AF_SW)) +
  geom_jitter()+
  theme_classic()

ggsave(paste0("swegen_subset_", Sys.Date(), ".png"), plot = p2, width = 7, height = 4)

#Linear regression plot
fit1<-lm(AF_SW ~AF_SC, data=data_alt)


ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(alpha=0.1) +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)),
         x = "Allele Frequency SCANB", y="Allele Frequency SWEGEN")+
    theme_classic()+
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14))
}

p3<-ggplotRegression(fit1)
ggsave(paste0("swegen_linreg_", Sys.Date(), ".png"), plot = p3, width = 7, height = 4)


#Same for more frequent variant for better resolution. First for 1 % of samples, i.e. 345 samples with genotype.
data2<-read_tsv("~/Documents/PhD/DetectSpliceVariant/Results/Result_Exploration/SweGen_Comparison/swegen_scanb345_variants.tsv", col_names=T)


data2 %>% 
  #The other 5000 have several AF for several nucleotides. but we have enough data points anyways, so i skip them. could be extracted tho.
  filter(alt_SW==alt_SC) %>% 
  filter(!grepl(",", AF_SW)) %>%
  mutate(AF_SW=as.numeric(AF_SW))-> data2_alt 

#Linear regression plot
fit2<-lm(AF_SW ~AF_SC, data=data2_alt)


p4<-ggplotRegression(fit2)
ggsave(paste0("swegen_linreg_345_", Sys.Date(), ".png"), plot = p4, width = 7, height = 4)

#That is MUCH better.

#Lets also see it with a cutoff of 100 samples, to have a minimum AF of 0.5 %
data3<-read_tsv("~/Documents/PhD/DetectSpliceVariant/Results/Result_Exploration/SweGen_Comparison/swegen_scanb100_variants.tsv", col_names=T)


data3 %>% 
  #The other 5000 have several AF for several nucleotides. but we have enough data points anyways, so i skip them. could be extracted tho.
  filter(alt_SW==alt_SC) %>% 
  filter(!grepl(",", AF_SW)) %>%
  mutate(AF_SW=as.numeric(AF_SW))-> data3_alt 

#Linear regression plot
fit3<-lm(AF_SW ~AF_SC, data=data3_alt)


p5<-ggplotRegression(fit3)
ggsave(paste0("swegen_linreg_100_", Sys.Date(), ".png"), plot = p5, width = 7, height = 4)


# We want to know if there is an overrepresentation of one variant genotype over the other. so we plot majorities.
data4<-read_tsv("~/Documents/PhD/DetectSpliceVariant/Results/Result_Exploration/SweGen_Comparison/swegen_scanb345_variants.tsv", col_names=T)


data4 %>% 
  #The other 5000 have several AF for several nucleotides. but we have enough data points anyways, so i skip them. could be extracted tho.
  filter(alt_SW==alt_SC) %>% 
  filter(!grepl(",", AF_SW)) %>%
  mutate(AF_SW=as.numeric(AF_SW))-> data4_alt 

#Linear regression plot
fit4<-lm(AF_SW ~AF_SC, data=data4_alt)



p6<- ggplot(data4_alt, aes(x = AF_SC, y = AF_SW, color=Majority_VarAllele)) + 
  geom_point(alpha=0.1) +
  stat_smooth(method = "lm", col = "red") +
  labs(title = paste("Adj R2 = ",signif(summary(fit4)$adj.r.squared, 5),
                     "Intercept =",signif(fit4$coef[[1]],5 ),
                     " Slope =",signif(fit4$coef[[2]], 5),
                     " P =",signif(summary(fit4)$coef[2,4], 5)),
       x = "Allele Frequency SCANB", y="Allele Frequency SWEGEN")+
  theme_classic()+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

ggsave(paste0("swegen_linreg_345_col_", Sys.Date(), ".png"), plot = p6, width = 7, height = 4)

#We dont seem to see any overrepresentation at all... 

# Lets read in read depth for the variants.
read_depth <- read_tsv("~/Documents/PhD/DetectSpliceVariant/Results/Result_Exploration/SweGen_Comparison/swegen_scanb_rd.tsv")

#Make comparison between "belly" and "line"
p7<-read_depth %>% 
  mutate(AF_SW=as.numeric(AF_SW)) %>% 
  drop_na() %>% 
  filter(AF_SC>= 0.1) %>% #Remove low allele frequency variants, as they are likely inaccurate due to location in exons.
  mutate(group=if_else(AF_SC/AF_SW >= 1.5, "down", "up")) %>% 
  ggplot(aes(x=as.factor(group), y = ReadDepth, fill=group))+
  geom_boxplot()+
  scale_y_log10(breaks=c(0, 0.1, 1, 10, 100, 1000))+
  theme_classic()+
  theme(legend.position="none")+
  xlab("Group")+
  ylab("Read Depth")+
  stat_summary(
    fun = mean,
    geom = "text",
    aes(label = sprintf("%.2f", ..y..)),
    position = position_nudge(x = 0.25, y = -0.1),
    color = "black",
    size = 4
  )

ggsave(paste0("swegen_readdepth_", Sys.Date(), ".png"), plot = p7, width = 7, height = 4)
