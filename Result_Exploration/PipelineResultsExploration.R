# Exploratory Analysis of Pipeline results
# 07.02.24
# Mirjam Müller

library(patchwork)
library(genetics)
library(tidyverse)
library(GenomicRanges)

setwd("/home/mirjam/Documents/PhD/Result_Exploration")

##########################################################################################################
#Lets do an overview of how many variants and events we have in those genes.
psi_counts<-read_delim("AS_counts.txt", col_names = FALSE)
genotype_counts<-read_delim("variant_counts.txt", col_names = FALSE)

psi_counts
psi_counts %>% 
  filter(X1>1) %>% #to remove empty files
  mutate(counts=X1-1) %>% #cause theres a file header.
  filter(X2!="total") %>%  #theres a sum line at the end
  mutate(gene=X2) %>% 
  select(!c(X1, X2)) %>% 
  summarize(mean_counts = mean(counts),
            SD_counts = sd(counts),
            median_counts = median(counts),
            IQR_counts = IQR(counts),
            min_counts = min(counts),
            max_counts = max(counts),
            total_counts=sum(counts))

  genotype_counts %>% 
  filter(X1>7) %>% #to remove empty files
  mutate(counts=X1-7) %>% #cause theres a file header.
  mutate(gene=X2) %>% 
  filter(gene!="total") %>% 
  select(!c(X1, X2)) %>% 
  summarize(mean_counts = mean(counts),
            SD_counts = sd(counts),
            median_counts = median(counts),
            IQR_counts = IQR(counts),
            min_counts = min(counts),
            max_counts = max(counts),
            total_counts=sum(counts))

#Not same amount of genes. check
genotype_counts %>% 
  filter(X1>7) %>% #to remove empty files
  mutate(counts=X1-7) %>% #cause theres a file header.
  mutate(gene=X2) %>% 
  filter(gene!="total") %>% 
  select(!c(X1, X2)) -> genotypes_alt

psi_counts %>% 
  filter(X1>1) %>% #to remove empty files
  mutate(counts=X1-1) %>% #cause theres a file header.
  filter(X2!="total") %>%  #theres a sum line at the end
  mutate(gene=X2) %>% 
  select(!c(X1, X2)) %>% 
  inner_join(genotypes_alt, by="gene")

##########################################################################################################
#Overview on exon variants

setwd("~/Documents/PhD/StatisticalComparison/")  
exon_variants<-read_tsv("variant_counts_exon.txt", col_names = F)

exon_variants %>% 
  mutate(counts=X1) %>%
  summarize(mean_counts = mean(counts),
            SD_counts = sd(counts),
            median_counts = median(counts),
            IQR_counts = IQR(counts),
            min_counts = min(counts),
            max_counts = max(counts),
            total_counts=sum(counts))


  
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
  labs(title="Variants", x="counts of events per gene", y="frequency")
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
  labs(title="AS events", x="counts of events per gene", y="frequency")

combined_plots<-plot1+plot2
print(combined_plots)
#both show exponential distribution. Very skewed.

##########################################################################################################
#Read in data.
setwd("/home/mirjam/Documents/PhD/Download_060224/")
genotype_files<-list.files("/home/mirjam/Documents/PhD/Download_060224", pattern = "genotypes.tsv", full.names = FALSE)
psi_files<-list.files("/home/mirjam/Documents/PhD/Download_060224", pattern = "PSI.tsv", full.names = FALSE)

test<-genotype_files[1:10]

# If we read them all in in one go, it dies. So we want to read them in in a loop and remove the entries we dont need as we go along!
# read_filter <- function(paths) {
#   df <- read_tsv(paths, skip=6, col_names = TRUE) #skip header lines.
#   #remove empty ones
#   if (nrow(df)==0){
#     return(NULL)
#   }
#   # Only want frequent variants.
#   df <- df %>% #remove variants which occur in less than 1 % of samples. I.e. filter less than 350
#     rowwise() %>%
#     mutate(count_specific = sum(c_across(everything()) == "1/1" | c_across(everything()) == "0/1")) %>%
#     mutate(Gene=str_extract(paths, "^[^_]+")) %>% 
#     filter(count_specific>350) %>% 
#     select(!count_specific) %>% 
#     pivot_longer(cols=!c("Location", "Gene"), names_to = "Sample", values_to="Gentoype")
#   # Return the transformed dataframe
#   return(df)
# }
# frequent_variants<- bind_rows(lapply(genotype_files, read_filter))
#
# #so we never have to do this again, write it in a file!
# setwd("/home/mirjam/Documents/PhD/")
# str(frequent_variants)
# frequent_variants %>% 
#   write_tsv("frequent_variants_070224.tsv")


##########################################################################################################
setwd("/home/mirjam/Documents/PhD/Result_Exploration")
variants <- read_tsv("frequent_variants_070224.tsv")
  #Do they follow a Hardy Weinberg distribution?
# Theres an R package for that!
# Count genotypes
str(variants)
head(variants) %>% 
  mutate(Genotype=Gentoype) %>% 
  select(!Gentoype) %>% 
  filter(Genotype!="NE") %>%  
  filter(Genotype!="ND") %>% 
  select(c(Location,Genotype)) %>% 
  mutate(Genotype=genotype(Genotype))-> temp_genotype

temp_genotype
result<-HWE.test(temp_genotype)
  
variants %>% 
  select(Location) %>% 
  distinct(Location) -> Locations

Locations<- list(Locations$Location)
Locations
variant_test<-variants
result_df<- data.frame(matrix(ncol=2, nrow=0))
colnames(result_df)<-c("Location", "Test_Result")
for (l in Locations[[1]]) {
  print(l)
  variant_test %>% 
    filter(Location==l) %>% 
    mutate(Genotype=genotype(Gentoype)) %>% 
    select(!c(Gentoype, Gene, Sample)) %>%  
    drop_na(Genotype)-> genotypes
  result<-capture.output(HWE.test(genotypes))
  new_row<-data.frame(l, result)
  result_df <-rbind(result_df, new_row)
}


str(result_df)
result_df
result_df %>% 
     write_tsv("HWresults_140224.tsv")

#Parsed the results with python. Read them back in.
hw_results <- read_tsv("HW_results_parsed.tsv")

#272 locations had no difference between alleles and therefore there could be no HW test and no p-values.
hw_results %>% 
  drop_na() %>% 
  mutate(p_adjusted= p.adjust(`p-value`, n=16624)) %>% 
  filter(p_adjusted >= 0.05)

# 6'786 significant p-values after correction for multiple testing
# 9'838 not significant after correction for multiple testing.
# This is still more significant test results than by chance.

#Plot the p-value distribution
hw_results %>%
  drop_na() %>% 
  mutate(p_adjusted= p.adjust(`p-value`, n=16624)) %>% 
  ggplot( aes(x=p_adjusted)) +
  geom_histogram(binwidth= 0.025, fill="#69b3a2", color="black", alpha=0.9) +
  ggtitle("Corrected p-values Hardy-Weinberg tests for variants found in more than 1% of the samples.") +
  theme_classic() +
  theme(
    plot.title = element_text(size=10)
  )+
  geom_vline(xintercept =0.05, color="red")


#Now looking at the variants again, we try to find the annotated ones. Which means another long reading in loop with filtering...


##########################################################################################################
# Comparing allele frequencies to known variants:
#We do ESR1 because otherwise too many variants...

#use same dataset of swegen variants as for the same comparison in master thesis, so we get a direct comparison of plots.
swegen38<-read_tsv("~/Documents/Windows_Subsystem/University/Master Thesis/R and outputs/Statistical_Test/ESR1/hglft_genome_27da8_d7c440.bed", col_names=FALSE)

swegen38 %>% 
  rename(X1="chrom") %>% 
  rename(X2="position") %>% 
  select(!c(X3, X5)) %>% 
  separate(X4, c("ref", "alt", "Allele_Frequency_Expected", "Observed_Expected"), sep="_") %>% 
  select(!Observed_Expected) %>% 
  mutate(Allele_Frequency_Expected=as.double(Allele_Frequency_Expected))->swegen_alt


#read only genotypes of ESR1

genotype<-read_tsv("~/Documents/PhD/Download_060224/ESR1_genotypes.tsv", skip=6)
genotype<- read_tsv("~/Documents/Windows_Subsystem/University/Master Thesis/R and outputs/Statistical_Test/ESR1/ESR1_genotype_table.tsv", skip=8)

genotype %>%
  #Split infostring into separate columns.
  separate(col=Location, c("other", "alt"), sep="\\(") %>%
  separate(col=other, c("chrom", "position", "ref", "empty"), sep="_") %>%
  select(!empty) %>%
  mutate(alt=str_replace(alt, "\\)", "")) %>%
  mutate(position=as.integer(position)+1) %>% #Because my file coordinates are zero based 
  mutate(hmza=rowSums(genotype=="1/1"| genotype=="2/2" | genotype =="3/3")) %>% 
  mutate(hetz=rowSums(genotype=="0/1"| genotype=="1/2"| genotype=="0/2"| genotype=="2/3" | genotype=="0/3")) %>% 
  mutate(hmzr=rowSums(genotype=="0/0")) %>% 
  mutate(total=hmza+hmzr+hetz) %>% 
  mutate("Ref"=(2*hmzr/(2*total))) %>% 
  mutate("Het"=hetz/(2*total)) %>% 
  mutate("Alt"=2*hmza/(2*total)) %>% 
  mutate('Allele_Frequency_observed'=(2*hmza + hetz)/(2*total)) %>% 
  select(c(chrom, position, ref, alt, total, Allele_Frequency_observed)) %>% 
  rename(total="Number_of_Samples")-> genotypes_expanded

#Exclude non coding variants
annotation<-read_tsv("~/Documents/Windows_Subsystem/University/Master Thesis/R and outputs/Database/gencode.v39.annotation.gff3", skip=7, col_names=FALSE)

annotation %>% 
  select(c(X1, X3, X4, X5)) %>% 
  filter(X1=="chr6") %>% 
  mutate("chrom"=X1) %>% 
  mutate("start"=X4) %>% 
  mutate("stop"=X5) %>% 
  mutate("Location"=X3) %>% 
  select(!c(X1, X3, X4, X5)) %>% 
  filter(Location=="exon")-> ann_interest

ann_interest      
rm(annotation)

gr_ann=makeGRangesFromDataFrame(ann_interest)
genotypes_expanded %>% 
  select(!c(ref, alt, Allele_Frequency_observed)) %>% 
  mutate(start=position) %>% 
  mutate(stop=position+1) %>% 
  makeGRangesFromDataFrame()-> gr_geno

coordinates<-as_tibble(findOverlaps(gr_geno, gr_ann))

genotypes_expanded

coordinates %>% 
  #mutate(ann_interest[subjectHits,1]) %>% 
  mutate(ann_interest[subjectHits,4]) %>% 
  mutate(genotypes_expanded[queryHits,1]) %>% 
  mutate(genotypes_expanded[queryHits,2]) %>% 
  mutate(genotypes_expanded[queryHits,3]) %>% 
  mutate(genotypes_expanded[queryHits,4]) %>% 
  left_join(genotypes_expanded) %>% 
  select(!c(queryHits, subjectHits)) %>% 
  distinct()-> genotypes_exons

genotypes_exons %>% 
  inner_join(swegen_alt, by=c("chrom", "position", "ref", "alt"))-> comparison_data

write_tsv(comparison_data, "SweGen_SCANB_ESR1.tsv")

comparison_data
ggplot(comparison_data, aes(x=Allele_Frequency_observed, y=Allele_Frequency_Expected)) +
  theme_classic()+
  geom_point(size=2)+
  xlab("Allele Frequency observed in SCAN-B")+
  ylab("Allele Frequency observed in SweGen")+
  theme(text = element_text(size = 15), axis.text=element_text(size=15))

#check read depth at variant locations
read_depth<-read_tsv("ESR1_read_depth.tsv")

read_depth %>% 
  select_all(~sub("^.*(S00\\d+).*", "\\1",.))-> temp

col_names_readdepth<-c("chrom", "start", "stop", "info", "quality", "strand", colnames(temp)[1:(length(colnames(temp))-1)])
#col_names_readdepth<-c("chrom", "start", "stop", "info", "quality", "strand", colnames(temp))
tail(col_names_readdepth)

read_depth<-read_tsv("ESR1_read_depth.tsv", skip=1,  col_names = col_names_readdepth)

read_depth %>% 
  select(-c(start, info, quality, strand)) %>% 
  rename("position"=stop) %>% #because 0 based bed file, Swegen is 1 based
  inner_join(comparison_data) %>% 
  select(chrom, position, ref, alt, Number_Samples, Allele_Frequency_observed, Allele_Frequency_Expected, everything()) -> comparison_w_readdepth

write_tsv(comparison_w_readdepth, "SWEGEN_SCANB_readdepth.tsv")

##########################################################################################################

#Remake same plot for all genes not just ESR1. But only include exonic variants
data<-read_tsv("~/Documents/PhD/DetectSpliceVariant/Results/Result_Exploration/SweGen_Comparison/swegen_scanb_variants.tsv", col_names=T)
data

data %>% 
  #The other 5000 have several AF for several nucleotides. but we have enough data points anyways, so i skip them. could be extracted tho.
  filter(alt_SW==alt_SC) %>% 
  filter(!grepl(",", AF_SW)) %>%
  mutate(AF_SW=as.numeric(AF_SW))-> data_alt 

data_alt %>% 
  ggplot(aes(x=AF_SC, y=AF_SW)) +
  geom_jitter()+
  theme_classic()

#Do this with a random subset
set.seed(12345)
selected_rows<-sample(nrow(data_alt), 8000)
random_subset<-data_alt[selected_rows,]

random_subset %>% 
  ggplot(aes(x=AF_SC, y=AF_SW)) +
  geom_jitter()+
  theme_classic()

#Linear regression plot
fit1<-lm(AF_SW ~AF_SC, data=data_alt)
summary(fit1)

ggplot(data_alt, aes(x = AF_SC, y = AF_SW)) + 
  geom_point(alpha=0.1) +
  stat_smooth(method = "lm", col = "red")+
  theme_classic()

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

ggplotRegression(fit1)

#Same for more frequent variant for better resolution. First for 1 % of samples, i.e. 345 samples with genotype.
data2<-read_tsv("~/Documents/PhD/DetectSpliceVariant/Results/Result_Exploration/SweGen_Comparison/swegen_scanb345_variants.tsv", col_names=T)
data2

data2 %>% 
  #The other 5000 have several AF for several nucleotides. but we have enough data points anyways, so i skip them. could be extracted tho.
  filter(alt_SW==alt_SC) %>% 
  filter(!grepl(",", AF_SW)) %>%
  mutate(AF_SW=as.numeric(AF_SW))-> data2_alt 

#Linear regression plot
fit2<-lm(AF_SW ~AF_SC, data=data2_alt)
summary(fit2)

ggplotRegression(fit2)

#That is MUCH better.

#Lets also see it with a cutoff of 100 samples, to have a minimum AF of 0.5 %
data3<-read_tsv("~/Documents/PhD/DetectSpliceVariant/Results/Result_Exploration/SweGen_Comparison/swegen_scanb100_variants.tsv", col_names=T)
data3

data3 %>% 
  #The other 5000 have several AF for several nucleotides. but we have enough data points anyways, so i skip them. could be extracted tho.
  filter(alt_SW==alt_SC) %>% 
  filter(!grepl(",", AF_SW)) %>%
  mutate(AF_SW=as.numeric(AF_SW))-> data3_alt 

#Linear regression plot
fit3<-lm(AF_SW ~AF_SC, data=data3_alt)
summary(fit3)

ggplotRegression(fit3)


# We want to know if there is an overrepresentation of one variant genotype over the other. so we plot majorities.
data4<-read_tsv("~/Documents/PhD/DetectSpliceVariant/Results/Result_Exploration/SweGen_Comparison/swegen_scanb345_variants.tsv", col_names=T)
data4

data4 %>% 
  #The other 5000 have several AF for several nucleotides. but we have enough data points anyways, so i skip them. could be extracted tho.
  filter(alt_SW==alt_SC) %>% 
  filter(!grepl(",", AF_SW)) %>%
  mutate(AF_SW=as.numeric(AF_SW))-> data4_alt 

#Linear regression plot
fit4<-lm(AF_SW ~AF_SC, data=data4_alt)
fit4$model %>% 
  left_join(data4_alt, by=c("AF_SW", "AF_SC"))

data4_alt
ggplot(data4_alt, aes(x = AF_SC, y = AF_SW, color=Majority_VarAllele)) + 
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

fit4$model

#We dont seem to see any overrepresentation at all... 

#Lets check what the variants at AF 0 are and if we can explain those somehow.
data4_alt %>% 
  filter(AF_SC<0.01) %>% 
  filter(AF_SW>0.05)# smaller than 1 % is still like 31'605 rows. Thats a lot! How do I check their expression?



# make psi score distribution plots.

#psi_files contains all the names to psi files.
# read_files <- function(paths) {
#   df <- read_tsv(paths, col_names = TRUE) 
#   #remove empty ones
#   if (nrow(df)==0){
#     return(NULL)
#   }
#   df %>% 
#     mutate(file_path=paths) %>% 
#     separate(file_path, into=c("gene", "file_end"), sep="_") %>% 
#     select(-file_end) %>% 
#     pivot_longer(cols=-c("Event", "gene"), names_to="Sample", values_to="PSI") %>% 
#     drop_na()-> df_alt
#   return(df_alt)
# }
# setwd("~/Documents/PhD/Download_060224/")
# psi_scores<- bind_rows(lapply(psi_files, read_files))



# ESR1<-read_tsv("../Download_060224/ESR1_PSI.tsv")
# 
# ESR1 %>%
#   mutate(file_path="ESR1_PSI.tsv") %>%
#   separate(file_path, into=c("gene", "file_end"), sep="_") %>%
#   select(-file_end) %>%
#   pivot_longer(cols=-c("Event", "gene"), names_to="Sample", values_to="PSI") %>%
#   drop_na() %>% 
  


psi_files
setwd("~/Documents/PhD/Result_Exploration/")
psi_scores<- read_tsv("PSI_summary.tsv")

#distribution plots


#how many genes have enough data?



