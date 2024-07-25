# Exploratory Analysis of Pipeline results
# 07.02.24
# Mirjam Müller

library(patchwork)
library(genetics)
library(GenomicRanges)
library(tidyverse)

setwd("/home/mirjam/Documents/PhD/DetectSpliceVariant/Results/Result_Exploration")

##########################################################################################################
#Lets do an overview of how many variants and events we have in those genes.
psi_counts<-read_delim("AS_counts.txt", col_names = FALSE)
#genotype_counts<-read_delim("variant_counts.txt", col_names = FALSE)
genotype_counts<-read_delim("exonic_counts.txt", col_names = FALSE)

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

  # genotype_counts %>%
  # filter(X1>7) %>% #to remove empty files
  # mutate(counts=X1-7) %>% #cause theres a file header.
  # mutate(gene=X2) %>%
  # filter(gene!="total") %>%
  # select(!c(X1, X2)) %>%
  # summarize(mean_counts = mean(counts),
  #           SD_counts = sd(counts),
  #           median_counts = median(counts),
  #           IQR_counts = IQR(counts),
  #           min_counts = min(counts),
  #           max_counts = max(counts),
  #           total_counts=sum(counts))


#use only exonic variants instead
genotype_counts
genotype_counts %>% 
  filter(X1>1) %>% #to remove empty files
  mutate(counts=X1-1) %>% #cause theres a file header.
  mutate(gene=X2) %>% 
  drop_na() %>%  # To remove total line at the end, doesnt have a gene, thus NA 
  select(!c(X1, X2)) %>% 
  summarize(mean_counts = mean(counts),
            SD_counts = sd(counts),
            median_counts = median(counts),
            IQR_counts = IQR(counts),
            min_counts = min(counts),
            max_counts = max(counts),
            total_counts=sum(counts)) %>% 



  
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
setwd("/home/mirjam/Documents/PhD/DetectSpliceVariant/Results/Result_Exploration/")
variants <- read_tsv("germline_genotypes.tsv", col_names =T)
#Do they follow a Hardy Weinberg distribution?
# Theres an R package for that!
# Count genotypes
# str(variants)
# head(variants) %>% 
#   mutate(Genotype=Gentoype) %>% 
#   select(!Gentoype) %>% 
#   filter(Genotype!="NE") %>%  
#   filter(Genotype!="ND") %>% 
#   select(c(Location,Genotype)) %>% 
#   mutate(Genotype=genotype(Genotype))-> temp_genotype

variants %>% 
  drop_na() %>% 
  select(Location) %>% 
  distinct(Location) -> Locations

Locations<- list(Locations$Location)

result_df<- data.frame(matrix(ncol=2, nrow=0))
colnames(result_df)<-c("Location", "Test_Result")

for (l in Locations[[1]]) {
  variants %>% 
    filter(Location==l) %>% 
    pivot_longer(cols=!c(Location, Gene), names_to="Sample", values_to="Genotype") %>% 
    filter(Genotype != "NE" & Genotype != "ND") %>% 
    drop_na() %>% 
    #{ print(head(.)); . } %>%
    mutate(Genotype=genotype(Genotype)) %>% 
    drop_na(Genotype) %>% 
    select(c(Genotype, Gene, Sample)) -> genotypes
  result<-capture.output(HWE.test(genotypes))
  new_row<-data.frame(l, result)
  result_df <-rbind(result_df, new_row)
}


str(result_df)
result_df
result_df %>% 
     write_tsv("HWresults_250724.tsv")

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

# Lets read in read depth for the variants.
read_depth <- read_tsv("~/Documents/PhD/DetectSpliceVariant/Results/Result_Exploration/SweGen_Comparison/swegen_scanb_rd.tsv")

#Make comparison between "belly" and "line"
read_depth %>% 
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
    

##########################################################################################################
# make psi score distribution plots.
setwd("~/Documents/PhD/DetectSpliceVariant/Results/Result_Exploration/PSI")

#One AS event type at a time, replace event abreviations, remove slice for anything thats not IR, CE.

## IR 
psi_scores<- read_tsv("IR_summary.tsv")

psi_scores %>% 
  slice_sample(n=10000000) %>% #10'000'000
  ggplot(aes(x=PSI))+
  geom_histogram()+
  scale_y_log10()+
  theme_classic()+
  ggtitle("PSI distribution IR events, subset 10'000'000")

#still need to do distribution plots that do not contain all events.

psi_scores %>% 
  slice_sample(n=10000000) %>% #10'000'000
  group_by(Event) %>% 
  filter(n_distinct(PSI) > 1) %>%
  ungroup() %>%
  ggplot(aes(x=PSI))+
  geom_histogram()+
  scale_y_log10()+
  theme_classic()+
  ggtitle("PSI distribution varying IR events, subset 10'000'000")

rm(psi_scores)

## CE

psi_scores<- read_tsv("CE_summary.tsv")

psi_scores %>% 
  slice_sample(n=10000000) %>% #10'000'000
  ggplot(aes(x=PSI))+
  geom_histogram()+
  scale_y_log10()+
  theme_classic()+
  ggtitle("PSI distribution CE events, subset 10'000'000")

psi_scores %>% 
  slice_sample(n=10000000) %>% #10'000'000
  group_by(Event) %>% 
  filter(n_distinct(PSI)>1) %>% 
  ungroup() %>% 
  ggplot(aes(x=PSI))+
  geom_histogram()+
  scale_y_log10()+
  theme_classic()+
  ggtitle("PSI distribution varying CE events, subset 10'000'000")

rm(psi_scores)
## AA
psi_scores<- read_tsv("AA_summary.tsv")

psi_scores %>% 
  ggplot(aes(x=PSI))+
  geom_histogram()+
  scale_y_log10()+
  theme_classic()+
  ggtitle("PSI distribution AA events")

rm(psi_scores)
## AD
psi_scores<- read_tsv("AD_summary.tsv")

psi_scores %>% 
  ggplot(aes(x=PSI))+
  geom_histogram()+
  scale_y_log10()+
  theme_classic()+
  ggtitle("PSI distribution AD events")

rm(psi_scores)


### Helena suggests we plot psi mean (per event) distribution as well. 
setwd("~/Documents/PhD/DetectSpliceVariant/Results/Result_Exploration/PSI")

#One AS event type at a time, replace event abreviations, remove slice for anything thats not IR, CE.

## IR 
psi_scores<- read_tsv("IR_summary.tsv")

psi_scores %>% 
  group_by(Event) %>% 
  summarize(average_psi = mean(PSI)) %>% 
  ggplot(aes(x=average_psi))+
  geom_histogram()+
  theme_classic()+
  ggtitle("Mean PSI distribution IR events")

rm(psi_scores)

## CE

psi_scores<- read_tsv("CE_summary.tsv")

psi_scores %>% 
  group_by(Event) %>% 
  summarize(average_psi = mean(PSI)) %>% 
  ggplot(aes(x=average_psi))+
  geom_histogram()+
  theme_classic()+
  ggtitle("Mean PSI distribution CE events")

rm(psi_scores)
## AA
psi_scores<- read_tsv("AA_summary.tsv")

psi_scores %>% 
  group_by(Event) %>% 
  summarize(average_psi = mean(PSI)) %>% 
  ggplot(aes(x=average_psi))+
  geom_histogram()+
  theme_classic()+
  ggtitle("Mean PSI distribution AA events")

rm(psi_scores)
## AD
psi_scores<- read_tsv("AD_summary.tsv")

psi_scores %>% 
  group_by(Event) %>% 
  summarize(average_psi = mean(PSI)) %>% 
  ggplot(aes(x=average_psi))+
  geom_histogram()+
  theme_classic()+
  ggtitle("Mean PSI distribution AD events")

rm(psi_scores)

#### AND MEDIAN PLOTS
setwd("~/Documents/PhD/DetectSpliceVariant/Results/Result_Exploration/PSI")

#Make it a 4 panel figure
library(gridExtra)
library(tidyverse)

## IR 
psi_scores<- read_tsv("IR_summary.tsv")

psi_scores %>% 
  group_by(Event) %>% 
  filter(n_distinct(PSI)>1) %>%
  summarize(average_psi = median(PSI)) %>% 
  write.table("IR_medians.tsv")

rm(psi_scores)
gc()
## CE

psi_scores<- read_tsv("CE_summary.tsv")

psi_scores %>% 
  select(Event, PSI) %>% 
  group_by(Event) %>% 
  filter(n_distinct(PSI)>1) %>% 
  summarize(average_psi = median(PSI)) %>% 
  write.table("CE_medians.tsv")

rm(psi_scores)
gc()

## AA
psi_scores<- read_tsv("AA_summary.tsv")

psi_scores %>% 
  group_by(Event) %>% 
  summarize(average_psi = median(PSI)) %>% 
  select(average_psi)%>% 
  write.table("AA_medians.tsv")

rm(psi_scores)
## AD
psi_scores<- read_tsv("AD_summary.tsv")

psi_scores %>% 
  group_by(Event) %>% 
  summarize(average_psi = median(PSI)) %>% 
  select(average_psi)%>% 
  write.table("AD_medians.tsv")

rm(psi_scores)

#make figure
IR_medians<- read.table("IR_medians.tsv")
CE_medians<- read.table("CE_medians.tsv")
AA_medians<- read.table("AA_medians.tsv")
AD_medians<- read.table("AD_medians.tsv")


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

grid.arrange(p1, p2, p3, p4, ncol=2)

##########################################################################################################

#Allele Frequency plots

setwd("/home/mirjam/Documents/PhD/DetectSpliceVariant/Results/Result_Exploration")

dbsnp <- read_tsv("germline_variants.tsv")

frequent<-read_tsv("frequent_variants_180724.tsv", skip=1)

#Take common data

dbsnp %>% 
  select(var_ID) %>% 
  distinct(var_ID) %>% 
  left_join(frequent, by = c("var_ID"="Location")) %>% 
  rename("Location"="var_ID") %>% 
  write_tsv("germline_genotypes.tsv")

frequent %>% 
  sample_n(252225) %>% 
  write_tsv("somatic_subset_genotypes.tsv")

rm(frequent)
rm(dbsnp)
gc()
#restart the session too otherwise it doesnt want to work.

germline <- read_tsv("germline_AF.tsv")
somatic <- read_tsv("somatic_AF.tsv")

germline %>% 
  mutate(Variant_Type="germline")-> germ_temp

somatic %>% 
  mutate(Variant_Type="somatic")-> som_temp

plot_data<- rbind(germ_temp, som_temp)

library(viridis)
library(hrbrthemes)
library(gridExtra)

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

grid.arrange(p1, p2, ncol=2)
