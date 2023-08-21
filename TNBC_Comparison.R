# TNBC Sequence Depth Plotting
# 18.08.23
# Needs to import all 238 sequencing files.

library(tidyverse)
library(viridis)

setwd("C:\\Users\\mirja\\Documents\\University\\PhD\\R\\TNBC_Comparison\\Sequence_Depths_TNBC")
 
seq_depths<- tibble(full_Path=list.files(full.names=TRUE, recursive=TRUE)) %>% 
  filter(str_detect(full_Path, pattern="_seq_deph.tsv")) %>%
  mutate(Sample=substring(full_Path,3,9)) %>% 
  mutate(df = map(full_Path, ~read_tsv(.x, show_col_types = F, col_types="cccccccc"))) %>% 
  unnest(df) %>% 
  select(!full_Path)

seq_depths

#This seems too short.
#But ive noticed when using samtools depth in python, that it just refuses to give an output if the depth is 0
# which means most likely, there is a depth of 0 for most variants (which would make sense)
# So maybe join it with the files from before we put them through samtools?
# And then every row that has NA gets a 0


#Went back to check with old file. Once merged, RNA, DNA and RNA filtered give a table with 
#about 55000 variants there. So this is probably fine.


seq_depths %>%
  mutate(depth=as.integer(depth)) %>% 
  #(depth<=75) %>%
  #filter(FPKM>10) %>% 
  ggplot(aes(x=Match, y=depth, fill=Match))+
  geom_boxplot()+
  theme_classic()+
  xlab("Matches between DNA data and RNA data.")+
  ylab("Read depth")+
  ylim(0,75)+
  ggtitle("Sequencing depth of substitutions found in")+
  scale_fill_discrete(name = "Substitutions found in", labels = c("DNA", "DNA + RNA filtered", "DNA + RNA unfiltered"))+
  theme(axis.text.x = element_blank())


#This doesnt show us super much, except that the filtered RNA variants (that are also found in DNA) have the highest average depth.
#This makes sense as one of our filter criteria is read depth. (var depth >5)

#Lets try and make a barplot as we did for FPKM.

seq_depths %>% 
  mutate(depth=as.integer(depth)) %>% 
  mutate(bins=case_when(depth>=0 & depth<=1 ~ "<1",
                        depth>1 & depth<=10 ~ "1-10",
                        depth>10 & depth<=20 ~ "11-20",
                        depth>20 & depth<=30 ~ "21-30",
                        depth>30 & depth<=40 ~ "31-40",
                        depth>40 & depth<=50 ~ "41-50",
                        depth>50 & depth<=60 ~ "51-60",
                        depth>60 & depth<=70 ~ "61-70",
                        depth>70 & depth<=80 ~ "71-80",
                        depth>80 & depth<=90 ~ "81-90",
                        depth>90 & depth<=100 ~ "91-100",
                        depth>100 ~ ">100"))  %>%
  mutate(bins=as.factor(bins)) %>% 
  mutate(bins=factor(bins,levels = c("<1","1-10", "11-20", "21-30", "31-40","41-50","51-60","61-70","71-80","81-90","91-100", ">100")))-> seq_depths_binned

ggplot(seq_depths_binned, aes(fill=factor(Match, levels=c("absent", "unfiltered", "Filtered")), x=bins))+
  theme_classic()+
  scale_fill_viridis(discrete = T, name = "Fraction of Variants found in", labels = c("DNA", "DNA + RNA unfiltered", "DNA + RNA filtered"))+
  xlab("Sequence Depth")+
  ylab("%Variants found")+
  #scale_color_brewer(palette = "Set1") + 
  #scale_fill_brewer(name = "Variants found in", labels = c("DNA", "DNA + RNA unfiltered", "DNA + RNA filtered"))+
  ggtitle("Variants found in DNA data, that are also found in RNA data")+
  geom_bar(position="fill")+
  theme(text = element_text(size = 15), axis.text=element_text(size=15))+
  #geom_text(aes(label = ..count..), size=5, stat = "count", position="fill", vjust =1.5)+
  scale_x_discrete(guide=guide_axis(n.dodge=2))


#Big differences in the 1-40 area - look at closer.

seq_depths %>% 
  mutate(depth=as.integer(depth)) %>% 
  filter(depth<30) %>%
  filter(depth>1) %>% 
  mutate(bins=case_when(depth>1 & depth<=3 ~ "1-3",
                        depth>3 & depth<=6 ~ "4-6",
                        depth>6 & depth<=9 ~ "7-9",
                        depth>9 & depth<=12 ~ "10-12",
                        depth>12 & depth<=15 ~ "13-15",
                        depth>15 & depth<=18 ~ "16-18",
                        depth>18 & depth<=21 ~ "19-21",
                        depth>21 & depth<=24 ~ "22-24",
                        depth>24 & depth<=27 ~ "25-27",
                        depth>27 & depth<=30 ~ "28-30"))  %>%
  mutate(bins=as.factor(bins)) %>% 
  mutate(bins=factor(bins,levels = c("1-3", "4-6", "7-9", "10-12","13-15","16-18","19-21","22-24","25-27","28-30")))-> seq_depths_binned

ggplot(seq_depths_binned, aes(fill=factor(Match, levels=c("absent", "unfiltered", "Filtered")), x=bins))+
  theme_classic()+
  scale_fill_viridis(discrete = T, name = "Fraction of Variants found in", labels = c("DNA", "DNA + RNA unfiltered", "DNA + RNA filtered"))+
  xlab("Sequence Depth")+
  ylab("%Variants found")+
  #scale_color_brewer(palette = "Set1") + 
  #scale_fill_brewer(name = "Variants found in", labels = c("DNA", "DNA + RNA unfiltered", "DNA + RNA filtered"))+
  ggtitle("Variants found in DNA data, that are also found in RNA data")+
  geom_bar(position="fill")+
  theme(text = element_text(size = 15), axis.text=element_text(size=15))+
  #geom_text(aes(label = ..count..), size=5, stat = "count", position="fill", vjust =1.5)+
  scale_x_discrete(guide=guide_axis(n.dodge=2))


seq_depths %>% 
  mutate(depth=as.integer(depth)) %>% 
  filter(depth<28) %>%
  filter(depth>20) %>% 
  mutate(bins=case_when(depth==21 ~"21",
                        depth>21 & depth<=22 ~ "22",
                        depth>22 & depth<=23 ~ "23",
                        depth>23 & depth<=24 ~ "24",
                        depth>24 & depth<=25 ~ "25",
                        depth>25 & depth<=26 ~ "26",
                        depth>26 & depth<=27 ~ "27"))  %>%
  mutate(bins=as.factor(bins))-> zoom_in

ggplot(zoom_in, aes(fill=factor(Match, levels=c("absent", "unfiltered", "Filtered")), x=bins))+
  theme_classic()+
  scale_fill_viridis(discrete = T, name = "Fraction of Variants found in", labels = c("DNA", "DNA + RNA unfiltered", "DNA + RNA filtered"))+
  xlab("Sequence Depth")+
  ylab("%Variants found")+
  #scale_color_brewer(palette = "Set1") + 
  #scale_fill_brewer(name = "Variants found in", labels = c("DNA", "DNA + RNA unfiltered", "DNA + RNA filtered"))+
  ggtitle("Variants found in DNA data, that are also found in RNA data")+
  geom_bar(position="fill")+
  theme(text = element_text(size = 15), axis.text=element_text(size=15))+
  #geom_text(aes(label = ..count..), size=5, stat = "count", position="fill", vjust =1.5)+
  scale_x_discrete(guide=guide_axis(n.dodge=2))
