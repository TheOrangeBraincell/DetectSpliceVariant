# 08.08.24
# HW_Plot.R
# Mirjam Müller
#
# Reads in all genotypes that can be associated with a dbsnp id. (This is done with a different script)
# Then does a HW test for each location separately and writes the results into a tsv
# these are then to be parsed with a python script. And then read into HW_Plot_PT2.R
# 
#
# IMPORTS
if (!require("genetics")) install.packages("genetics")
if (!require("tidyverse")) install.packages("tidyverse")

library(genetics)
library(tidyverse)

#Set up input file, germline variants
args <- commandArgs(trailingOnly=T)

input_file <- args[1]

germline_variants <- read_tsv(input_file, col_names=T)

germline_variants %>% 
  #drop_na() %>% 
  select(Location) %>% 
  distinct(Location)-> Locations

Locations <- list(Locations$Location)

result_df<- data.frame(matrix(ncol=2, nrow=0))

for (l in Locations[[1]]) {
  germline_variants %>% 
    filter(Location == l) %>% 
    pivot_longer(cols = !c(Location, Gene), names_to = "Sample", values_to = "Genotype") %>% 
    filter(Genotype != "NE" & Genotype != "ND") %>% 
    drop_na() -> temp
  
  if (nrow(temp) == 0) {
    cat(paste0("No test for location ", l, "\n"))
    next
  } else {
    temp %>% 
      mutate(Genotype = genotype(Genotype)) %>% 
      drop_na(Genotype) %>% 
      select(c(Genotype, Gene, Sample)) -> genotypes
    
    # Check if there's more than one genotype type
    if (length(unique(genotypes$Genotype)) > 1) {
      result <- capture.output(HWE.test(genotypes))
      
      # Combine the output into a single string
      result_string <- paste(result, collapse = " ")
      
      # Use regex to extract the p-value
      p_value <- str_extract(result_string, "(?<=p-value = )[0-9.eE+-]+")
      
      if (!is.na(p_value)) {
        new_row <- data.frame(Location = l, P_Value = as.numeric(p_value))
        result_df <- rbind(result_df, new_row)
      } else {
        cat(paste0("No valid p-value for location ", l, "\n"))
      }
    } else {
      cat(paste0("Only one genotype at location ", l, " - skipping HWE test\n"))
    }
  }
}


# colnames(result_df)<-c("Location", "Test_Result")
# 
# for (l in Locations[[1]]) {
#   germline_variants %>% 
#     filter(Location==l) %>% 
#     pivot_longer(cols=!c(Location, Gene), names_to="Sample", values_to="Genotype") %>% 
#     filter(Genotype != "NE" & Genotype != "ND") %>% 
#     drop_na() -> temp
#     { 
#       if (nrow(temp) == 0) {
#         cat(paste0("No test for location ",l, "\n"))
#         next
#       } else{
#         temp %>% 
#           mutate(Genotype=genotype(Genotype)) %>% 
#           drop_na(Genotype) %>% 
#           select(c(Genotype, Gene, Sample)) -> genotypes
#         result<-capture.output(HWE.test(genotypes))
#         new_row<-data.frame(l, result)
#         result_df <-rbind(result_df, new_row)
#       }
#     }
# }

#Write output file
output_file <- paste0("HW_results_", Sys.Date(), ".tsv")
result_df %>% 
  write_tsv(output_file)

