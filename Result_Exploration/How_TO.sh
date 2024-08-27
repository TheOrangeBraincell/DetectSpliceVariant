# 08.08.2024
# Mirjam Müller
#
# Notes on how to run all the scripts needed for the data exploration plots:
conda init
conda activate base

#1. Summary table
wc -l ~/Wednesday/H730/noncoding/projects/Mirjam_Outputs_Pipeline/Genotype_Tables/* > exonic_variant_counts.txt
#Remove all whitespace at the beginning of the line in output file
sed 's/^[[:space:]]*//' exonic_variant_counts.txt > temp.txt
# format second column
sed 's/[[:space:]]/\//' temp.txt | cut -d"/" -f1,5 | cut -d"_" -f1| sed 's/\// /' > exonic_counts.txt

wc -l ~/Wednesday/H730/noncoding/projects/Mirjam_Outputs_Pipeline/AS_Events/* > AS_counts.txt
#Remove all whitespace at the beginning of the line in output file
sed 's/^[[:space:]]*//' AS_counts.txt > temp.txt
# format second column
sed 's/[[:space:]]/\//' temp.txt | cut -d"/" -f1,5 | cut -d"_" -f1| sed 's/\// /' > AS_counts.txt


Rscript R_Scripts/Summary_Table.R AS_counts.txt exonic_counts.txt

# 2. HW script

conda activate DetectSpliceVariants
#First run germline script
python Data_Exploration/germline_variants.py ~/Wednesday/H730/noncoding/projects/Mirjam_Outputs_Pipeline/Genotype_Tables/ germline_variants_270824.tsv

#Note that this step requires 00-All.vcf.gz and corresponding indexed file to be in same folder.
#If you lose it, remake it with: tabix -p vcf 00-All.vcf.gz

#I also only want to do the HW thing for frequent variants. So I also need to run frequent_variants.py. Then we merge those output files after.
python Data_Exploration/frequent_variants.py ~/Wednesday/H730/noncoding/projects/Mirjam_Outputs_Pipeline/Genotype_Tables/ frequent_variants_270824.tsv

#Merge the two! 
python Data_Exploration/frequent_germline.py germline_variants_270824.tsv frequent_variants_270824.tsv germline_genotypes_270824.tsv

conda activate base
# Run HW prep script
Rscript R_Scripts/HW_Prep.R germline_genotypes_270824.tsv > HW_Prep_log.txt  2>&1

#Check for number of tests
n=$(wc -l < HW_results_2024-08-27.tsv)

#Then run HW plot script, insert number of tests instead of 16624
Rscript R_Scripts/HW_Plot.R HW_results_2024-08-27.tsv $n


# 3. Swegen
# This one needs a lot of files and pre steps.
#Unzipping Swegen file
#gunzip swegen_frequencies_fixploidy_GRCh38_20190204.vcf.gz
#Its binary
#conda activate bcftools
#bcftools convert -O v -o swegen_variants.vcf swegen_frequencies_fixploidy_GRCh38_20190204.vcf

#conda activate DetectSpliceVariants
#Comparing Swegen file to our variants: first all variants:
python Data_Exploration/SweGen_comp.py swegen_variants.vcf ~/Wednesday/H730/noncoding/projects/Mirjam_Outputs_Pipeline/Genotype_Tables/ swegen_scanb_variants.tsv

#Now with thresholds:
python Data_Exploration/SweGen_comp.py swegen_variants.vcf ~/Wednesday/H730/noncoding/projects/Mirjam_Outputs_Pipeline/Genotype_Tables/ swegen_scanb345_variants.tsv 345
python Data_Exploration/SweGen_comp.py swegen_variants.vcf ~/Wednesday/H730/noncoding/projects/Mirjam_Outputs_Pipeline/Genotype_Tables/ swegen_scanb100_variants.tsv 100
#Zipping Swegen file again, cause she big
#gzip swegen_frequencies_fixploidy_GRCh38_20190204.vcf

# For the read depth plot:
python Data_Exploration/Swegen_readdepth.py ~/Wednesday/H730/noncoding/projects/Mirjam_Outputs_Pipeline/Read_Depth/ swegen_scanb_variants.tsv swegen_scanb_rd.tsv

#Finally plotting. Note that all file names are hardcoded in this one, because theres too many
#The required files are:
# - swegen_scanb_variants.tsv
# - swegen_scanb345_variants.tsv
# - swegen_scanb100_variants.tsv
# - swegen_scanb_rd.tsv

conda activate base

Rscript R_Scripts/Swegen.R


# 4. Allele Frequency

n=$(wc -l < germline_genotypes_270824.tsv)

# Take same size subset from frequent variants
Rscript R_Scripts/Somatic_Subset.R frequent_variants_130824.tsv $n

# Run Allele Frequency Scripts in python for both datasets
python Data_Exploration/Genotype_AF.py ~/Wednesday/H730/noncoding/projects/Mirjam_Outputs_Pipeline/Genotype_Tables/ somatic_subset_genotypes.tsv somatic_AF.tsv
python Data_Exploration/Genotype_AF.py ~/Wednesday/H730/noncoding/projects/Mirjam_Outputs_Pipeline/Genotype_Tables/ germline_genotypes_270824.tsv germline_AF.tsv

# Run allele frequency plot script in R
Rscript R_Scripts/Allele_Freq.R germline_AF.tsv somatic_AF.tsv


# 5. PSI

# Create a PSI summary table per AS event with Python
python Data_Exploration/Make_PSI_Table.py -f ~/Wednesday/H730/noncoding/projects/Mirjam_Outputs_Pipeline/PSI_Tables -o PSI_summary.tsv

# Make separate summary tables for each AS type event
python Data_Exploration/AStype_summary.py PSI_summary.tsv AA_summary.tsv AA
python Data_Exploration/AStype_summary.py PSI_summary.tsv AD_summary.tsv AD
python Data_Exploration/AStype_summary.py PSI_summary.tsv CE_summary.tsv CE
python Data_Exploration/AStype_summary.py PSI_summary.tsv IR_summary.tsv IR

#Make median tables
Rscript R_Scripts/PSI_Prep.R AA
Rscript R_Scripts/PSI_Prep.R AD
Rscript R_Scripts/PSI_Prep.R CE
Rscript R_Scripts/PSI_Prep.R IR

# Make four panel median plot
Rscript R_Scripts/PSI_Plot.R 

#This requires the summary files to be in the same folder.