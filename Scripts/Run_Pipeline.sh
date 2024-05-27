##22.09.23
# Mirjam Müller
# Subprocesses Pipeline

#To be run for every gene of the gene list. Core part of pipeline.
#Inputs: 1. gene, 2. output directory, 3. refseq, 4. gencode, 5. parameter string
start_pipeline=`date +%s`
#Initialize log file for gene.
echo "#Log file for ${1}" > ${2}Log_Files_genes/${1}_log.txt

#Start run of AS and variants parallel, as AS needs 2 cores and variants only needs 1.
python Scripts/Identify_AS.py -o ${2}AS_Events/${1}_AS_events.tsv -n $1 -g $4 -r $3 >temp_AS_${1}.txt &
python Scripts/vcf_location_table.py -s vcf_file_list.txt -o ${2}Variant_Locations/${1}_locations.tsv -n ${1} -r gene_ranges.bed >temp_Var_${1}.txt &
wait #So both processes are done before we proceed.

#Check if scripts ran successfully, if they didnt, exit Run_Pipeline.sh.
tail -n 1 temp_AS_${1}.txt | grep -q 'Run time' || { echo "Crash: Script Identify_AS.py terminated with an error message, stop pipeline run for gene $1"; exit 1; }
tail -n 1 temp_Var_${1}.txt | grep -q 'Run time' || { echo "Crash: Script vcf_location_table.py terminated with an error message, stop pipeline run for gene $1"; exit 1; }

#Add temp log files into the general gene log file!
cat temp_AS_${1}.txt >>${2}Log_Files_genes/${1}_log.txt
cat temp_Var_${1}.txt >>${2}Log_Files_genes/${1}_log.txt
#Delete temporary log files.
rm temp_AS_${1}.txt
rm temp_Var_${1}.txt 

# Now Psiscores, this one uses 3 cores as is
python Scripts/PsiScores.py -i ${2}AS_Events/${1}_AS_events.tsv -o ${2}PSI_Tables/${1}_PSI.tsv -s bam_file_list.txt -is "${5}" >> temp_PSI_${1}.txt &
#at the same time we can filter out the utr variant locations with the 1 remaining core. 
python Scripts/utr_variants.py -v ${2}Variant_Locations/${1}_locations.tsv -o ${2}Variant_Locations/${1}_locations_noutr.tsv -n $1 -g $4 -r $3 >> temp_UTR_${1}.txt &
#wait for both processes to finish
wait

#Check if scripts ran successfully, if they didnt, exit Run_Pipeline.sh
tail -n 1 temp_PSI_${1}.txt | grep -q 'Run time' || { echo "Crash: Script PsiScores.py terminated with an error message, stop pipeline run for gene $1"; exit 1; }
tail -n 1 temp_UTR_${1}.txt | grep -q 'Run time' || { echo "Crash: Script utr_variants.py terminated with an error message, stop pipeline run for gene $1"; exit 1; }

#Add temp log files to general gene log file
cat temp_PSI_${1}.txt >> ${2}Log_Files_genes/${1}_log.txt
cat temp_UTR_${1}.txt >> ${2}Log_Files_genes/${1}_log.txt
#delete temporary log files
rm temp_PSI_${1}.txt
rm temp_UTR_${1}.txt


#The first result exploration showed that only the genotypes for exon variants are reliable. So we keep only those.
python Scripts/Exon_Variants.py -v ${2}Variant_Locations/${1}_locations_noutr.tsv -g $4 -r $3 -o ${2}Exon_Variants/${1}_exonvar.tsv >> ${2}Log_Files_genes/${1}_log.txt

#Check if scripts ran successfully, if they didnt, exit Run_Pipeline.sh
tail -n 1 ${2}Log_Files_genes/${1}_log.txt | grep -q 'Run time' || { echo "Crash: Script Exon_Variants.py terminated with an error message, stop pipeline run for gene $1"; exit 1; }


#When that one is done we do genotypes.
#Bedtools multicov needs the processes to be split into 4, because it can only handle 1021 samples at a time (1021 bam files)
#That means we split the bam file list so that theres less than that number, and then run it parallel for those 4.
#Mind that this needs to be adjusted if we run it for more samples.
how_many=$((`wc -l < bam_file_list.txt` / 1021 +1)) #thats into how many files.
number_lines=$((`wc -l < bam_file_list.txt`/$how_many +1))
split -l $number_lines bam_file_list.txt ${1}_split_
#make bed file out of locations
echo ""> ${1}_variants.bed;

tail -n +2 ${2}Exon_Variants/${1}_exonvar.tsv| awk -F'\t' '{split($1, arr, "_"); start=arr[2]; end=start+1; print arr[1] "\t" start "\t" end "\t" $1}' > ${1}_variants.bed
for file in ${1}_split_*; do Scripts/bedtools.sh $file ${1}_variants.bed temp_${file}.tsv > log_${file}.tsv & done
wait
#Add header to read depth output
cat ${1}_split_* | tr '\n' '\t' | paste > ${2}Read_Depth/${1}_read_depth.tsv
#add locations
cut -f1-6 temp_${1}_split_aa.tsv > ${1}_locations.txt
for file in temp_${1}_split_*; do cut -f7- $file > ${1}_cut_${file}; done
paste ${1}_cut_* > ${1}_counts.txt
paste ${1}_locations.txt ${1}_counts.txt >> ${2}Read_Depth/${1}_read_depth.tsv

wait 
#remove location and counts file
rm ${1}_locations.txt
rm ${1}_counts.txt
#remove split files.
rm ${1}_split_*
#remove bed file
rm ${1}_variants.bed;
#remove temp files
rm temp_${1}_split_*
rm ${1}_cut_*
#remove log files
rm log_${1}_split_*

#When that one is done, we need to read off the read depth table to generate the genotypes.
python Scripts/genotype.py -v ${2}Exon_Variants/${1}_exonvar.tsv -r ${2}Read_Depth/${1}_read_depth.tsv -o ${2}Genotype_Tables/${1}_genotypes.tsv -g $1 >>${2}Log_Files_genes/${1}_log.txt

#Check if scripts ran successfully, if they didnt, exit Run_Pipeline.sh
tail -n 1 ${2}Log_Files_genes/${1}_log.txt | grep -q 'Run time' || { echo "Crash: Script genotype.py terminated with an error message, stop pipeline run for gene $1"; exit 1; }


end_pipeline=`date +%s`
runtime=$((end_pipeline-start_pipeline))
echo "Runtime for gene ${1} is ${runtime} seconds" >> ${2}Log_Files_genes/${1}_log.txt
