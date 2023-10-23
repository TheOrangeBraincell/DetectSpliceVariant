##22.09.23
# Mirjam Karlsson-Müller
# Subprocesses Pipeline

#To be run for every gene of the gene list. Core part of pipeline.
#Inputs: 1. gene, 2. output directory, 3. refseq, 4. gencode, 5. parameter string

#Initialize log file for gene.
echo "#Log file for ${1}" > Log_Files_genes/${1}_log.txt

#Start run of AS and variants parallel, as AS needs 2 cores and variants only needs 1.
python Scripts/Identify_AS.py -o ${2}AS_Events/${1}_AS_events.tsv -n $1 -g $4 -r $3 >temp_AS_${1}.txt &
python Scripts/vcf_location_table.py -s vcf_file_list.txt -o ${2}Variant_Locations/${1}_locations.tsv -n ${1} -r gene_ranges.bed >temp_Var_${1}.txt &
wait #So both processes are done before we proceed.

#Add temp log files into the general gene log file!
cat temp_AS_${1}.txt >>Log_Files_genes/${1}_log.txt
cat temp_Var_${1}.txt >>Log_Files_genes/${1}_log.txt
#Delete temporary log files.
rm temp_AS_${1}.txt
rm temp_Var_${1}.txt 

# Now Psiscores, this one uses 3 cores as is
python Scripts/PsiScores.py -i ${2}AS_Events/${1}_AS_events.tsv -o ${2}PSI_Tables/${1}_PSI.tsv -s bam_file_list.txt -is "${5}" >> Log_Files_genes/${1}_log.txt
#at the same time we can filter out the utr variant locations with the 1 remaining core. 
python Scripts/utr_variants.py -v ${2}Variant_Locations/${1}_locations.tsv -o ${2}Variant_Locations/${1}_locations_noutr.tsv -n $1 -g $4 -r $3
#wait for both processes to finish
wait

#When that one is done we do genotypes.
#Bedtools multicov needs the processes to be split into 4, because it can only handle 1021 samples at a time (1021 bam files)
#That means we split the bam file list so that theres less than that number, and then run it parallel for those 4.
#Mind that this needs to be adjusted if we run it for more samples.
number_lines=$((`wc -l < bam_file_list.txt` / 1021 +1))
split -l $number_lines bam_file_list.txt #this needs to be adjusted to any number of samples instead.
#make bed file out of locations
echo ""> ${1}_variants.bed;
cat ${2}Variant_Locations/${1}_locations_noutr.tsv | grep -v "^#" | grep -v "^Location"| cut -f1 | while read var_ID; do echo $var_ID | awk -F '_' '{print $1"\t"$2"\t"$2+1"\t"$1"_"$2"_"$3"_"$4"\t.\t+"}' >> ${1}_variants.bed; done

for file in xa*; do Scripts/bedtools.sh $file ${1}_variants.bed ${2}Read_Depth/${1}_read_depth.tsv > ${1}_rd_log.txt; done
wait 
#remove split files.
rm xa*
#remove bed file
rm ${1}_variants.bed;

#When that one is done, we need to read off the read depth table to generate the genotypes.
python Scripts/genotype.py -v ${2}Variant_Locations/${1}_locations_noutr.tsv -r ${2}Read_Depth/${1}_read_depth.tsv -o ${2}Genotype_Tables/${1}_genotypes.tsv -g $1


