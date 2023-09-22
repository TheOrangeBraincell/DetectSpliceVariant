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

# Now Psiscores, this one uses 3 cores as is, so the genotypes have to wait.
python Scripts/PsiScores.py -i ${2}AS_Events/${1}_AS_events.tsv -o ${2}PSI_Tables/${1}_PSI.tsv -s bam_file_list.txt -is $5 >> Log_Files_genes/${1}_log.txt

#When that one is done we do genotypes.
#First Read depth script to make read depth tables.
./Scripts/Read_Depth.sh $1 ${2}Read_Depth_Tables/${1}_read_depth.tsv
#When that one is done, we need to read off the read depth table to generate the genotypes.
python genotype.py


