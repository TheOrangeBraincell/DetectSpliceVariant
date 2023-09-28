#!/bin/bash

start=`date +%s`

# $1 is the gene that we are doing it for, $2 the bam file list, $3 the output folder.
#make the bed file.
echo "" > ${1}_variants.bed;
cat ${3}Variant_Locations/${1}_locations.tsv | grep -v "^#" | grep -v "^Location"| cut -f1 | while read var_ID; do echo $var_ID | awk -F '_' '{print $1"\t"$2"\t"$2+1"\t"$1"_"$2"_"$3"_"$4"\t.\t+"}' >> ${1}_variants.bed; done
#cat ${1}_locations.tsv | grep -v "^#" | grep -v "^Location"| cut -f1 | while read var_ID; do echo $var_ID | awk -F '_' '{print $1"\t"$2"\t"$2+1"\t"$1"_"$2"_"$3"_"$4"\t.\t+"}' >> ${1}_variants.bed; done

#split into 3 files, for "parallelization"
split -l $(($(wc -l <${1}_variants.bed)/3)) ${1}_variants.bed
#xaa, xab, xac

#Initialize output file.
echo "#Read Depth for ${1}" > ${3}Read_Depth/${1}_read_depth.tsv;
#echo "#Read Depth for ${1}" > Read_Depth/${1}_read_depth.tsv;

#For each split file, do the loop.
#cat ${2} | while read file; do prefix=$(echo ${file%/l*});
#id=$(echo ${prefix#*Data/}); samtools bedcov xaa $file | cut -f4,7 | while read var_id count; do printf "$id\t$var_id\t$count\n" >>Read_Depth/${1}_read_depth.tsv; done ; done &

#cat ${2}  | while read file; do prefix=$(echo ${file%/l*});
#id=$(echo ${prefix#*Data/}); samtools bedcov xab $file | cut -f4,7 | while read var_id count; do printf "$id\t$var_id\t$count\n" >>Read_Depth/${1}_read_depth.tsv; done ; done &

#cat ${2}  | while read file; do prefix=$(echo ${file%/l*});
#id=$(echo ${prefix#*Data/}); samtools bedcov xac $file | cut -f4,7 | while read var_id count; do printf "$id\t$var_id\t$count\n" >>Read_Depth/${1}_read_depth.tsv; done ; done &


#Iterate through bam files in file  list (input 2). The samtools bedcov output doesnt show the sample name. So add that.
cat ${2} | while read file; do prefix=$(echo ${file%/l*});
id=$(echo ${prefix#*Data/}); samtools bedcov xaa $file | cut -f4,7 | while read output; do printf "$id\t$output\n" >>${3}Read_Depth/${1}_read_depth.tsv; done ; done &

cat ${2}  | while read file; do prefix=$(echo ${file%/l*});
id=$(echo ${prefix#*Data/}); samtools bedcov xab $file | cut -f4,7 | while read output; do printf "$id\t$output\n" >>${3}Read_Depth/${1}_read_depth.tsv; done ; done &

cat ${2}  | while read file; do prefix=$(echo ${file%/l*});
id=$(echo ${prefix#*Data/}); samtools bedcov xac $file | cut -f4,7 | while read output; do printf "$id\t$output\n" >>${3}Read_Depth/${1}_read_depth.tsv; done ; done &

wait &&
rm ${1}_variants.bed
rm xaa
rm xab
rm xac

end=`date +%s`
runtime=$((end-start))
echo $runtime