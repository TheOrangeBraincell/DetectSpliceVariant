#!/bin/bash

# $1 is the gene that we are doing it for, $2 the bam file list
#make the bed file.
echo "" > ${1}_variants.bed;
cat ${1}_locations.tsv | grep -v "^#" | grep -v "^Location"| cut -f1 | while read var_ID; do echo $var_ID | awk -F '_' '{print $1"\t"$2"\t"$2+1"\t"$1"_"$2"_"$3"_"$4"\t.\t+"}' >> ${1}_variants.bed; done

#Initialize output file.
echo "#Read Depth for ${1}" > Read_Depth/${1}_read_depth.tsv;
cat $2 | while read file; do prefix=$(echo ${file%/l*});
id=$(echo ${prefix#*BAM/}); echo \#${id} >>Read_Depth/${1}_read_depth.tsv; samtools bedcov -j ${1}_variants.bed $file >>Read_Depth/${1}_read_depth.tsv; done

# 