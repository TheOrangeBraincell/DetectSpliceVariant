## 19.09.23
# Mirjam Karlsson-Müller
# Script to create read depth table for variant locations.

#Input is the gene name, used to create the output file!

echo "#Read Depth Table"> ${1}_read_depth.tsv
cat $1_locations.tsv | grep -v "^#" | grep -v "^Location"| cut -f1 | while read var_ID; do coordinate=$( echo $var_ID | awk -F '_' '{print $1":"$2+1"-"$2+1}'); cat bam_file_list.txt | while read file; do prefix=$(echo ${file%/l*}); id=$(echo ${prefix#*Data/}); count=$(samtools depth -a -s -r $coordinate $file | cut -f3); printf $var_ID"\t"$id"\t"$count"\n">>${1}_read_depth.tsv; done; done

