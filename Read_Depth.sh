#!/bin/bash -l
## 19.09.23
# Mirjam Karlsson-Müller
# Script to create read depth table for variant locations.

#Input is the gene name, used to create the output file!

start=`date +%s`
num_procs=3 #because thats how it is for the rest of the subparts of the pipeline
num_jobs="\j"
echo "#Read Depth Table"> ${1}_read_depth.tsv;

run_script() {
    echo $1;
    echo $var_ID;
    coordinate=$( echo $var_ID | awk -F '_' '{print $1":"$2+1"-"$2+1}');
    cat bam_file_list.txt | while read file; do prefix=$(echo ${file%/l*}); id=$(echo ${prefix#*BAM/}); count=$(samtools depth -a -s -r $coordinate $file | cut -f3); printf $var_ID"\t"$id"\t"$count"\n">>${1}_read_depth.tsv; done;
}

cat $1_locations.tsv | grep -v "^#" | grep -v "^Location"| cut -f1 | while read var_ID;
do while (( ${num_jobs@P} >= num_procs));
do wait -n;
done;
run_script $1 &
done;

end=`date +%s`
runtime=$((end-start))
echo $runtime