#!/bin/bash

#input one = bam file list, input two = bed file, 3 = output file
start=`date +%s`

files=$(cat $1 | tr '\n' ' ')

echo $files

bedtools multicov -p -bams $files -bed $2 > $3

end=`date +%s`
runtime=$((end-start))
echo "Finished bedtools"
echo $runtime