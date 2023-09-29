# Detect Splice Variants
*Mirjam Karlsson-Müller*

This is a pipeline develped to analyze the effect of genetic variants in breast cancer patients on mRNA splicing, based on the patients RNA Sequencing data.

#### Abbrevations
* AS = Alternative Splicing
* PSI = Percent Spliced In
* CE = Casette Exon
* IR = Intron Retention
* AD = Alternative Donor
* AA = Alternative Acceptor

## Steps

1. Using GENCODE and RefSeq annotation files, AS events are identified. 
2. By parsing all samples variant calling files, all variant locations passing filtering in at least one sample are collected.
3. Scoring AS events for each sample, based on the aligned reads of the sample.
4. Determining genotypes at variant locations for samples without a variant based on their read depth at that position.


The RNA Sequencing data is part of the Sweden Cancerome Analysis Network - Breast (SCAN-B). We are using their alignmentand variant calling files. Most of the pipeline is written in python (v.3.9.7), added to with shell scripts (v.5.0.16) and to be run in a linux console environment.

Detailed information with Screenshots on how to acquire the same RefSeq and GENCODE annotation files, can be found in the "How_To_Get_Database_Tables" folder. Besides those files, the gene annotation in gff format has also been used to determine the average insert size among samples.

## Installations

To run the pipeline, pysam (v.0.19.0) and samtools (v.0.17) have to be installed via conda. 


