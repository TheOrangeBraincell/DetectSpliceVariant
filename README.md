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

To run the pipeline, pysam (v0.19.0) and samtools (v1.15.1) have to be installed via conda. 

## Running the Pipeline

To run the Pipeline, you require the Scripts/ folder and its contents, as well as the annotation files described above. You also require a list of sample names and a list of genes of interest as input. Both need to have one sample or gene respectively, per line.

To start a run, you run the Scripts/Run_Pipeline.py script as follows. Make sure to have your conda environment with the installations activated.

```
Initializes Pipeline run to find variants associated with alternative splicing.

optional arguments:
  -h, --help            show this help message and exit
  --samples SAMPLES, -s SAMPLES
                        file containing the names of all samples that the pipeline is to be run for.
                        One line per sample.
  --genes GENES, -g GENES
                        File containing all HGNC gene names, the pipeline should run for. One gene
                        name per line.
  --out OUT, -o OUT     Output directory, containing pipeline outputs in subfolders.
  --gencode GENCODE, -gc GENCODE
                        tsv file containing bed file information on annotated exons from GENCODE39 as
                        well as gene names. Instructions on how to find the tables, can be found on
                        github.
  --refseq REFSEQ, -rs REFSEQ
                        tsv file containing bed file information on annotated exons from RefSeq as
                        well as gene names. Instructions on how to find the tables, can be found on
                        github.
  --vcf VCF, -vcf VCF   Folder containing vcf files for all samples.
  --bam BAM, -bam BAM   Folder containing bam files for all samples.
  --gff GFF, -gff GFF   gff3 annotation file from gencode (here used v39).
  --cores CORES, -n CORES
                        how many cores youd like the pipeline to use.
```

The script will create the file tree needed for the pipeline, as long as it finds its prerequisites in the folder. Otherwise it will exit and ask for what it could not find.
It will start a pipeline run per gene in the input list of gene, using 3 cores per gene run. Thus it will run as many genes parallel as the number of cores given to it divided by 3.

