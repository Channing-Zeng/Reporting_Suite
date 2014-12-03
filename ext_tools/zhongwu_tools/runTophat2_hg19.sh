#!/bin/bash

# Usage $0 sample fastq1 fastq2
#/users/kdld047/work/NGS/NGS/tophat-2.0.3.Linux_x86_64/tophat -r 20 --mate-std-dev 35 -p 12 -G ~/work/NGS/NGS/genomes/mm10/Annotation/Genes/genes.gtf -o tophat2_out_mouse --rg-id $1 --rg-sample $1 ~/work/NGS/NGS/genomes/mm10/Sequence/BowtieIndex/bowtie2/genome $2 $3 2> mouse_tophat2.err

#export PATH=/group/cancer_informatics/tools_resources/NGS/bin:$PATH

#/users/kdld047/work/NGS/NGS/tophat-2.0.10.Linux_x86_64/tophat -r 25 --mate-std-dev 40 -p 8 -G /ngs/reference_data/genomes/Hsapiens/hg19/rnaseq-2013-09-25/genes.gtf -o tophat2.0.4_out_hg19 --transcriptome-index /ngs/reference_data/genomes/Hsapiens/hg19/rnaseq/tophat/hg19_transcriptome --rg-id $1 --rg-sample $1 /ngs/reference_data/genomes/Hsapiens/hg19/bowtie2/hg19 $2 $3 2> hg19_tophat2.err
/users/kdld047/work/NGS/NGS/tophat-2.0.10.Linux_x86_64/tophat -r 25 --mate-std-dev 40 -p 8 -G ~/work/NGS/NGS/genomes/hg19/Annotation/Genes/genes.gtf -o . --transcriptome-index ~/work/NGS/NGS/genomes/hg19/transcriptome --rg-id $1 --rg-sample $1 ~/work/NGS/NGS/genomes/hg19/Sequence/Bowtie2Index/genome $2 $3 2> human_tophat2.err 

mv accepted_hits.bam ${1}_tophat2.0.10.bam
samtools index ${1}_tophat2.0.10.bam

