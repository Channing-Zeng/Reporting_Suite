#!/bin/bash

# Usage: runBWA_mouse.sh 1.fastq 2.fastq sample

~/local/bwa-0.6.2/bwa aln -t 16 -k 2 -l 30 -n 4 -R 5 /users/kdld047/work/NGS/NGS/genomes/mm10_complete/BWAIndex/mm10_complete.fa $1 > ${3}_1.sai 2> 1err
~/local/bwa-0.6.2/bwa aln -t 16 -k 2 -l 30 -n 4 -R 5 /users/kdld047/work/NGS/NGS/genomes/mm10_complete/BWAIndex/mm10_complete.fa $2 > ${3}_2.sai 2> 2err
~/local/bwa-0.6.2/bwa sampe -P -r "@RG\tID:$3\tSM:$3" /users/kdld047/work/NGS/NGS/genomes/mm10_complete/BWAIndex/mm10_complete.fa ${3}_1.sai ${3}_2.sai $1 $2 > $3.sam 2> sampe.err 
samtools view -Sb -F 0x4 $3.sam > $3.bam
samtools sort $3.bam ${3}_sorted
samtools index ${3}_sorted.bam
samtools view ${3}_sorted.bam | mapping_stat_BWA.pl $3 > stat.txt

module load java
# Remove duplicates by Picard.  Time: 64.63min, Mem: 22G (21,250,637,824)
java -jar ~/local/picard-tools-1.70/MarkDuplicates.jar INPUT=${3}_sorted.bam OUTPUT=${3}_sorted_NoDup.bam COMMENT="Remove Duplicates Using Picard" REMOVE_DUPLICATES=true METRICS_FILE=duplicates.txt ASSUME_SORTED=true
samtools view ${3}_sorted_NoDup.bam | mapping_stat_BWA.pl ${3}_NoDup > stat_NoDup.txt
