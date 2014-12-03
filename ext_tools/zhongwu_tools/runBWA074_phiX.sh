#!/bin/bash

# Usage: runBWA.sh 1.fastq 2.fastq sample

/users/kdld047/local/bwa-0.7.4/bwa mem -t 8 -P -c 20 -R "@RG\tID:$3\tSM:$3" /group/cancer_informatics/tools_resources/NGS/genomes/PhiX/BWAIndex/PhiX.fa $1 $2 > $3.sam 2> bwa.err 
samtools view -Sb -F 0x4 $3.sam > $3.bam
samtools sort $3.bam ${3}_sorted
samtools index ${3}_sorted.bam

mapping_stat_BWA.pl $3 $3.sam > stat_sam_$3.txt
#rm $3.sam
rm $3.bam

