#!/bin/bash

#module load bcbio-nextgen/0.7.6

/opt/az/local/bcbio-nextgen/stable/0.7.6/tooldir/bin/STAR --genomeDir /ngs/reference_data/genomes/Hsapiens/hg19/star --readFilesIn $1 $2 --readFilesCommand zcat --runThreadN 8 --outFileNamePrefix $3 --outReadsUnmapped Fastx --outFilterMultimapNmax 10 --outSAMstrandField intronMotif

samtools view -Sb -@ 15 ${3}Aligned.out.sam > $3.bam
samtools sort -@ 15 -m 2G $3.bam ${3}_sorted
samtools index ${3}_sorted.bam

rm ${3}Aligned.out.sam $3.bam
