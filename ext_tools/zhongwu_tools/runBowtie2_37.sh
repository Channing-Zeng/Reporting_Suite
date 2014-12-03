#!/bin/bash

REF=/ngs/reference_data/genomes/Hsapiens/GRCh37/bowtie2/GRCh37.fa

# Usage: runBowtie2_hg19.sh 1.fastq 2.fastq sample
/group/cancer_informatics/tools_resources/NGS/bin/bowtie2 -p 8 --sensitive-local --dovetail --quiet --rg-id $3 --rg "SM:$3" -x $REF -1 $1 -2 $2 > ${3}.bowtie2.sam
/group/cancer_informatics/tools_resources/NGS/bin/mapping_stat_Bowtie.pl -s $3 ${3}.bowtie2.sam > stat_mapping_${3}_bowtie2.txt
samtools view -Sb -F 0x4 ${3}.bowtie2.sam > ${3}.bowtie2.bam
samtools sort ${3}.bowtie2.bam ${3}.bowtie2.sorted
samtools index ${3}.bowtie2.sorted.bam

/opt/az/sun/java/jdk64_1.6.0_24/bin/java -Xmx8g -jar /opt/az/local/gatk/GenomeAnalysisTK-2013.2-2.5.2-0-g3ae1219/GenomeAnalysisTK.jar -T RealignerTargetCreator -I ${3}.bowtie2.sorted.bam -R $REF -o ${3}.bowtie2.intervals -allowPotentiallyMisencodedQuals
/opt/az/sun/java/jdk64_1.6.0_24/bin/java -Xmx8g -jar /opt/az/local/gatk/GenomeAnalysisTK-2013.2-2.5.2-0-g3ae1219/GenomeAnalysisTK.jar -T IndelRealigner -I ${3}.bowtie2.sorted.bam -R $REF -targetIntervals ${3}.bowtie2.intervals -allowPotentiallyMisencodedQuals --out ${3}.bowtie2.sorted.realign.bam

rm ${3}.bowtie2.sam ${3}.bowtie2.bam
