#!/bin/bash

# Usage: runBowtie2_mm10.sh 1.fastq 2.fastq sample
/group/cancer_informatics/tools_resources/NGS/bin/bowtie2 -p 8 --sensitive-local --dovetail --quiet --rg-id $3 --rg "SM:$3" -x /group/cancer_informatics/tools_resources/NGS/genomes/mm10/Sequence/BowtieIndex/bowtie2/genome -1 $1 -2 $2 > ${3}_bowtie2.sam
/group/cancer_informatics/tools_resources/NGS/bin/mapping_stat_Bowtie.pl -s $3 ${3}_bowtie2.sam > ${3}_stat_mapping.txt
samtools view -Sb -F 0x4 ${3}_bowtie2.sam > ${3}_bowtie2.bam
samtools sort ${3}_bowtie2.bam ${3}_bowtie2_sorted
samtools index ${3}_bowtie2_sorted.bam

/opt/az/sun/java/jdk64_1.6.0_24/bin/java -Xmx8g -jar /opt/az/local/gatk/GenomeAnalysisTK-2013.2-2.5.2-0-g3ae1219/GenomeAnalysisTK.jar -T RealignerTargetCreator -I ${3}_bowtie2_sorted.bam -R /group/cancer_informatics/tools_resources/NGS/genomes/mm10/Sequence/WholeGenomeFasta/genome.fa -o ${3}_bowtie2.intervals -allowPotentiallyMisencodedQuals
/opt/az/sun/java/jdk64_1.6.0_24/bin/java -Xmx8g -jar /opt/az/local/gatk/GenomeAnalysisTK-2013.2-2.5.2-0-g3ae1219/GenomeAnalysisTK.jar -T IndelRealigner -I ${3}_bowtie2_sorted.bam -R /group/cancer_informatics/tools_resources/NGS/genomes/mm10/Sequence/WholeGenomeFasta/genome.fa -targetIntervals ${3}_bowtie2.intervals -allowPotentiallyMisencodedQuals --out ${3}_bowtie2_sorted.realign.bam

rm ${3}_bowtie2.sam ${3}_bowtie2.bam
