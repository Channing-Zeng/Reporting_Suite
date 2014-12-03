#!/bin/bash

# Usage: runBWA074_hg19_ion.sh 1.fastq sample target_bed
PATH=$PATH:/group/cancer_informatics/tools_resources/NGS/bin
sample=$2

/users/kdld047/local/bwa-0.7.4/bwa mem -t 6 -P -c 20 -R "@RG\tID:${sample}\tSM:${sample}" /users/kdld047/work/NGS/NGS/genomes/hg19/Sequence/BWAIndex/genome.fa $1 > ${sample}.sam 2> bwa.err 
samtools view -Sb -F 0x4 ${sample}.sam > ${sample}.bam
samtools sort ${sample}.bam ${sample}_sorted
samtools index ${sample}_sorted.bam

/group/cancer_informatics/tools_resources/NGS/bin/mapping_stat_BWA.pl ${sample} ${sample}.sam > stat_sam_${sample}.txt
rm ${sample}.sam
rm ${sample}.bam

/opt/az/sun/java/jdk64_1.6.0_24/bin/java -Xmx8g -jar /opt/az/local/gatk/GenomeAnalysisTK-2013.2-2.5.2-0-g3ae1219/GenomeAnalysisTK.jar -T RealignerTargetCreator -I ${sample}_sorted.bam -R /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Sequence/WholeGenomeFasta/genome.fa -o ${sample}.intervals -allowPotentiallyMisencodedQuals
/opt/az/sun/java/jdk64_1.6.0_24/bin/java -Xmx8g -jar /opt/az/local/gatk/GenomeAnalysisTK-2013.2-2.5.2-0-g3ae1219/GenomeAnalysisTK.jar -T IndelRealigner -I ${sample}_sorted.bam -R /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Sequence/WholeGenomeFasta/genome.fa -targetIntervals ${sample}.intervals -allowPotentiallyMisencodedQuals --out ${sample}_sorted.realign.bam
samtools index ${sample}_sorted.realign.bam
if [ $3 ]
    then
	/group/cancer_informatics/tools_resources/NGS/bin/checkVar.pl -c 1 -s 2 -e 3 -S 2 -E 3 -g 4 -x 0 -Q 1 -f 0.0075 -N ${sample} -b ${sample}_sorted.realign.bam $3 > ${sample}_vars.txt
	/group/cancer_informatics/tools_resources/NGS/bin/teststrandbias.R ${sample}_vars.txt > ${sample}_vars.txt.t
	mv ${sample}_vars.txt.t ${sample}_vars.txt
	#/group/cancer_informatics/tools_resources/NGS/bin/checkCov.pl -c 1 -s 2 -e 3 -S 2 -E 3 -g 4 -N ${sample} -b ${sample}_sorted.realign.bam -d 1:10:50:100:500:1000:5000:10000:50000 $3 > ${sample}_cov.txt
	/group/cancer_informatics/tools_resources/NGS/bin/checkCov.pl -c 1 -s 2 -e 3 -S 2 -E 3 -g 4 -N ${sample} -b ${sample}_sorted.realign.bam -d 1:5:10:25:50:100:500:1000:5000:10000:50000 $3 > ${sample}_cov.txt
	/group/cancer_informatics/tools_resources/NGS/bin/var2vcf.pl ${sample}_vars.txt > ${sample}_vars.vcf
	/opt/az/oracle/java/jdk1.7.0_11/bin/java -Xmx4g -jar /group/cancer_informatics/tools_resources/NGS/snpEff/snpEff.jar eff -c /group/cancer_informatics/tools_resources/NGS/snpEff/snpEff.config -d -v -canon hg19 ${sample}_vars.vcf > ${sample}_vars.eff.vcf
	/opt/az/oracle/java/jdk1.7.0_11/bin/java -Xmx4g -jar /group/cancer_informatics/tools_resources/NGS/snpEff/SnpSift.jar annotate -v /ngs/cancer_informatics/GenomeData/human/dbsnp_latest.vcf ${sample}_vars.eff.vcf > ${sample}_vars.eff.dbsnp.vcf
	/opt/az/oracle/java/jdk1.7.0_11/bin/java -Xmx4g -jar /group/cancer_informatics/tools_resources/NGS/snpEff/SnpSift.jar annotate -v /ngs/cancer_informatics/GenomeData/human/CosmicCodingMuts_latest.vcf ${sample}_vars.eff.dbsnp.vcf > ${sample}_vars.eff.dbsnp.cosmic.vcf
fi

rm ${sample}_sorted.bam
touch runBWA074_hg19.done

#samtools view ${sample}.bam | mapping_stat_BWA.pl ${sample}_bam > stat_bam.txt
#samtools view ${sample}_sorted.bam | mapping_stat_BWA.pl ${sample}_sorted > stat_sorted_bam.txt

# Remove duplicates by Picard.  Time: 64.63min, Mem: 22G (21,250,637,824)
#module load java
#java -jar ~/local/picard-tools-1.70/MarkDuplicates.jar INPUT=${sample}_sorted.bam OUTPUT=${sample}_sorted_NoDup.bam COMMENT="Remove Duplicates Using Picard" REMOVE_DUPLICATES=true METRICS_FILE=duplicates.txt ASSUME_SORTED=true
#samtools view ${sample}_sorted_NoDup.bam | mapping_stat_BWA.pl ${sample}_NoDup > stat_NoDup.txt
