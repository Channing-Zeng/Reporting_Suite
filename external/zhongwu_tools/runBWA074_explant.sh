#!/bin/bash

# Usage: runBowtie2_explant.sh 1.fastq 2.fastq sample bed
fastq1=$1
fastq2=$2
sample=$3
bed=$4
#/group/cancer_informatics/tools_resources/NGS/bin/bowtie2 -p 10 --sensitive-local --dovetail --quiet --rg-id $sample --rg "SM:$sample" -x /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Sequence/Bowtie2Index/genome -1 $fastq1 -2 $fastq2 > ${sample}_hg19.sam
/users/kdld047/local/bwa-0.7.4/bwa mem -t 8 -P -c 20 -R "@RG\tID:$sample\tSM:$sample" /users/kdld047/work/NGS/NGS/genomes/hg19/Sequence/BWAIndex/genome.fa $fastq1 $fastq2 > ${sample}_hg19.sam 2> bwa_hg19.err
/group/cancer_informatics/tools_resources/NGS/bin/mapping_stat_BWA.pl $sample ${sample}_hg19.sam > stat_mapping_hg19_$sample.txt
samtools view -Sb -F 0x4 ${sample}_hg19.sam > ${sample}_hg19.bam
samtools sort ${sample}_hg19.bam ${sample}_hg19_sorted
samtools index ${sample}_hg19_sorted.bam

/opt/az/sun/java/jdk64_1.6.0_24/bin/java -Xmx8g -jar /opt/az/local/gatk/GenomeAnalysisTK-2013.2-2.5.2-0-g3ae1219/GenomeAnalysisTK.jar -T RealignerTargetCreator -I ${sample}_hg19_sorted.bam -R /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Sequence/WholeGenomeFasta/genome.fa -o ${sample}_hg19_sorted.intervals -allowPotentiallyMisencodedQuals
/opt/az/sun/java/jdk64_1.6.0_24/bin/java -Xmx8g -jar /opt/az/local/gatk/GenomeAnalysisTK-2013.2-2.5.2-0-g3ae1219/GenomeAnalysisTK.jar -T IndelRealigner -I ${sample}_hg19_sorted.bam -R /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Sequence/WholeGenomeFasta/genome.fa -targetIntervals ${sample}_hg19_sorted.intervals -allowPotentiallyMisencodedQuals --out ${sample}_hg19_sorted.realign.bam

rm ${sample}_hg19.sam ${sample}_hg19.bam ${sample}_hg19_sorted.ba*

#/group/cancer_informatics/tools_resources/NGS/bin/bowtie2 -p 10 --sensitive-local --dovetail --quiet --rg-id $sample --rg "SM:$sample" -x /group/cancer_informatics/tools_resources/NGS/genomes/mm10/Sequence/BowtieIndex/bowtie2/genome -1 $fastq1 -2 $fastq2 > ${sample}_mm10.sam
/users/kdld047/local/bwa-0.7.4/bwa mem -t 8 -P -c 20 -R "@RG\tID:$sample\tSM:$sample" /users/kdld047/work/NGS/NGS/genomes/mm10/Sequence/BWAIndex/genome.fa $fastq1 $fastq2 > ${sample}_mm10.sam 2> bwa_mm10.err
/group/cancer_informatics/tools_resources/NGS/bin/mapping_stat_BWA.pl $sample ${sample}_mm10.sam > stat_mapping_mm10_$sample.txt
samtools view -Sb -F 0x4 ${sample}_mm10.sam > ${sample}_mm10.bam
samtools sort ${sample}_mm10.bam ${sample}_mm10_sorted
samtools index ${sample}_mm10_sorted.bam

/opt/az/sun/java/jdk64_1.6.0_24/bin/java -Xmx8g -jar /opt/az/local/gatk/GenomeAnalysisTK-2013.2-2.5.2-0-g3ae1219/GenomeAnalysisTK.jar -T RealignerTargetCreator -I ${sample}_mm10_sorted.bam -R /group/cancer_informatics/tools_resources/NGS/genomes/mm10/Sequence/WholeGenomeFasta/genome.fa -o ${sample}_mm10_sorted.intervals -allowPotentiallyMisencodedQuals
/opt/az/sun/java/jdk64_1.6.0_24/bin/java -Xmx8g -jar /opt/az/local/gatk/GenomeAnalysisTK-2013.2-2.5.2-0-g3ae1219/GenomeAnalysisTK.jar -T IndelRealigner -I ${sample}_mm10_sorted.bam -R /group/cancer_informatics/tools_resources/NGS/genomes/mm10/Sequence/WholeGenomeFasta/genome.fa -targetIntervals ${sample}_mm10_sorted.intervals -allowPotentiallyMisencodedQuals --out ${sample}_mm10_sorted.realign.bam

/group/cancer_informatics/tools_resources/NGS/bin/compareSAM_AS.pl -s _AS ${sample}_hg19_sorted.realign.bam ${sample}_mm10_sorted.realign.bam

rm ${sample}_mm10.sam ${sample}_mm10.bam ${sample}_mm10_sorted.ba*

samtools view -h ${sample}_hg19_sorted.realign.bam | pickLine -v reads_sample2_AS.txt | samtools view -Sb - > ${sample}_hg19_sorted.realign.disamb.bam
samtools index ${sample}_hg19_sorted.realign.disamb.bam
rm ${sample}_hg19_sorted.bam*
rm ${sample}_mm10_sorted.bam*

rawsample=raw_$sample
if [ $bed ]
    then
        /group/cancer_informatics/tools_resources/NGS/bin/checkCov.pl -c 1 -s 2 -e 3 -S 2 -E 3 -g 4 -N $sample -b ${sample}_hg19_sorted.realign.disamb.bam -d 1:5:10:25:50:100:500:1000:5000:10000:50000 $bed > ${sample}_cov.txt
        /group/cancer_informatics/tools_resources/NGS/bin/checkVar.pl -c 1 -s 2 -e 3 -S 2 -E 3 -g 4 -x 0 -f 0.0075 -Q 10 -N $sample -b ${sample}_hg19_sorted.realign.disamb.bam $bed > ${sample}_vars.txt
        /group/cancer_informatics/tools_resources/NGS/bin/teststrandbias.R ${sample}_vars.txt > ${sample}_vars.txt.t
        mv ${sample}_vars.txt.t ${sample}_vars.txt
        /group/cancer_informatics/tools_resources/NGS/bin/var2vcf.pl ${sample}_vars.txt > ${sample}_vars.vcf
        /opt/az/oracle/java/jdk1.7.0_11/bin/java -Xmx4g -jar /group/cancer_informatics/tools_resources/NGS/snpEff/snpEff.jar eff -c /group/cancer_informatics/tools_resources/NGS/snpEff/snpEff.config -d -v -canon hg19 ${sample}_vars.vcf > ${sample}_vars.eff.vcf
        /opt/az/oracle/java/jdk1.7.0_11/bin/java -Xmx4g -jar /group/cancer_informatics/tools_resources/NGS/snpEff/SnpSift.jar annotate -v /ngs/cancer_informatics/GenomeData/human/dbsnp137_00_all.vcf ${sample}_vars.eff.vcf > ${sample}_vars.eff.dbsnp.vcf
        /opt/az/oracle/java/jdk1.7.0_11/bin/java -Xmx4g -jar /group/cancer_informatics/tools_resources/NGS/snpEff/SnpSift.jar annotate -v /ngs/cancer_informatics/GenomeData/human/CosmicCodingMuts_latest.vcf ${sample}_vars.eff.dbsnp.vcf > ${sample}_vars.eff.dbsnp.cosmic.vcf

        /group/cancer_informatics/tools_resources/NGS/bin/checkCov.pl -c 1 -s 2 -e 3 -S 2 -E 3 -g 4 -N $rawsample -b ${sample}_hg19_sorted.realign.bam -d 1:5:10:25:50:100:500:1000:5000:10000:50000 $bed > ${rawsample}_cov.txt
        /group/cancer_informatics/tools_resources/NGS/bin/checkVar.pl -c 1 -s 2 -e 3 -S 2 -E 3 -g 4 -x 0 -f 0.0075 -Q 10 -N $rawsample -b ${sample}_hg19_sorted.realign.bam $bed > ${rawsample}_vars.txt
        /group/cancer_informatics/tools_resources/NGS/bin/teststrandbias.R ${rawsample}_vars.txt > ${rawsample}_vars.txt.t
        mv ${rawsample}_vars.txt.t ${rawsample}_vars.txt
        /group/cancer_informatics/tools_resources/NGS/bin/var2vcf.pl ${rawsample}_vars.txt > ${rawsample}_vars.vcf
        /opt/az/oracle/java/jdk1.7.0_11/bin/java -Xmx4g -jar /group/cancer_informatics/tools_resources/NGS/snpEff/snpEff.jar eff -c /group/cancer_informatics/tools_resources/NGS/snpEff/snpEff.config -d -v -canon hg19 ${rawsample}_vars.vcf > ${rawsample}_vars.eff.vcf
        /opt/az/oracle/java/jdk1.7.0_11/bin/java -Xmx4g -jar /group/cancer_informatics/tools_resources/NGS/snpEff/SnpSift.jar annotate -v /ngs/cancer_informatics/GenomeData/human/dbsnp137_00_all.vcf ${rawsample}_vars.eff.vcf > ${rawsample}_vars.eff.dbsnp.vcf
        /opt/az/oracle/java/jdk1.7.0_11/bin/java -Xmx4g -jar /group/cancer_informatics/tools_resources/NGS/snpEff/SnpSift.jar annotate -v /ngs/cancer_informatics/GenomeData/human/CosmicCodingMuts_latest.vcf ${rawsample}_vars.eff.dbsnp.vcf > ${rawsample}_vars.eff.dbsnp.cosmic.vcf
fi

