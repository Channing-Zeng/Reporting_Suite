#!/bin/bash

# Usage: runBWA074_hg19.sh 1.fastq 2.fastq sample target_bed
PATH=$PATH:/group/cancer_informatics/tools_resources/NGS/bin
REF=/ngs/reference_data/genomes/Hsapiens/GRCh37/bwa/GRCh37.fa

/users/kdld047/local/bwa-0.7.4/bwa mem -t 8 -P -c 20 -R "@RG\tID:$3\tSM:$3" $REF $1 $2 > $3.sam 2> bwa.err 
samtools view -Sb -F 0x4 $3.sam > $3.bam
samtools sort $3.bam ${3}.sorted
samtools index ${3}.sorted.bam

/group/cancer_informatics/tools_resources/NGS/bin/mapping_stat_BWA.pl $3 $3.sam > stat_sam_$3.txt
rm $3.sam
rm $3.bam*

/opt/az/sun/java/jdk64_1.6.0_24/bin/java -Xmx8g -jar /opt/az/local/gatk/GenomeAnalysisTK-2013.2-2.5.2-0-g3ae1219/GenomeAnalysisTK.jar -T RealignerTargetCreator -I ${3}.sorted.bam -R $REF -o ${3}.intervals -allowPotentiallyMisencodedQuals
/opt/az/sun/java/jdk64_1.6.0_24/bin/java -Xmx8g -jar /opt/az/local/gatk/GenomeAnalysisTK-2013.2-2.5.2-0-g3ae1219/GenomeAnalysisTK.jar -T IndelRealigner -I ${3}.sorted.bam -R $REF -targetIntervals ${3}.intervals -allowPotentiallyMisencodedQuals --out ${3}.sorted.realign.bam --maxReadsForRealignment 40000
samtools index ${3}.sorted.realign.bam
if [ $4 ]
    then
	/group/cancer_informatics/tools_resources/NGS/bin/checkVar.pl -c 1 -s 2 -e 3 -S 2 -E 3 -g 4 -x 0 -Q 1 -f 0.004 -N $3 -b ${3}.sorted.realign.bam $4 > ${3}_vars.txt
	/group/cancer_informatics/tools_resources/NGS/bin/teststrandbias.R ${3}_vars.txt > ${3}_vars.txt.t
	mv ${3}_vars.txt.t ${3}_vars.txt
	#/group/cancer_informatics/tools_resources/NGS/bin/checkCov.pl -c 1 -s 2 -e 3 -S 2 -E 3 -g 4 -N $3 -b ${3}.sorted.realign.bam -d 1:10:50:100:500:1000:5000:10000:50000 $4 > ${3}_cov.txt
	/group/cancer_informatics/tools_resources/NGS/bin/checkCov.pl -c 1 -s 2 -e 3 -S 2 -E 3 -g 4 -N $3 -b ${3}.sorted.realign.bam -d 1:5:10:25:50:100:500:1000:5000:10000:50000 $4 > ${3}_cov.txt
	/group/cancer_informatics/tools_resources/NGS/bin/var2vcf.pl ${3}_vars.txt > ${3}_vars.vcf
	/opt/az/oracle/java/jdk1.7.0_11/bin/java -Xmx4g -jar /group/cancer_informatics/tools_resources/NGS/snpEff/snpEff.jar eff -c /group/cancer_informatics/tools_resources/NGS/snpEff/snpEff.config -d -v -canon hg19 ${3}_vars.vcf > ${3}_vars.eff.vcf
	/opt/az/oracle/java/jdk1.7.0_11/bin/java -Xmx4g -jar /group/cancer_informatics/tools_resources/NGS/snpEff/SnpSift.jar annotate -v /ngs/cancer_informatics/GenomeData/human/dbsnp_latest.vcf ${3}_vars.eff.vcf > ${3}_vars.eff.dbsnp.vcf
	/opt/az/oracle/java/jdk1.7.0_11/bin/java -Xmx4g -jar /group/cancer_informatics/tools_resources/NGS/snpEff/SnpSift.jar annotate -v /ngs/cancer_informatics/GenomeData/human/CosmicCodingMuts_latest.vcf ${3}_vars.eff.dbsnp.vcf > ${3}_vars.eff.dbsnp.cosmic.vcf
fi

rm ${3}.sorted.bam*
touch runBWA074_hg19.done

#samtools view $3.bam | mapping_stat_BWA.pl ${3}_bam > stat_bam.txt
#samtools view ${3}.sorted.bam | mapping_stat_BWA.pl ${3}.sorted > stat.sorted_bam.txt

# Remove duplicates by Picard.  Time: 64.63min, Mem: 22G (21,250,637,824)
#module load java
#java -jar ~/local/picard-tools-1.70/MarkDuplicates.jar INPUT=${3}.sorted.bam OUTPUT=${3}.sorted_NoDup.bam COMMENT="Remove Duplicates Using Picard" REMOVE_DUPLICATES=true METRICS_FILE=duplicates.txt ASSUME_SORTED=true
#samtools view ${3}.sorted_NoDup.bam | mapping_stat_BWA.pl ${3}_NoDup > stat_NoDup.txt
