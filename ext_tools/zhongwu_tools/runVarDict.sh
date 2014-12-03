#!/bin/bash

# Usage: runVarDict.sh 1.fastq 2.fastq sample target_bed freq
PATH=$PATH:/group/cancer_informatics/tools_resources/NGS/bin
JBIN=/opt/az/oracle/java/jdk1.7.0_11/bin
SBIN=/group/cancer_informatics/tools_resources/NGS/bin
SNPEFF=/group/cancer_informatics/tools_resources/NGS/snpEff
REF=/ngs/reference_data/genomes/Hsapiens/hg19
GATK=/opt/az/local/gatk/GenomeAnalysisTK-2013.2-2.5.2-0-g3ae1219
SAMPLE=$3
BED=$4
FREQ=$5
BEDBASE=`basename $BED`

if [ ! $FREQ ]
    then
        FREQ=0.01
fi

if [ ! -e ${SAMPLE}[._]sorted.realign.bam ]
    then
	/users/kdld047/local/bwa-0.7.4/bwa mem -t 8  -P -c 20 -R "@RG\tID:$SAMPLE\tSM:$SAMPLE" $REF/bwa/hg19.fa $1 $2 > $SAMPLE.sam 2> bwa.err 
	samtools view -Sb -F 0x4 $SAMPLE.sam > $SAMPLE.bam
	samtools sort $SAMPLE.bam ${SAMPLE}.sorted
	samtools index ${SAMPLE}.sorted.bam

	$SBIN/mapping_stat_BWA.pl $SAMPLE $SAMPLE.sam > stat_sam_$SAMPLE.txt
	rm $SAMPLE.sam
	rm $SAMPLE.bam*

	$JBIN/java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -I ${SAMPLE}.sorted.bam -R $REF/seq/hg19.fa -o ${SAMPLE}.intervals -allowPotentiallyMisencodedQuals
	$JBIN/java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -T IndelRealigner -I ${SAMPLE}.sorted.bam -R $REF/seq/hg19.fa -targetIntervals ${SAMPLE}.intervals -allowPotentiallyMisencodedQuals --out ${SAMPLE}.sorted.realign.bam --maxReadsForRealignment 40000
	samtools index ${SAMPLE}.sorted.realign.bam
	rm ${SAMPLE}.sorted.bam*
fi

BAM=`ls ${SAMPLE}[._]sorted.realign.bam`

if [ $BED ]
    then
	$SBIN/splitBed.pl $BED 8
	if [ ! -e done.vardict ]
	    then
	    for n in {1..8}; do
		minivardict.sh $BAM ${SAMPLE} $BEDBASE $n $FREQ &
	    done
	    waitVardict.pl vardict 8
	    cat ${SAMPLE}_vars.txt.[1-9] > ${SAMPLE}_vars.txt
	    $SBIN/teststrandbias.R ${SAMPLE}_vars.txt > ${SAMPLE}_vars.txt.t
	    mv ${SAMPLE}_vars.txt.t ${SAMPLE}_vars.txt
	    $SBIN/var2vcf_valid.pl -S -f $FREQ ${SAMPLE}_vars.txt > ${SAMPLE}_vars.vcf
	    $JBIN/java -Xmx4g -jar $SNPEFF/snpEff.jar eff -c $SNPEFF/snpEff.config -d -v -canon hg19 ${SAMPLE}_vars.vcf > ${SAMPLE}_vars.eff.vcf
	    $JBIN/java -Xmx4g -jar $SNPEFF/SnpSift.jar annotate -v $REF/variation/dbsnp_latest.vcf ${SAMPLE}_vars.eff.vcf > ${SAMPLE}_vars.eff.dbsnp.vcf
	    $JBIN/java -Xmx4g -jar $SNPEFF/SnpSift.jar annotate -v $REF/variation/clinvar_latest.vcf ${SAMPLE}_vars.eff.dbsnp.vcf > ${SAMPLE}_vars.eff.dbsnp.clin.vcf
	    $JBIN/java -Xmx4g -jar $SNPEFF/SnpSift.jar annotate -v $REF/variation/CosmicCodingMuts_latest.vcf ${SAMPLE}_vars.eff.dbsnp.clin.vcf > ${SAMPLE}_vars.eff.dbsnp.clin.cosmic.vcf
	    rm ${SAMPLE}_vars.txt.[1-9]
	    rm vardict.done.*
	    touch done.vardict
	fi

	# For coverage
	if [ ! -e done.checkCov ]
	    then
	    #for n in {1..8}; do
		#minicheckCov.sh $BAM ${SAMPLE} $BEDBASE $n &
	    #done
	    #waitVardict.pl checkCov 8
	    #cat ${SAMPLE}_cov.txt.1 > ${SAMPLE}_cov.txt
	    #cat ${SAMPLE}_cov.txt.[2-9] | grep -v Sample >> ${SAMPLE}_cov.txt
	    #rm ${SAMPLE}_cov.txt.[1-9]
	    #rm checkCov.done.*
	    checkCov.pl -c 1 -s 2 -e 3 -S 2 -E 3 -g 4 -N $SAMPLE -b $BAM -d 1:5:10:25:50:100:500:1000:5000:10000:50000 $BED > ${SAMPLE}_cov.txt
	    touch done.checkCov
	fi

	rm snpEff_summary.html
	rm snpEff_genes.txt

	# For alignment summary
	if [ ! -e done.alignment_sum ]
	    then
	    for n in {1..8}; do
		minicheckSNV.sh $BAM ${SAMPLE} $BEDBASE $n &
	    done
	    waitVardict.pl alignment_sum 8
	    cat ${SAMPLE}_alignment_summary.txt.1 > ${SAMPLE}_alignment_summary.txt
	    cat ${SAMPLE}_alignment_summary.txt.[2-9] | grep -v Sample >> ${SAMPLE}_alignment_summary.txt
	    rm ${SAMPLE}_alignment_summary.txt.[1-9]
	    rm alignment_sum.done.*
	    touch done.alignment_sum
	fi
	rm ${BEDBASE}.[1-8]

	#$SBIN/checkSNV.pl -c 1 -s 2 -e 3 -S 2 -E 3 -g 4 -N $SAMPLE -b ${SAMPLE}.sorted.realign.bam $BED > ${SAMPLE}_alignment_summary.txt
	#$SBIN/checkCov.pl -c 1 -s 2 -e 3 -S 2 -E 3 -g 4 -N $SAMPLE -b ${SAMPLE}.sorted.realign.bam -d 1:5:10:25:50:100:500:1000:5000:10000:50000 $BED > ${SAMPLE}_cov.txt

	#/group/cancer_informatics/tools_resources/NGS/bin/teststrandbias.R ${3}_vars.txt > ${3}_vars.txt.t
	#mv ${3}_vars.txt.t ${3}_vars.txt
	#/group/cancer_informatics/tools_resources/NGS/bin/checkCov.pl -c 1 -s 2 -e 3 -S 2 -E 3 -g 4 -N $3 -b ${3}.sorted.realign.bam -d 1:10:50:100:500:1000:5000:10000:50000 $4 > ${3}_cov.txt
	#/group/cancer_informatics/tools_resources/NGS/bin/var2vcf.pl ${SAMPLE}_vars.txt > ${SAMPLE}_vars.vcf
	#/opt/az/oracle/java/jdk1.7.0_11/bin/java -Xmx4g -jar /group/cancer_informatics/tools_resources/NGS/snpEff/snpEff.jar eff -c /group/cancer_informatics/tools_resources/NGS/snpEff/snpEff.config -d -v -canon hg19 ${SAMPLE}_vars.vcf > ${SAMPLE}_vars.eff.vcf
	#/opt/az/oracle/java/jdk1.7.0_11/bin/java -Xmx4g -jar /group/cancer_informatics/tools_resources/NGS/snpEff/SnpSift.jar annotate -v /ngs/cancer_informatics/GenomeData/human/dbsnp_latest.vcf ${SAMPLE}_vars.eff.vcf > ${SAMPLE}_vars.eff.dbsnp.vcf
	#/opt/az/oracle/java/jdk1.7.0_11/bin/java -Xmx4g -jar /group/cancer_informatics/tools_resources/NGS/snpEff/SnpSift.jar annotate -v /ngs/cancer_informatics/GenomeData/human/CosmicCodingMuts_latest.vcf ${SAMPLE}_vars.eff.dbsnp.vcf > ${SAMPLE}_vars.eff.dbsnp.cosmic.vcf
fi

touch runBWA074_hg19.done

#samtools view $3.bam | mapping_stat_BWA.pl ${3}_bam > stat_bam.txt
#samtools view ${3}_sorted.bam | mapping_stat_BWA.pl ${3}_sorted > stat_sorted_bam.txt

# Remove duplicates by Picard.  Time: 64.63min, Mem: 22G (21,250,637,824)
#module load java
#java -jar ~/local/picard-tools-1.70/MarkDuplicates.jar INPUT=${3}_sorted.bam OUTPUT=${3}_sorted_NoDup.bam COMMENT="Remove Duplicates Using Picard" REMOVE_DUPLICATES=true METRICS_FILE=duplicates.txt ASSUME_SORTED=true
#samtools view ${3}_sorted_NoDup.bam | mapping_stat_BWA.pl ${3}_NoDup > stat_NoDup.txt
