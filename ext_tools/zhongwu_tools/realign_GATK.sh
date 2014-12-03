#!/bin/bash

# The script will perform local realignment using GATK
# Usage: realign_GATK.sh bam_file

base=`basename $1 .bam`
/opt/az/sun/java/jdk64_1.6.0_24/bin/java -Xmx8g -jar /opt/az/local/gatk/GenomeAnalysisTK-2013.2-2.5.2-0-g3ae1219/GenomeAnalysisTK.jar -T RealignerTargetCreator -I $1 -R /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Sequence/WholeGenomeFasta/genome.fa -o ${base}.intervals -allowPotentiallyMisencodedQuals

/opt/az/sun/java/jdk64_1.6.0_24/bin/java -Xmx8g -jar /opt/az/local/gatk/GenomeAnalysisTK-2013.2-2.5.2-0-g3ae1219/GenomeAnalysisTK.jar -T IndelRealigner -I $1 -R /group/cancer_informatics/tools_resources/NGS/genomes/hg19/Sequence/WholeGenomeFasta/genome.fa -targetIntervals ${base}.intervals -allowPotentiallyMisencodedQuals --out ${base}.realign.bam
