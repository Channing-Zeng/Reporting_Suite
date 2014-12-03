#!/bin/bash

# Usage: vartect_hg19.sh bam sample target_bed

/group/cancer_informatics/tools_resources/NGS/bin/checkVar.pl -c 1 -s 2 -e 3 -S 2 -E 3 -g 4 -x 0 -f 0.0075 -N $2 -b $1 $3 > ${2}_vars.txt
/group/cancer_informatics/tools_resources/NGS/bin/teststrandbias.R ${2}_vars.txt > ${2}_vars.txt.t
mv ${2}_vars.txt.t ${2}_vars.txt
/group/cancer_informatics/tools_resources/NGS/bin/var2vcf.pl ${2}_vars.txt > ${2}_vars.vcf
/opt/az/oracle/java/jdk1.7.0_11/bin/java -Xmx4g -jar /group/cancer_informatics/tools_resources/NGS/snpEff/snpEff.jar eff -c /group/cancer_informatics/tools_resources/NGS/snpEff/snpEff.config -d -v -canon hg19 ${2}_vars.vcf > ${2}_vars.eff.vcf
/opt/az/oracle/java/jdk1.7.0_11/bin/java -Xmx4g -jar /group/cancer_informatics/tools_resources/NGS/snpEff/SnpSift.jar annotate -v /ngs/cancer_informatics/GenomeData/human/dbsnp137_00_all.vcf ${2}_vars.eff.vcf > ${2}_vars.eff.dbsnp.vcf
/opt/az/oracle/java/jdk1.7.0_11/bin/java -Xmx4g -jar /group/cancer_informatics/tools_resources/NGS/snpEff/SnpSift.jar annotate -v /ngs/cancer_informatics/GenomeData/human/CosmicCodingMuts_v64_260313_noLimit.vcf ${2}_vars.eff.dbsnp.vcf > ${2}_vars.eff.dbsnp.cosmic.vcf
/group/cancer_informatics/tools_resources/NGS/bin/checkSNV.pl -c 1 -s 2 -e 3 -S 2 -E 3 -g 4 -N $2 -b $1 $3 > ${2}_alignment_summary.txt
touch done.vartect
