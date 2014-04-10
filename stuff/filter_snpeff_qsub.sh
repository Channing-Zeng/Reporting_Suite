#!/bin/bash
#$ -V
# #$ -q ngs.q
#$ -q batch.q
#$ -l huge_ram=1
#$ -S /bin/bash
#$ -N filter_annotate
#$ -pe smp 1
#$ -l mem_free=8G
#$ -cwd
#$ -b y
#$ -j y
# Annotation script that takes 1-3 inputs, first being the vcf file name,
# second being an indicator if the vcf is from bcbio's ensemble pipeline ('true' if true) and
# third being 'RNA' if the vcf is from the rna-seq mutect pipeline

set -x

source /etc/profile.d/modules.sh
module load bcbio-nextgen/0.7.6
#export PATH=$PATH:/ngs/RDI/PROGRAMS/snpEff3.5/scripts
#export SnpEff=/ngs/RDI/PROGRAMS/snpEff3.5/

# uncomment in Waltham:
export PATH=$PATH:/group/ngs/src/snpEff/snpEff3.5/scripts
export SnpEff=/group/ngs/src/snpEff/snpEff3.5/
export PERL5LIB=$PERL5LIB:/opt/az/local/bcbio-nextgen/stable/0.7.6/tooldir/lib/perl5/site_perl

sample=$(basename $1 .vcf)
ensemble=$2
RNA=$3

# check if samples need to be split (output of bcbio ensemble pipeline)
if [ "$ensemble" == 'true' ]; then
    vcf-subset -c $(echo $sample | sed s/-ensemble//g) -e $1 > $1".tmp.vcf"
    mv $1 $(echo $1 | sed 's/vcf/combined.vcf/g')
    mv $1".tmp.vcf" $1
fi

# for RNA-seq, filter mutect (only RNA-seq variant caller for now) calls by PASS
if [ "$RNA" == 'RNA' ]; then
    grep -v REJECT $1 > $sample".PASS.vcf"
    sample=$sample".PASS"
    inputvcf=$sample".vcf"
else
    inputvcf=$1
fi

java -jar "$SnpEff"SnpSift.jar annotate -v /ngs/reference_data/genomes/Hsapiens/hg19/variation/dbsnp_137.vcf $inputvcf > $sample".dbsnp.vcf"
if [ -f $sample".vcf" ]; then
    rm $sample".vcf"
fi
sample=$sample".dbsnp"

# annotate COSMIC
java -jar "$SnpEff"SnpSift.jar annotate -v /ngs/reference_data/genomes/Hsapiens/hg19/variation/cosmic-v67_20131024-hg19.vcf $sample".vcf" > $sample".COSMIC.vcf"
rm $sample".vcf"
sample=$sample".COSMIC"

# annotate by dbnsfp
java -jar "$SnpEff"SnpSift.jar dbnsfp -f Uniprot_acc,1000Gp1_AF,\
SIFT_score,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,MutationTaster_score,MutationTaster_pred,\
MutationAssessor_score,MutationAssessor_pred,FATHMM_score,ESP6500_AA_AF,ESP6500_EA_AF,Ensembl_geneid,Ensembl_transcriptid \
-v /ngs/reference_data/genomes/Hsapiens/hg19/dbNSF/dbNSFP2.3/dbNSFP2.3.txt.gz $sample".vcf" > $sample".dbnsf.vcf"
rm $sample".vcf"
sample=$sample".dbnsf"

java -Xmx4g -jar "$SnpEff"snpEff.jar eff -dataDir /ngs/reference_data/genomes/Hsapiens/hg19/snpeff -cancer -noLog -1 -i vcf -o vcf hg19 $sample".vcf" > $sample".snpeff.vcf"
rm $sample".vcf"
sample=$sample".snpeff"

# for RNA-seq, add snpeff and RNA editing sites
if [ "$RNA" == 'RNA' ];then
    vcfannotate -b /ngs/reference_data/genomes/Hsapiens/hg19/variation/Human_AG_all_hg19_INFO.bed -k RNA_editing_site $sample".vcf" > $sample".RNAeditSites.vcf"
    rm $sample".vcf"
    sample=$sample".RNAeditSites"
fi

# add genesets
#java -jar "$SnpEff"SnpSift.jar geneSets \
#-v /ngs/reference_data/genomes/Hsapiens/hg19/variation/msigdb.v4.0.symbols.gmt $sample".vcf" > $sample".msigDB.vcf"
#rm $sample".vcf"
#sample=$sample".msigDB"
# generate tab delimited output file

cat $sample".vcf" | vcfEffOnePerLine.pl \
| java -jar "$SnpEff"SnpSift.jar extractFields - \
CHROM POS ID CNT GMAF REF ALT QUAL FILTER TYPE \
"EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].GENE" \
dbNSFP_SIFT_score dbNSFP_Polyphen2_HVAR_score dbNSFP_Polyphen2_HVAR_pred \
dbNSFP_LRT_score dbNSFP_LRT_pred dbNSFP_MutationTaster_score dbNSFP_MutationTaster_pred \
dbNSFP_MutationAssessor_score dbNSFP_MutationAssessor_pred dbNSFP_FATHMM_score \
dbNSFP_Ensembl_geneid dbNSFP_Ensembl_transcriptid dbNSFP_Uniprot_acc \
dbNSFP_1000Gp1_AC dbNSFP_1000Gp1_AF dbNSFP_ESP6500_AA_AF dbNSFP_ESP6500_EA_AF KGPROD PM PH3 \
AB AC AF DP FS GC HRun HaplotypeScore \
MQ0 QA QD ReadPosRankSum set > $sample".tsv"
#MSigDb \

set +x