#!/bin/bash

# Usage: processMiSeq.sh fastq_dir bed freq

ls $1/*R[12]*fastq.gz | pairFastq.pl -n "([^\/]*)_[RS]\d+[_\.].*" > samples.txt # for MiSeq file names
#ls $1/*R[12]*fastq.gz | pairFastq.pl -n "([^\/]*)_R[12].*" > samples.txt # for HiSeq file names

#for i in `getcols -c 2 samples.txt`; do countReads.pl -n "([^\/]*)_S\d+" $i; done > read_counts.txt & # for MiSeq processing
for i in `getcols -c 2 samples.txt`; do countReads.pl -n "([^\/]*)_[RS]\d+[_\.].*" $i; done > read_counts.txt &

#sample2align.pl -b $2 samples.txt | addQueue2.pl > runAlign_hg19.sh
sample2align.pl -p /group/cancer_informatics/tools_resources/NGS/bin/runVarDict.sh -b $2 -f $3 samples.txt > runAlign_hg19.sh
bash runAlign_hg19.sh 

mkdir PhiX
cd PhiX
#sample2align.pl -p /group/cancer_informatics/tools_resources/NGS/bin/runBWA074_phiX.sh ../samples.txt | addQueue.pl > runAlign_PhiX.sh
sample2align.pl -p /group/cancer_informatics/tools_resources/NGS/bin/runBWA074_phiX.sh ../samples.txt > runAlign_PhiX.sh
bash runAlign_PhiX.sh
cd ..

#vcf2txt.pl -p 5 -q 25 */*cos* | pickLine -c 35 -l PASS:TRUE | pickLine -c 13 -l ACAA1 -v | getcols -c 1..6:13:9:11:10:16:17..19:21..24:26:29:30:31..36 > variants_filtered.txt
