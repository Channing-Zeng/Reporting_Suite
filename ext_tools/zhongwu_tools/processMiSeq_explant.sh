#!/bin/bash

# Usage: processMiSeq.sh fastq_dir bed

ls $1/*fastq.gz | pairFastq.pl -n "([^\/]*)_[SL]\d+" > samples.txt

for i in `getcols -c 2 samples.txt`; do countReads.pl -n "([^\/]*)_[SL]\d+" $i; done > read_counts.txt &

sample2align.pl -b $2 -p /group/cancer_informatics/tools_resources/NGS/bin/runBWA074_explant.sh samples.txt | addQueue2.pl > runAlign_explant.sh
bash runAlign_explant.sh 

mkdir PhiX
cd PhiX
sample2align.pl -p /group/cancer_informatics/tools_resources/NGS/bin/runBWA074_phiX.sh ../samples.txt | addQueue.pl > runAlign_PhiX.sh
bash runAlign_PhiX.sh
cd ..

#vcf2txt.pl -p 5 -q 25 */*cos* | pickLine -c 35 -l PASS:TRUE | pickLine -c 13 -l ACAA1 -v | getcols -c 1..6:13:9:11:10:16:17..19:21..24:26:29:30:31..36 > variants_filtered.txt
