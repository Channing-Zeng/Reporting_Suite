#!/bin/bash
source /etc/profile.d/modules.sh; module load python/64_2.7.3; module load java; module load bedtools; module load samtools;
python /group/ngs/src/ngscat/ngscat.py --bams sample1-ready.bam --bed /Users/vlad/vagrant/variantannotation/test/data/sample1.bed --out NGSCat --reference /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa --saturation y
rm -- "/Users/vlad/vagrant/variantannotation/post_for_bcbio.sh"
