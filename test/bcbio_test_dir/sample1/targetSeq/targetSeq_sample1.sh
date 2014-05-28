#!/bin/bash
source /etc/profile.d/modules.sh; module load python/64_2.7.3; module load java; module load bedtools; module load samtools;
python /group/ngs/src/targetcov.py --bam sample1-ready.bam --bed /Users/vlad/vagrant/variantannotation/test/data/sample1.bed --nt=4 -o targetSeq
rm -- "/Users/vlad/vagrant/variantannotation/post_for_bcbio.sh"
