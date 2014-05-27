#!/bin/bash
source /etc/profile.d/modules.sh; module load python/64_2.7.3; module load java; module load bedtools; module load samtools;
python /group/ngs/bin/InDelFilter.py sample1-mutect.vcf > sample1-mutect.filtered.vcf
