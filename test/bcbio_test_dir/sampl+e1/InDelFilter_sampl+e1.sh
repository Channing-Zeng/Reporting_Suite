#!/bin/bash
source /etc/profile.d/modules.sh; module load python/64_2.7.3; module load java; module load bedtools; module load samtools;
python /group/ngs/bin/InDelFilter.py "sampl+e1-mutect.vcf" > "sampl+e1-mutect.filtered.vcf"
rm -- "/Users/vlad/vagrant/variantannotation/test/bcbio_test_dir/sampl+e1/InDelFilter_sampl+e1.sh"
