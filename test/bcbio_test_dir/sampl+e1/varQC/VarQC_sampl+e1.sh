#!/bin/bash
source /etc/profile.d/modules.sh; module load python/64_2.7.3; module load java; module load bedtools; module load samtools;
python /group/ngs/src/varqc.py --var "sampl+e1-mutect.filtered.vcf" -o varQC
rm -- "/Users/vlad/vagrant/variantannotation/test/bcbio_test_dir/sampl+e1/varQC/VarQC_sampl+e1.sh"
