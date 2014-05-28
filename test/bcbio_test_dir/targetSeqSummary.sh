#!/bin/bash
source /etc/profile.d/modules.sh; module load python/64_2.7.3; module load java; module load bedtools; module load samtools;
python /gpfs/group/ngs/src/ngs_reporting/targetcov_summary.py test/bcbio_test_dir test/samples.txt targetSeq
rm -- "/Users/vlad/vagrant/variantannotation/post_for_bcbio.sh"
