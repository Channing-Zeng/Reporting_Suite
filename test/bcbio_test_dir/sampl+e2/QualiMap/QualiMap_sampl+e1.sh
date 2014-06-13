#!/bin/bash
source /etc/profile.d/modules.sh; module load python/64_2.7.3; module load java; module load bedtools; module load samtools;
/group/ngs/src/qualimap/qualimap bamqc -nt 8 --java-mem-size=24G -nr 5000 -bam "sampl+e1-ready.bam" -outdir QualiMap -gff "/ngs/tmp/bed_fixed_for_qualimap.bed" -c -gd HUMAN
rm -- "/Users/vlad/vagrant/variantannotation/test/bcbio_test_dir/sampl+e1/QualiMap/QualiMap_sampl+e1.sh"
