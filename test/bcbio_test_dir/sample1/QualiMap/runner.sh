#!/bin/bash
source /etc/profile.d/modules.sh; module load python/64_2.7.3; module load java; module load bedtools; module load samtools;
/group/ngs/src/qualimap/qualimap bamqc -nt 8 --java-mem-size=24G -nr 5000 -bam ../sample1-ready.bam -outdir QualiMap -gff test/data/sample1.bed -c -gd HUMAN
