#!/bin/bash
source /etc/profile.d/modules.sh; module load python/64_2.7.3; module load java; module load bepython /group/ngs/src/varannotate.py --var ../sample1-mutect.filtered.vcf --bam ../sample1-ready.bam --nt=4 -o annotatione}module load samtools;
python /group/ngs/src/varannotate.py --var ../sample1-mutect.filtered.vcf --bam ../sample1-ready.bam --nt=4 -o annotation
