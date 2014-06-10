#!/bin/bash
source /etc/profile.d/modules.sh
module load python/64_2.7.3 java bedtools samtools
python /group/ngs/src/ngscat/ngscat.py --bams "$1" --bed "$2" --out ngscat --saturation y
