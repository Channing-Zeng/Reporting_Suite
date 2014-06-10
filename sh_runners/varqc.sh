#!/bin/bash
source /etc/profile.d/modules.sh
module load python/64_2.7.3 java bedtools samtools
python /group/ngs/src/varqc.py --var "$0" -o varqc