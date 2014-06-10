#!/bin/bash
source /etc/profile.d/modules.sh
module load python/64_2.7.3 java bedtools samtools
python /gpfs/group/ngs/src/ngs_reporting/varqc_summary.py "$1" "$2" varqc