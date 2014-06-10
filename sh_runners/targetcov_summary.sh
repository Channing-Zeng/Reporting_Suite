#!/bin/bash
source /etc/profile.d/modules.sh
module load python/64_2.7.3 java bedtools samtools
python /gpfs/group/ngs/src/ngs_reporting/targetcov_summary.py "$0" "$1" targetcov