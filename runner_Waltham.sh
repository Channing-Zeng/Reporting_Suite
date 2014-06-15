#!/bin/bash
date >&2
source /etc/profile.d/modules.sh >&2
module load python/64_2.7.3 java bedtools samtools >&2
eval $1 && date >&2