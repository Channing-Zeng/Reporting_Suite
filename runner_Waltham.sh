#!/bin/bash
date >&2
source /etc/profile.d/modules.sh >&2
module unload python
module load python/64_2.7.3 java perl bedtools samtools >&2
module unload gcc
module load gcc/4.8.3
echo >&2
echo "$@" >&2
echo >&2
echo >&2
eval $@
echo >&2
date >&2