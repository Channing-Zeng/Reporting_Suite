#!/bin/bash
#set -x
date >&2
hostname >&2
#source /etc/profile.d/modules.sh >&2
module unload python >&2 2>&2
module unload gcc >&2 2>&2
module load gcc/4.9.2 python/64_2.7.3 java perl bedtools samtools >&2 2>&2
echo >&2
echo "$@" >&2
echo >&2
echo >&2
eval $@
echo >&2
date >&2
#set +x