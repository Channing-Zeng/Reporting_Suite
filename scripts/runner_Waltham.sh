#!/bin/bash
#set -x
date >&2
hostname >&2
source /etc/profile.d/modules.sh >&2
module load gcc/4.9.2 sge bedtools/2.24.0 bcbio bedops >&2 2>&2
echo >&2
echo "$2" >&2
echo >&2
echo >&2
eval $2
echo >&2
date >&2
touch $1 >&2
#set +x
