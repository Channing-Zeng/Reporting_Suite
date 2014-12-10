#!/bin/bash
date >&2
hostname >&2
source /etc/profile.d/modules.sh >&2
module load java perl bcbio-nextgen >&2
echo >&2
echo "$@" >&2
echo >&2
echo >&2
eval $@
echo >&2
date >&2