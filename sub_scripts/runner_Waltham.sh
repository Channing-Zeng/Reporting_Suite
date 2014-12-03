#!/bin/bash
date >&2
module unload python >&2; module load python/64_2.7.3 >&2
module unload gcc >&2; module load gcc/4.8.3 >&2
module load java perl bedtools samtools >&2
echo >&2
echo "$@" >&2
echo >&2
echo >&2
eval $@
echo >&2
date >&2