#!/bin/bash
date >&2
module unload python >&2; module load python/2.7.3 >&2
module load java perl bedtools samtools bcbio-nextgen vcftools >&2
echo >&2
echo "$@" >&2
echo >&2
echo >&2
eval $@
echo >&2
date >&2