#!/bin/bash

# Usage: runVarDict_SNP.sh sample bam snp
vardict.pl -c 1 -S 2 -E 2 -N $1 -b $2 -p -D $3 > ${1}_snp_vars.txt
