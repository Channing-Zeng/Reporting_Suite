#!/bin/bash

# Usage: bcl2fastq.sh inputdir samplesheet mask

DIR=$1
SHEET=$2
MASK=$3
BCLDIR=/group/cancer_informatics/tools_resources/NGS/bcl2fastq-1.8.4/bcl2fastq

if [ ! $MASK ]
    then
        MASK=Y100n,I6n,Y100n
fi

/scratch/cancer_informatics/tools_resources/perl/perl-5.8.9/bin/perl $BCLDIR/bin/configureBclToFastq.pl --input-dir $DIR/Data/Intensities/BaseCalls --output-dir $DIR/Unalign --sample-sheet $SHEET --no-eamss --mismatches 1 --adapter-sequence $BCLDIR/share/bcl2fastq-1.8.4/adapters/TruSeq_r1.fa --adapter-sequence $BCLDIR/share/bcl2fastq-1.8.4/adapters/TruSeq_r2.fa --use-bases-mask $MASK --fastq-cluster-count 0 --force

cd Unalign/
nohup make -j 10

