#!/bin/bash
../varannotate.py \
    test/system_info.yaml test/run_info_annotation.yaml -t 10 \
    --var test/data/sample1.vcf --bam test/data/sample1.bam -o test/results_anno