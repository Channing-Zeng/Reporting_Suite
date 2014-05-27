#!/bin/bash
../varannotate.py \
    system_info.yaml run_info_annotation.yaml -t 10 \
    --var data/sample1.vcf --bam data/sample1.bam -o results_anno