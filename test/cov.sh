#!/bin/bash
../targetcov.py \
    system_info.yaml run_info_cov.yaml -t 10 \
    --bam data/sample1.bam --bed data/sample1.bed -o results_cov

