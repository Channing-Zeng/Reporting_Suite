#!/bin/bash
../varqc.py \
    system_info.yaml run_info_varqc.yaml -t 10 \
    --var data/sample1.vcf -o results_qc