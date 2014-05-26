#!/bin/bash
../varqc.py \
    test/system_info.yaml test/run_info_varqc.yaml -t 10 \
    --var test/data/sample1.vcf -o test/results_qc