#!/bin/bash
cd ..
/group/ngs/src/qualimap/qualimap bamqc -nt 8 --java-mem-size=24G -nr 5000 -bam "$1" -outdir qualimap -gff "$2" -c -gd HUMAN