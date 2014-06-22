#!/bin/bash

/users/kdld047/work/NGS/NGS/tophat-2.0.3.Linux_x86_64/tophat -r 20 --mate-std-dev 35 -p 12 -G ~/work/NGS/NGS/genomes/mm10/Annotation/Genes/genes.gtf -o tophat2_out_mouse --rg-id $1 --rg-sample $1 ~/work/NGS/NGS/genomes/mm10/Sequence/BowtieIndex/bowtie2/genome $2 $3 2> mouse_tophat2.err

