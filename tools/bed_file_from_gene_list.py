#!/usr/bin/env python

import sys
if len(sys.argv) < 2:
    sys.stderr.write('Usage: ' + __file__ + ' genes.txt [regions.bed] > regions_for_genes.txt')
    sys.stderr.write('    the script filters regions to have gene names in genes.txt')
    sys.exit(1)

genes = set(open(sys.argv[1]).read().split())

exons_fpath = '/ngs/reference_data/genomes/Hsapiens/hg19/bed/Exons/Exons.with_genes.bed'
if len(sys.argv) > 2:
    exons_fpath = sys.argv[2]

exons = [l.split('\t') for l in open(exons_fpath).read().split('\n')]
for ts in exons:
    if len(ts) < 4:
        pass
    else:
        g = ts[3]
        if g in genes:
            print '\t'.join(ts)
