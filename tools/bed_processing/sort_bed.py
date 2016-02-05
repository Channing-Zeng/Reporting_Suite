#!/usr/bin/env python
import bcbio_postproc  # do not remove it: checking for python version and adding site dirs inside

from os.path import abspath, dirname, realpath, join
import sys

from source.logger import critical
from source.targetcov.Region import SortableByChrom
from tools.bed_processing.get_chr_lengths import get_chr_lengths


class Region(SortableByChrom):
    def __init__(self, chrom, start, end, other_fields, ref_order=None, genome=None):
        SortableByChrom.__init__(self, chrom, genome)
        self.start = start
        self.end = end
        self.ref_order = ref_order
        self.other_fields = tuple(other_fields)

    def get_key(self):
        return self.ref_order, self.start, self.end, self.other_fields


def main(seq_fpath, genome):
    regions = []

    sys.stderr.write('Genome: ' + genome + '\n')
    chr_lengths = get_chr_lengths(seq_fpath, genome, silence=True)
    chr_order = [chr for (chr, l) in chr_lengths]
    for l in sys.stdin:
        if not l.strip():
            continue
        if l.strip().startswith('#'):
            sys.stdout.write(l)
            continue

        fs = l.strip().split('\t')
        chrom = fs[0]
        start = int(fs[1])
        end = int(fs[2])
        other_fields = fs[3:]
        regions.append(Region(chrom, start, end, other_fields, chr_order.index(chrom), genome))

    for region in sorted(regions, key=lambda r: r.get_key()):
        fs = [region.chrom, str(region.start), str(region.end)]
        fs.extend(region.other_fields)
        sys.stdout.write('\t'.join(fs) + '\n')

    sys.stderr.write('Sorted ' + str(len(regions)) + ' regions.\n')


if __name__ == '__main__':
    if len(sys.argv) <= 2:
        critical('Usage: ' + __file__ + ' path_to_.fa genome_build_name')

    seq_fpath = sys.argv[1]
    genome_build = sys.argv[2]
    main(seq_fpath, genome_build)