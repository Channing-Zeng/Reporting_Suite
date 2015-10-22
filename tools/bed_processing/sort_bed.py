#!/usr/bin/env python
import bcbio_postproc  # do not remove it: checking for python version and adding site dirs inside

from os.path import abspath, dirname, realpath, join
import sys
from source.targetcov.Region import SortableByChrom


class Region(SortableByChrom):
    def __init__(self, chrom, start, end, other_fields, genome=None):
        SortableByChrom.__init__(self, chrom, genome)
        self.start = start
        self.end = end
        self.other_fields = tuple(other_fields)

    def get_key(self):
        return SortableByChrom.get_key(self), self.start, self.end, self.other_fields


def main(args):
    regions = []

    genome = None
    if len(args) > 0:
        genome = args[0]
    # if len(sys.argv) > 2:
    #     genome = sys.argv[2]
        sys.stderr.write('Genome: ' + genome + '\n')

    for l in sys.stdin:
        if not l.strip():
            continue
        if l.strip().startswith('#'):
            sys.stdout.write(l)

        fs = l[:-1].split('\t')
        chrom = fs[0]
        start = int(fs[1])
        end = int(fs[2])
        other_fields = fs[3:]
        regions.append(Region(chrom, start, end, other_fields, genome))

    sys.stderr.write('Found ' + str(len(regions)) + ' regions.\n')

    for region in sorted(regions, key=lambda r: r.get_key()):
        fs = [region.chrom, str(region.start), str(region.end)]
        fs.extend(region.other_fields)
        sys.stdout.write('\t'.join(fs) + '\n')


if __name__ == '__main__':
    main(sys.argv[1:])