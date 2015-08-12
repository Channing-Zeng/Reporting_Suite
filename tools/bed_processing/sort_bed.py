#!/usr/bin/env python

import __check_python_version  # do not remove it: checking for python version and adding site dirs inside

from os.path import abspath, dirname, realpath, join
import sys


class SortableByChrom:
    def __init__(self, chrom, genome=None):
        self.chrom = chrom
        self.genome = genome
        self._chrom_key = self.__make_chrom_key()

    def __make_chrom_key(self):
        CHROMS = [('Y', 23), ('X', 24), ('M', 0)]
        if self.genome == 'mm10':
            CHROMS = [('Y', 23), ('X', 24), ('M', 25)]
        for i in range(22, 0, -1):
            CHROMS.append((str(i), i))

        chr_remainder = self.chrom
        if self.chrom.startswith('chr'):
            chr_remainder = self.chrom[3:]
        for (c, i) in CHROMS:
            if chr_remainder == c:
                return i
            elif chr_remainder.startswith(c):
                return i + 25

        sys.stderr.write('Cannot parse chromosome ' + self.chrom + '\n')
        return None

    def get_key(self):
        return self._chrom_key


class Region(SortableByChrom):
    def __init__(self, chrom, start, end, other_fields, genome=None):
        SortableByChrom.__init__(self, chrom)
        self.chrom = chrom
        self.start = start
        self.end = end
        self.other_fields = tuple(other_fields)
        self.genome = genome

    def get_key(self):
        return self._chrom_key, self.start, self.end, self.other_fields


def main():
    regions = []

    genome = None
    if len(sys.argv) > 1:
        genome = sys.argv[1]

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
    main()