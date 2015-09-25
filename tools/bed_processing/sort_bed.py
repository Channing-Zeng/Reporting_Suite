#!/usr/bin/env python

import __check_python_version  # do not remove it: checking for python version and adding site dirs inside

from os.path import abspath, dirname, realpath, join
import sys


HG19_CHROMS = [('X', 23), ('Y', 24), ('M', 0), ('Un', 25)]
for i in range(22, 0, -1):
    HG19_CHROMS.append((str(i), i))

MM10_CHROMS = [('X', 22), ('Y', 23), ('M', 24), ('Un', 25)]
for i in range(21, -1, -1):
    MM10_CHROMS.append((str(i), i))
#
# HG19_CHROMS_plus = [('Un_', 49)]
# for c, i in HG19_CHROMS:
#     HG19_CHROMS_plus.append((str(i) + '_', i + 25))
#     HG19_CHROMS_plus.append((str(i), i))
# HG19_CHROMS = HG19_CHROMS_plus
#
# MM10_CHROMS_plus = [('Un_', 49)]
# for c, i in MM10_CHROMS:
#     MM10_CHROMS_plus.append((str(i) + '_', i + 25))
#     MM10_CHROMS_plus.append((str(i), i))
# MM10_CHROMS = MM10_CHROMS_plus
#
# print str(HG19_CHROMS)


class SortableByChrom:
    def __init__(self, chrom, genome):
        self.chrom = chrom
        self.genome = genome
        self._chrom_key = self.__make_chrom_key(genome)

    def __make_chrom_key(self, genome):
        chroms = HG19_CHROMS if self.genome != 'mm10' else MM10_CHROMS

        chr_remainder = self.chrom
        if self.chrom.startswith('chr'):
            chr_remainder = self.chrom[3:]
        for (c, i) in chroms:
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
        SortableByChrom.__init__(self, chrom, genome)
        self.start = start
        self.end = end
        self.other_fields = tuple(other_fields)

    def get_key(self):
        return self._chrom_key, self.start, self.end, self.other_fields


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