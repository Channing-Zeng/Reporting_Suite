#!/usr/bin/env python

from os.path import abspath, dirname, realpath, join
from site import addsitedir
project_dir = abspath(dirname(dirname(realpath(__file__))))
addsitedir(join(project_dir))
addsitedir(join(project_dir, 'ext_modules'))
import sub_scripts.__check_python_version  # checking for python version and adding site dirs inside

import sys


class Region:
    def __init__(self, start, end):
        self.start = start
        self.end = end


class Gene:
    def __init__(self, name, chrom, strand=None):
        self.name = name
        self.chrom = chrom
        self.__chrom_key = self.__make_chrom_key()
        self.strand = strand
        self.regions = []
        self.start = None

    def __make_chrom_key(self):
        CHROMS = [('Y', 23), ('X', 24), ('M', 0)]
        for i in range(22, 0, -1):
            CHROMS.append((str(i), i))

        chrom = self.chrom
        if chrom.startswith('chr'):
            chrom = chrom[3:]
        for (c, i) in CHROMS:
            if chrom == c:
                return i
            if chrom.startswith(c):
                return i + 24

        sys.stderr.write('Cannot parse chromosome ' + self.chrom + '\n')
        return None

    def get_key(self):
        return self.__chrom_key, self.start, self.name

    def sort_regions(self):
        self.regions = sorted(self.regions, key=lambda r: (r.start, r.end))
        self.start = self.regions[0].start

    def merge_regions(self):
        non_overlapping_regions = [self.regions[0]]

        for r in self.regions[1:]:
            if r.start > non_overlapping_regions[-1].end:
                non_overlapping_regions.append(r)
            else:
                non_overlapping_regions[-1].end = r.end

        self.regions = non_overlapping_regions
        return non_overlapping_regions


def main():
    if len(sys.argv) <= 1:
        sys.exit('Usage: ' + __file__ + ' bed_file')

    gene_by_name = dict()

    with open(sys.argv[1]) as inp:
        for l in inp:
            if not l:
                pass
            elif l.startswith('#'):
                sys.stdout.write(l)
            else:
                fields = l[:-1].split('\t')
                if not len(fields) != 4 or len(fields) != 6:
                    sys.exit('Incorrect number of fields: ' + str(len(fields)) +
                             ' (' + ' | '.join(fields) + '). Should be 4 of 6.')
                chrom, start, end, gname = fields[:4]
                # start, end = int(start), int(end)

                gene = gene_by_name.get(gname)
                if gene is None:
                    strand = fields[5] if len(fields) == 6 else None
                    gene = Gene(gname, chrom, strand)
                    gene_by_name[gname] = gene

                gene.regions.append(Region(int(start), int(end)))

    genes = gene_by_name.values()
    for gene in genes:
        gene.sort_regions()

    genes = sorted(genes, key=lambda r: r.get_key())
    for gene in genes:
        for r in gene.merge_regions():
            sys.stdout.write(gene.chrom + '\t' + str(r.start) + '\t' +
                             str(r.end) + '\t' + gene.name)
            if gene.strand:
                sys.stdout.write('\t.\t' + gene.strand)
            sys.stdout.write('\n')


if __name__ == '__main__':
    main()