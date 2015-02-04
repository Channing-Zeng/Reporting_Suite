#!/usr/bin/env python

from os.path import abspath, dirname, realpath, join
from site import addsitedir
project_dir = abspath(dirname(dirname(realpath(__file__))))
addsitedir(join(project_dir))
addsitedir(join(project_dir, 'ext_modules'))
import sub_scripts.__check_python_version  # do not remove it: checking for python version and adding site dirs inside

import sys


class Exon:
    def __init__(self, start, end, biotype=None):
        self.start = start
        self.end = end
        self.biotype = biotype


class Gene:
    def __init__(self, name, chrom, strand=None):
        self.name = name
        self.chrom = chrom
        self.__chrom_key = self.__make_chrom_key()
        self.strand = strand
        self.biotype = None
        self.start = None
        self.end = None

        self.regions = []

    def __make_chrom_key(self):
        CHROMS = [('Y', 23), ('X', 24), ('M', 0)]
        for i in range(22, 0, -1):
            CHROMS.append((str(i), i))

        chr_remainder = self.chrom
        if self.chrom.startswith('chr'):
            chr_remainder = self.chrom[3:]
        for (c, i) in CHROMS:
            if chr_remainder == c:
                return i
            elif chr_remainder.startswith(c):
                return i + 24

        sys.stderr.write('Cannot parse chromosome ' + self.chrom + '\n')
        return None

    def get_key(self):
        return self.__chrom_key, self.start, self.end, self.name

    def sort_regions(self):
        self.regions = sorted(self.regions, key=lambda r: (r.start, r.end))
        if self.start is None or self.end is None:
            if not self.regions:
                return None  # no coordinates and no exons to infer coordinates
            else:
                self.start = self.regions[0].start
                self.end = self.regions[-1].end
        return self.regions

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
        sys.exit('Usage: ' + __file__ + ' bed_file [--exons]')

    running_with_exons = len(sys.argv) >= 3
    
    gene_by_chrom_and_name = dict()

    with open(sys.argv[1]) as inp:
        for l in inp:
            if not l:
                pass
            elif l.startswith('#'):
                sys.stdout.write(l)
            else:
                fields = l[:-1].split('\t')

                if len(fields) < 4:
                    sys.exit('Incorrect number of fields: ' + str(len(fields)) +
                             ' (' + ' | '.join(fields) + '). Should be >= 4.')

                else:
                    chrom, start, end, gname = fields[:4]
                    start, end = int(start), int(end)
                    strand = fields[5] if len(fields) >= 6 else None
                    (feature, biotype) = fields[6:8] if len(fields) >= 8 else (None, None)

                    gene = gene_by_chrom_and_name.get((chrom, gname))
                    if gene is None:
                        gene = Gene(gname, chrom, strand)
                        gene_by_chrom_and_name[(chrom, gname)] = gene

                    if feature == 'gene':
                        gene.biotype = biotype
                        gene.start = start
                        gene.end = end

                    elif feature is None or feature == 'exon':
                        gene.regions.append(Exon(int(start), int(end), biotype))

    genes = []
    for gene in gene_by_chrom_and_name.values():
        if gene.sort_regions() is not None:
            genes.append(gene)

    final_regions = []
    for gene in sorted(genes, key=lambda g: g.get_key()):
        final_regions.append((gene.chrom, gene.start, gene.end, gene.name, gene.strand, 'gene', gene.biotype))
        for r in gene.merge_regions():
            final_regions.append((gene.chrom, r.start, r.end, gene.name, gene.strand, 'exon', r.biotype))

    for chrom, start, end, gname, strand, feature, biotype in sorted(final_regions):
        sys.stdout.write('\t'.join([chrom, str(start), str(end), gname, '.', strand or '.', feature, biotype or '.']) + '\n')


if __name__ == '__main__':
    main()