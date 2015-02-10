#!/usr/bin/env python

from os.path import abspath, dirname, realpath, join
from site import addsitedir
project_dir = abspath(dirname(dirname(realpath(__file__))))
addsitedir(join(project_dir))
addsitedir(join(project_dir, 'ext_modules'))
import sub_scripts.__check_python_version  # do not remove it: checking for python version and adding site dirs inside

import sys


class Exon:
    def __init__(self, start, end, biotype=None, feature=None):
        self.start = start
        self.end = end
        self.biotype = biotype
        self.feature = feature


class Gene:
    def __init__(self, name, chrom, strand=None):
        self.name = name
        self.chrom = chrom
        self.__chrom_key = self.__make_chrom_key()
        self.strand = strand
        self.biotype = None
        self.start = None
        self.end = None
        self.feature = 'gene'
        self.met_gene_feature = False  # some bed files can contain 'gene' features, so we can take start and end from them

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
        if len(self.regions) == 0:
            sys.stderr.write('Error: sub-regions of ' + str(self) + ' is 0' + '\n')
            sys.exit(1)

        non_overlapping_regions = [self.regions[0]]

        for r in self.regions[1:]:
            if r.start > non_overlapping_regions[-1].end:
                non_overlapping_regions.append(r)
            else:
                prev_r = non_overlapping_regions[-1]
                prev_r.end = r.end
                prev_r.biotype = merge_fields(prev_r.biotype, r.biotype)
                prev_r.feature = merge_fields(prev_r.feature, r.feature)

        self.regions = non_overlapping_regions
        return non_overlapping_regions

    def __repr__(self):
        return self.chrom + ':' + str(self.start) + '-' + str(self.end) + ',' + str(self.name)

    def __str__(self):
        return self.__repr__()


def merge_fields(consensus_field, other_field):
    if not consensus_field:
        consensus_field = other_field
    else:
        consensus_field = ','.join(set(consensus_field.split(',')) | set(other_field.split(',')))
    return consensus_field


def main():
    if len(sys.argv) <= 1:
        sys.exit('Usage: ' + __file__ + ' bed_file')

    summarize_by_genes = True

    gene_by_chrom_and_name = dict()

    i = 0
    total_genes_inp = 0
    with open(sys.argv[1]) as inp:
        for l in inp:
            if not l:
                pass
            elif l.startswith('#'):
                sys.stdout.write(l)
            else:
                fields = l[:-1].split('\t')

                if len(fields) < 3:
                    sys.exit('Incorrect number of fields: ' + str(len(fields)) +
                             ' (' + ' | '.join(fields) + '). Should be >= 3.')
                else:
                    if len(fields) <= 3 and summarize_by_genes is True:
                        summarize_by_genes = False
                        sys.stderr.write('3 columns in BED; no summarizing by genes\n')

                    chrom, start, end = fields[:3]
                    start, end = int(start), int(end)
                    gname = fields[3] if len(fields) >= 4 else '.'
                    strand = fields[5] if len(fields) >= 6 else None
                    (feature, biotype) = fields[6:8] if len(fields) >= 8 else (None, None)

                    gene = gene_by_chrom_and_name.get((chrom, gname))
                    if gene is None:
                        gene = Gene(gname, chrom, strand)
                        gene_by_chrom_and_name[(chrom, gname)] = gene

                    if feature == 'gene':  # in fact 'gene' features in BED files are optional
                        total_genes_inp += 1

                        if gene.met_gene_feature:
                            # miltiple records for gene, picking the lowest start and the biggest end
                            gene.start = min(gene.start, start)
                            gene.end = max(gene.end, end)
                            gene.biotype = merge_fields(gene.biotype, biotype)
                            assert gene.strand == strand, 'Prev gene strand is ' + gene.strand + ', new strand is ' + strand

                        else:
                            gene.start = start
                            gene.end = end
                            gene.biotype = biotype
                            gene.met_gene_feature = True

                    elif feature is None or feature == 'CDS' or feature == 'exon':
                        gene.regions.append(Exon(int(start), int(end), biotype, feature))
            i += 1
            if i % 1000 == 0:
                sys.stderr.write('processed ' + str(i) + ' lines\n')
                sys.stderr.flush()
    sys.stderr.write('Processed ' + str(i) + ' lines, found ' + str(total_genes_inp) + ' genes, ' + str(len(gene_by_chrom_and_name)) + ' uniq gene names.\n')
    sys.stderr.write('\n')

    sys.stderr.write('Sorting regions...\n')
    genes = []
    for gene in gene_by_chrom_and_name.values():
        if gene.sort_regions() is not None:
            genes.append(gene)

    sys.stderr.write('Merging regions...\n')
    final_regions = []
    for gene in sorted(genes, key=lambda g: g.get_key()):
        if summarize_by_genes:
            final_regions.append((gene.chrom, gene.start, gene.end, gene.name, gene.strand, gene.feature, gene.biotype))

        for r in gene.merge_regions():
            final_regions.append((gene.chrom, r.start, r.end, gene.name, gene.strand, r.feature, r.biotype))

    sys.stderr.write('Merged, regions after merge: ' + str(len(final_regions)) + ', saving...\n')

    for chrom, start, end, gname, strand, feature, biotype in sorted(final_regions):
        sys.stdout.write('\t'.join([chrom, str(start), str(end)]))
        if summarize_by_genes:
            sys.stdout.write('\t' + '\t'.join([gname, '.', strand or '.', feature or '.', biotype or '.']))
        sys.stdout.write('\n')
    sys.stderr.write('Saved\n')


if __name__ == '__main__':
    main()