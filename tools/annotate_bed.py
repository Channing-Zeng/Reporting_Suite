#!/usr/bin/env python

from os.path import abspath, dirname, realpath, join, exists
from site import addsitedir
project_dir = abspath(dirname(dirname(realpath(__file__))))
addsitedir(join(project_dir))
addsitedir(join(project_dir, 'ext_modules'))
import sub_scripts.__check_python_version  # do not remove it: checking for python version and adding site dirs inside

import subprocess
import sys
import os

from source.utils import OrderedDefaultDict
from collections import defaultdict, OrderedDict


""" Input: Any BED file
    Ouput:
        BED file with regions from input, followed by symbol from best gene overlap from Ensembl.
        Regions can be duplicated, in case if they overlap multiple genes. For each gene, only one record.
        If a region do not overlap any gene, it gets output once in a 3-col line (no symbol is provided).

    Usage: python annotate_bed Input_BED_file work_dir [Reference_BED_file] [bedtools_tool_path] > Annotated_BED_file
"""


def _read_args(args):
    if len(args) < 2:
        log('Usage:\n')
        log('  ' + __file__ + ' Input_BED_file work_dir [Reference_BED_file] [bedtools_tool_path] > Annotated_BED_file\n')
        sys.exit(1)

    input_bed_fpath = abspath(args[0])
    log('Input: ' + input_bed_fpath)

    work_dirpath = abspath(args[1])
    log('Working directory: ' + work_dirpath)
    if not exists(work_dirpath):
        os.mkdir(work_dirpath)

    ref_bed_fpath = '/ngs/reference_data/genomes/Hsapiens/hg19/bed/Exons/Exons.with_genes.bed'
    bedtools = 'bedtools'
    if len(args) > 2:
        if args[2]:
            ref_bed_fpath = abspath(args[2])
            log('Over-set reference fpath: ' + ref_bed_fpath)

    if len(args) > 3:
        bedtools = args[3]
        log('Over-set bedtools: ' + bedtools)
    log()

    return input_bed_fpath, work_dirpath, ref_bed_fpath, bedtools


def log(msg=''):
    sys.stderr.write(msg + '\n')


class Region:
    def __init__(self, chrom, start, end, symbol=None, exon=None, strand=None, feature=None, biotype=None):
        self.chrom = chrom
        self.__chrom_key = self.__make_chrom_key()
        self.start = start
        self.end = end
        self.symbol = symbol
        self.exon = exon
        self.strand = strand
        self.feature = feature
        self.biotype = biotype
        self.total_merged = 0

    def __str__(self):
        fs = [self.chrom, '{}'.format(self.start), '{}'.format(self.end),
              self.symbol or '.', self.exon or '.', self.strand or '.',
              self.feature or '.', self.biotype or '.']

        return '\t'.join(fs) + '\n'

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
        return self.__chrom_key, self.start, self.end, self.symbol


def main():
    input_bed_fpath, work_dirpath, ref_bed_fpath, bedtools = _read_args(sys.argv[1:])

    ref_bed_no_genes_fpath, ref_bed_genes_fpath = _split_reference(work_dirpath, ref_bed_fpath)

    log('Annotating based on CDS and exons...')
    annotated, off_targets = _annotate(bedtools, input_bed_fpath, ref_bed_no_genes_fpath)

    if off_targets:
        off_target_fpath = _save_regions(off_targets, join(work_dirpath, 'off_target_1.bed'))
        log('Saved off target1 to ' + str(off_target_fpath))

        log()
        log('Trying to annotate based on genes rather than CDS and exons...')
        annotated_2, off_targets = _annotate(bedtools, off_target_fpath, ref_bed_genes_fpath)

        for a in annotated_2:
            a.feature = 'UTR/Intron'
        annotated.extend(annotated_2)

        annotated.extend(off_targets)

    log()
    log('Saving annotated regions...')
    for region in sorted(annotated, key=lambda r: r.get_key()):
        sys.stdout.write(str(region))

        # for r, overlap_size in overlaps:
        #     sys.stdout.write('\t' + '\t'.join([
        #         r.chrom, '{:,}'.format(r.start), '{:,}'.format(r.end), r.gene, r.exon, str(r.strand), r.feature, r.biotype,
        #         str(overlap_size),
        #         '{:.2f}%'.format(100.0 * overlap_size / (r.end - r.start))
        #     ]))
        # sys.stdout.write('\n')
    log('Done.')


def merge_fields(consensus_field, other_field):
    if not consensus_field:
        consensus_field = other_field
    else:
        consensus_field = ','.join(set(consensus_field.split(',')) | set(other_field.split(',')))
    return consensus_field


def _resolve_ambiguities(annotated_by_loc_by_gene):

    annotated = []
    for (chrom, start, end), overlaps_by_gene in annotated_by_loc_by_gene.iteritems():
        for g_name, overlaps in overlaps_by_gene.iteritems():
            consensus = Region(chrom, start, end, symbol=g_name, exon='', strand='', feature='', biotype='')
            for r, overlap_size in overlaps:
                if consensus.strand:
                    assert consensus.strand == r.strand, 'Consensus strand is ' + \
                         consensus.strand + ', region strand is ' + r.strand
                else:
                    consensus.strand = r.strand
                consensus.exon = merge_fields(consensus.exon, r.exon)
                consensus.feature = merge_fields(consensus.feature, r.feature)
                consensus.biotype = merge_fields(consensus.biotype, r.biotype)
                consensus.total_merged += 1

            annotated.append(consensus)

    return annotated


def _annotate(bedtools, bed_fpath, ref_fpath):
    cmdline = 'cut -f1,2,3 ' + bed_fpath
    sys.stderr.write(cmdline)
    p = subprocess.Popen(cmdline.split(), stdout=subprocess.PIPE)
    cmdline = '{bedtools} intersect -a - -b {ref_fpath} -wao'.format(**locals())
    sys.stderr.write(' | ' + cmdline + '\n')
    log()
    p = subprocess.Popen(cmdline.split(), stdin=p.stdout, stdout=subprocess.PIPE)

    total_lines = 0
    total_uniq_lines = 0
    total_annotated = 0
    total_uniq_annotated = 0

    met = set()

    annotated_by_loc_by_gene = OrderedDefaultDict(lambda: defaultdict(list))
    off_targets = list()

    for l in p.stdout:
        a_chr, a_start, a_end, e_chr, e_start, e_end, e_gene, e_exon, e_strand, \
            e_feature, e_biotype, overlap_size = l.strip().split('\t')
        assert e_chr == '.' or a_chr == e_chr, str((a_chr + ', ' + e_chr))
        total_lines += 1
        if (a_chr, a_start, a_end) not in met:
            total_uniq_lines += 1

        if e_chr != '.':
            annotated_by_loc_by_gene[(a_chr, int(a_start), int(a_end))][e_gene].append((
                Region(e_chr, int(e_start), int(e_end), e_gene, e_exon, e_strand, e_feature, e_biotype),
                int(overlap_size)))

            total_annotated += 1
            if (a_chr, a_start, a_end) not in met:
                total_uniq_annotated += 1

        else:
            off_targets.append(Region(a_chr, int(a_start), int(a_end)))

        met.add((a_chr, a_start, a_end))

    log('Total intersections, including off-target: ' + str(total_lines))
    log('Total uniq regions in intersections, including off-target: ' + str(total_uniq_lines))
    log('Total annotated regions: ' + str(total_annotated))
    log('Total uniq annotated regions: ' + str(total_uniq_annotated))
    log('Total off target regions: ' + str(len(off_targets)))
    log()

    log('Resolving ambiguities...')
    annotated = _resolve_ambiguities(annotated_by_loc_by_gene)

    return annotated, off_targets


def _save_regions(regions, fpath):
    with open(fpath, 'w') as off_target_f:
        for r in regions:
            off_target_f.write(str(r))

    return fpath


def _split_reference(work_dirpath, ref_bed_fpath):
    log('Splitting reference file into genes and non-genes:')
    ref_bed_no_genes_fpath = join(work_dirpath, os.path.basename(ref_bed_fpath) + '__no_whole_genes')
    ref_bed_genes_fpath = join(work_dirpath, os.path.basename(ref_bed_fpath) + '__whole_genes')
    with open(ref_bed_no_genes_fpath, 'w') as out:
        cmdline = 'grep -wv Gene {ref_bed_fpath} | grep -wv Multi_Gene'.format(**locals())
        log(cmdline + ' > ' + ref_bed_no_genes_fpath)
        subprocess.call(cmdline.split(), stdout=out)
    with open(ref_bed_genes_fpath, 'w') as out:
        cmdline = 'grep -w Gene {ref_bed_fpath}'.format(**locals())
        log(cmdline + ' > ' + ref_bed_genes_fpath)
        subprocess.call(cmdline.split(), stdout=out)
    log()

    return ref_bed_no_genes_fpath, ref_bed_genes_fpath


if __name__ == '__main__':
    main()