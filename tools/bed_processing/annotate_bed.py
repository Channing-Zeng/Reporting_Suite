#!/usr/bin/env python
import __check_python_version  # do not remove it: checking for python version and adding site dirs inside

import subprocess
import sys
import os
import tempfile
import shutil
from os.path import abspath, dirname, realpath, join, exists, basename, splitext
from source.calling_process import call

from source.logger import critical
from source.main import read_opts_and_cnfs
from source.prepare_args_and_cnf import check_system_resources
from source.prepare_args_and_cnf import check_genome_resources
from source.targetcov.Region import SortableByChrom
from source.targetcov.bam_and_bed_utils import intersect_bed
from source.utils import OrderedDefaultDict
from collections import defaultdict, OrderedDict
from os.path import getsize
from source.file_utils import verify_file, adjust_path
from pybedtools import BedTool


usage = """
    Input: Any BED file
    Output:
        BED file with regions from input, followed by symbol from best gene overlap from Ensembl.
        Regions can be duplicated, in case if they overlap multiple genes. For each gene, only one record.
        If a region do not overlap any gene, it gets output once in a 3-col line (no symbol is provided).

    Usage: python %s Input_BED_file [Reference_BED_file] [bedtools_tool_path] > Annotated_BED_file
""" % __file__


def _read_args(args):
    if len(args) < 1:
        sys.exit(['Usage:',
        '  ' + __file__ + ' Input_BED_file [Reference_BED_file] [bedtools_tool_path] > Annotated_BED_file'])

    input_bed_fpath = abspath(args[0])
    log('Input: ' + input_bed_fpath)

    # work_dirpath = abspath(args[1])
    work_dirpath = tempfile.mkdtemp()
    log('Creating a temporary working directory ' + work_dirpath)
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


def main():
    if len(sys.argv[1]) < 0:
        critical('Usage: ' + __file__ + ' Input_BED_file -o Annotated_BED_file')
    input_bed_fpath = verify_file(sys.argv[1], is_critical=True, description='Input BED file for ' + __file__)

    cnf = read_opts_and_cnfs(
        description='Annotating BED file based on Ensemble exons and genes annotations.',
        extra_opts=[
            (['--reference'], dict(
                dest='reference')
            ),
        ],
        required_keys=['output_file'],
        file_keys=['reference'],
        key_for_sample_name=None,
        fpath_for_sample_name=input_bed_fpath,
        output_fpath_not_dir=True
    )
    check_system_resources(cnf)
    check_genome_resources(cnf)

    exons_fpath = adjust_path(cnf.exons) if cnf.exons else adjust_path(cnf.genome.exons)
    if not verify_file(exons_fpath, 'Ensemble file with exons w/genes for annotate_bed'):
        critical('Ensemble file with exons w/genes for annotate_bed is required')

    genes_exons_bed, no_genes_exons_bed = _split_reference(cnf, exons_fpath)
    bed = BedTool(input_bed_fpath).cut([0, 1, 2])

    log('Annotating based on CDS and exons...')

    annotated, off_targets = _annotate(cnf, bed, no_genes_exons_bed)

    if off_targets:
        off_target_bed = BedTool([(r.chrom, r.start, r.end) for r in off_targets])
        # off_target_fpath = _save_regions(off_targets, join(work_dirpath, 'off_target_1.bed'))
        # log('Saved off target1 to ' + str(off_target_fpath))
        log()
        log('Trying to annotate based on genes rather than CDS and exons...')
        annotated_2, off_targets = _annotate(cnf, off_target_bed, genes_exons_bed)

        for a in annotated_2:
            a.feature = 'UTR/Intron/Decay'
        annotated.extend(annotated_2)

        annotated.extend(off_targets)

    log()
    log('Saving annotated regions to ' + str(cnf.output_file))
    with open(cnf.output_file, 'w') as out:
        for region in sorted(annotated, key=lambda r: r.get_key()):
            out.write(str(region))

        # for r, overlap_size in overlaps:
        #     sys.stdout.write('\t' + '\t'.join([
        #         r.chrom, '{:,}'.format(r.start), '{:,}'.format(r.end), r.gene, r.exon, str(r.strand), r.feature, r.biotype,
        #         str(overlap_size),
        #         '{:.2f}%'.format(100.0 * overlap_size / (r.end - r.start))
        #     ]))
        # sys.stdout.write('\n')
    log('Done.')


def log(msg=''):
    sys.stderr.write(msg + '\n')


class Region(SortableByChrom):
    def __init__(self, chrom, start, end, genome_build=None, gene_symbol=None, exon=None, strand=None, feature=None, biotype=None):
        SortableByChrom.__init__(self, chrom, genome_build)
        self.chrom = chrom
        self.start = start
        self.end = end
        self.symbol = gene_symbol
        self.exon = exon
        self.strand = strand
        self.feature = feature
        self.biotype = biotype
        self.total_merged = 0

    def __str__(self):
        fs = [self.chrom, '{}'.format(self.start), '{}'.format(self.end), self.symbol or '.']
        return '\t'.join(fs) + '\n'

    def get_key(self):
        return SortableByChrom.get_key(self), self.start, self.end, self.symbol


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
            consensus = Region(chrom, start, end, gene_symbol=g_name, exon='', strand='', feature='', biotype='')
            for r, overlap_size in overlaps:
                if consensus.strand:
                    # RefSeq has exons from different strands with the same gene name (e.g. CTAGE4 for hg19),
                    # Such pair of exons may overlap with a single region, so taking strand from the first one
                    if consensus.strand != r.strand:
                        log('Warning: different strands between consensus and next region (gene: ' + g_name + ')')
                    #assert consensus.strand == r.strand, 'Consensus strand is ' + \
                    #     consensus.strand + ', region strand is ' + r.strand
                else:
                    consensus.strand = r.strand
                consensus.exon = merge_fields(consensus.exon, r.exon)
                consensus.feature = merge_fields(consensus.feature, r.feature)
                consensus.biotype = merge_fields(consensus.biotype, r.biotype)
                consensus.total_merged += 1

            annotated.append(consensus)

    return annotated


def _annotate(cnf, bed, ref_bed):
        # off_targets = []
        # with open(bed_fpath) as f:
        #     for l in f:
        #         if l.startswith('#'):
        #             continue
        #         a_chr, a_start, a_end = l.strip().split('\t')[:3]
        #         off_targets.append(Region(a_chr, int(a_start), int(a_end)))
        # return [], off_targets

    # intersect_bed(cnf, '<(cut -f1,2,3 ' + bed_fpath + ')', ref_fpath, output_fpath=)
    # cmdline = '{bedtools} intersect -a - -b {ref_fpath} -wao'.format(**locals())
    # call(cnf, cmdline, output_fpath=)
    # p = subprocess.Popen(cmdline.split(), stdin=p.stdout, stdout=subprocess.PIPE)
    intersection = bed.intersect(ref_bed, wao=True)

    total_lines = 0
    total_uniq_lines = 0
    total_annotated = 0
    total_uniq_annotated = 0

    met = set()

    annotated_by_loc_by_gene = OrderedDefaultDict(lambda: defaultdict(list))
    off_targets = list()

    for fs in intersection:
        a_chr, a_start, a_end, e_chr, e_start, e_end, e_gene, e_exon, e_strand, \
            e_feature, e_biotype, overlap_size = fs
        assert e_chr == '.' or a_chr == e_chr, str((a_chr + ', ' + e_chr))
        total_lines += 1
        if (a_chr, a_start, a_end) not in met:
            total_uniq_lines += 1

        if e_chr == '.':
            off_targets.append(Region(a_chr, int(a_start), int(a_end)))
        else:
            total_annotated += 1
            if (a_chr, a_start, a_end) not in met:
                total_uniq_annotated += 1

            annotated_by_loc_by_gene[(a_chr, int(a_start), int(a_end))][e_gene].append((
                Region(chrom=e_chr, start=int(e_start), end=int(e_end), genome_build=cnf.genome.name,
                       gene_symbol=e_gene, exon=e_exon, strand=e_strand, feature=e_feature, biotype=e_biotype),
                int(overlap_size)))

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


def _split_reference(cnf, exons_bed_fpath):
    log('Splitting reference file into genes and non-genes:')
    exons_bed = BedTool(exons_bed_fpath)
    genes_exons_bed = exons_bed.filter(lambda x: x[6] in ['Gene'])
    no_genes_exons_bed = exons_bed.filter(lambda x: x[6] not in ['Gene', 'Multi_Gene'])

    # ref_bed_no_genes_fpath = join(cnf.work_dir, basename(exons_fpath) + '__no_whole_genes')
    # ref_bed_genes_fpath = join(cnf.work_dir, basename(exons_fpath) + '__whole_genes')
    # cmdline = 'grep -wv Gene {exons_fpath} | grep -wv Multi_Gene'.format(**locals())

    # if cnf.reuse_intermediate and verify_file(ref_bed_no_genes_fpath, silent=True, is_critical=False):
    #     log('Reusing existing ' + ref_bed_no_genes_fpath)
    # else:
    #     with open(ref_bed_no_genes_fpath, 'w') as out:
    #         cmdline = 'grep -wv Gene {ref_bed_fpath} | grep -wv Multi_Gene'.format(**locals())
    #         log(cmdline + ' > ' + ref_bed_no_genes_fpath)
    #         subprocess.call(cmdline, shell=True, stdout=out)
    # if verify_file(ref_bed_genes_fpath, silent=True, is_critical=False):
    #     log('Reusing existing ' + ref_bed_genes_fpath)
    # else:
    #     with open(ref_bed_genes_fpath, 'w') as out:
    #         cmdline = 'grep -w Gene {ref_bed_fpath}'.format(**locals())
    #         log(cmdline + ' > ' + ref_bed_genes_fpath)
    #         subprocess.call(cmdline, shell=True, stdout=out)
    # log()

    return genes_exons_bed, no_genes_exons_bed


if __name__ == '__main__':
    main()