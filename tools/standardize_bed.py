#!/usr/bin/env python

from os.path import abspath, dirname, realpath, join, exists, basename
from site import addsitedir
project_dir = abspath(dirname(dirname(realpath(__file__))))
addsitedir(join(project_dir))
addsitedir(join(project_dir, 'ext_modules'))
import sub_scripts.__check_python_version  # do not remove it: checking for python version and adding site dirs inside

import sys
import os
import annotate_bed
import subprocess
from source.file_utils import add_suffix, _remove_files
from source.utils import human_sorted


""" Input: Any BED file
    Output:
        BED file with regions from input and with exactly 4 columns (or exactly 8 columns for BEDs with primers info).
        The forth column ("gene symbol") is generated with annotate_bed.py script. If no annotation is found
        gene symbol is set to "pseudo_gene_%d". If more than one annotation is found for region only one is remained
        (either from key gene list or just first one).
        Output BED is sorted using "sort -k1,1V -k2,2n -k3,3n", GRCh37 chr names remain the same, header remain the same.

    Usage: python standardize_bed Input_BED_file work_dir
                  [Key_genes_file] [Reference_BED_file] [bedtools_tool_path]
                  > standardized_BED_file
"""

VERBOSE = False
CLEANUP = True


def _read_args(args):
    if len(args) < 2:
        log('Usage:\n')
        log('  ' + __file__ + ' Input_BED_file work_dir [Key_genes_file] [Reference_BED_file] [bedtools_tool_path] > standardized_BED_file\n')
        sys.exit(1)

    input_bed_fpath = abspath(args[0])
    log('Input: ' + input_bed_fpath)

    work_dirpath = abspath(args[1])
    log('Working directory: ' + work_dirpath)
    if not exists(work_dirpath):
        os.makedirs(work_dirpath)

    key_genes_fpath = '/gpfs/ngs/oncology/Analysis/dev/Dev_0075_PanelExomeFMComparison/bed/az_key_genes.txt'
    if len(args) > 2:
        if args[2]:
            key_genes_fpath = abspath(args[2])
            log('Over-set key genes fpath: ' + key_genes_fpath)

    ref_bed_fpath = '/ngs/reference_data/genomes/Hsapiens/hg19/bed/Exons/Exons.with_genes.bed'
    if len(args) > 3:
        if args[3]:
            ref_bed_fpath = abspath(args[3])
            log('Over-set reference fpath: ' + ref_bed_fpath)

    bedtools = 'bedtools'
    if len(args) > 4:
        bedtools = args[4]
        log('Over-set bedtools: ' + bedtools)
    log()

    return input_bed_fpath, work_dirpath, key_genes_fpath, ref_bed_fpath, bedtools


def log(msg=''):
    sys.stderr.write(msg + '\n')


def err(msg=''):
    log('Error: ' + msg + '\n')
    sys.exit(1)


class BedParams:
    GRCh_to_hg = {'MT': 'chrM', 'X': 'chrX', 'Y': 'chrY'}
    for i in range(1, 23):
        GRCh_to_hg[str(i)] = 'chr' + str(i)
    hg_to_GRCh = {v: k for k, v in GRCh_to_hg.items()}

    def __init__(self, header=list(), GRCh_names=None, n_cols_needed=None):
        self.header = header
        self.GRCh_names = GRCh_names
        self.n_cols_needed = n_cols_needed

    @staticmethod
    def calc_n_cols_needed(line):
        columns = line.strip().split('\t')
        if len(columns) < 8:
            return 4
        elif columns[6].isdigit() and columns[7].isdigit():  # primers info
            return 8
        else:
            return 4


class Region:
    GRCh_names = False
    n_cols_needed = 4

    def __init__(self, bed_line):
        entries = bed_line.strip().split('\t')
        self.chrom = entries[0]
        self.start = int(entries[1])
        self.end = int(entries[2])
        self.symbol = entries[3] if len(entries) > 3 else '{0}:{1}-{2}'.format(self.chrom, self.start, self.end)
        self.rest = entries[4:]

    def __str__(self):
        fs = [BedParams.hg_to_GRCh[self.chrom] if self.GRCh_names else self.chrom,
              str(self.start), str(self.end), self.symbol] + (self.rest if self.n_cols_needed > 4 else [])
        return '\t'.join(fs) + '\n'

    def get_key(self):
        return '\t'.join([self.chrom, str(self.start), str(self.end)])

    def __lt__(self, other):
        # special case: chrM goes to the end
        if self.chrom != other.chrom and (self.chrom == 'chrM' or other.chrom == 'chrM'):
            return True if other.chrom == 'chrM' else False
        sorted_pair = human_sorted([self.get_key(), other.get_key()])
        if sorted_pair[0] == self.get_key() and sorted_pair[1] != self.get_key():
            return True
        return False

    def __eq__(self, other):
        return self.get_key() == other.get_key()

    def __hash__(self):
        return hash(self.get_key())


def _preprocess(bed_fpath, work_dirpath):
    bed_params = BedParams()
    output_fpath = __intermediate_fname(work_dirpath, bed_fpath, 'prep')
    log('preprocessing: ' + bed_fpath + ' --> ' + output_fpath)
    with open(bed_fpath, 'r') as in_f:
        with open(output_fpath, 'w') as out_f:
            for line in in_f:
                if line.startswith('#') or line.startswith('track') or line.startswith('browser'):  # header
                    bed_params.header.append(line)
                else:
                    cur_ncn = BedParams.calc_n_cols_needed(line)
                    if bed_params.n_cols_needed is not None and cur_ncn != bed_params.n_cols_needed:
                        err('number and type of columns should be the same on all lines!')
                    bed_params.n_cols_needed = cur_ncn
                    if line.startswith('chr'):
                        if bed_params.GRCh_names is not None and bed_params.GRCh_names:
                            err('mixing of GRCh and hg chromosome names!')
                        bed_params.GRCh_names = False
                        out_f.write(line)
                    elif line.split('\t')[0] in BedParams.GRCh_to_hg:  # GRCh chr names
                        if bed_params.GRCh_names is not None and not bed_params.GRCh_names:
                            err('mixing of GRCh and hg chromosome names!')
                        bed_params.GRCh_names = True
                        out_f.write('\t'.join([BedParams.GRCh_to_hg[line.split('\t')[0]]] + line.split('\t')[1:]))
                    else:
                        err('incorrect chromosome name!')
    return output_fpath, bed_params


def _annotate(bed_fpath, work_dirpath, ref_bed_fpath, bedtools):
    output_fpath = __intermediate_fname(work_dirpath, bed_fpath, 'ann')
    log('annotating: ' + bed_fpath + ' --> ' + output_fpath)
    annotate_bed_py = sys.executable + ' ' + annotate_bed.__file__
    cmdline = '{annotate_bed_py} {bed_fpath} {work_dirpath} {ref_bed_fpath} {bedtools}'.format(**locals())
    __call(cmdline, output_fpath)
    return output_fpath


def _postprocess(input_fpath, annotated_fpath, bed_params, key_genes_fpath):
    '''
    1. Sorts.
    1. Chooses appropriate number of columns (4 or 8 for BEDs with primers).
    2. Removes duplicates.
    '''
    log('postprocessing (sorting, cutting, removing duplicates)')

    key_genes = []
    with open(key_genes_fpath, 'r') as f:
        for line in f:
            key_genes.append(line.strip())

    input_regions = set()  # we want only unique regions
    with open(input_fpath) as f:
        for line in f:
            input_regions.add(Region(line))
    annotated_regions = []
    with open(annotated_fpath) as f:
        for line in f:
            annotated_regions.append(Region(line))

    # starting to output result
    for line in bed_params.header:
        sys.stdout.write(line)
    Region.GRCh_names = bed_params.GRCh_names
    Region.n_cols_needed = bed_params.n_cols_needed

    annotated_regions.sort()
    i = 0
    prev_symbol, prev_chr, pseudo_count = "", "", 0
    for in_region in sorted(input_regions):
        ready_region = in_region
        assert annotated_regions[i] == in_region, str(in_region) + ' != ' + str(annotated_regions[i]) + '(i=%d)' % i
        if annotated_regions[i].symbol != '.':
            ready_region.symbol = annotated_regions[i].symbol
        else:
            if prev_chr != ready_region.chrom or not prev_symbol.startswith("pseudo_gene"):
                pseudo_count += 1
            ready_region.symbol = "pseudo_gene_%d" % pseudo_count
        i += 1
        while i < len(annotated_regions) and annotated_regions[i] == in_region:  # processing duplicates
            if annotated_regions[i].symbol != '.' and annotated_regions[i].symbol in key_genes:
                ready_region.symbol = annotated_regions[i].symbol
            i += 1
        sys.stdout.write(str(ready_region))  # automatically output correct number of columns and GRCh/hg names
        prev_chr = ready_region.chrom
        prev_symbol = ready_region.symbol

def __intermediate_fname(work_dir, fname, suf):
    output_fname = add_suffix(fname, suf)
    return join(work_dir, basename(output_fname))


def __call(cmdline, output_fpath=None):
    stdout = open(output_fpath, 'w') if output_fpath else None
    stderr = None if VERBOSE else open('/dev/null', 'w')
    if VERBOSE:
        log(cmdline)
    ret_code = subprocess.call(cmdline, shell=True, stdout=stdout, stderr=stderr, stdin=None)
    return ret_code


def __sort_bed(fpath):
    work_dirpath = dirname(fpath)
    output_fpath = __intermediate_fname(work_dirpath, fpath, 'sort')
    cmdline = 'sort -k1,1V -k2,2n -k3,3n ' + fpath  # TODO: -k1,1V not working on Mac
    __call(cmdline, output_fpath)
    return fpath


def main():
    input_bed_fpath, work_dirpath, key_genes_fpath, ref_bed_fpath, bedtools = _read_args(sys.argv[1:])

    preprocessed_fpath, bed_params = _preprocess(input_bed_fpath, work_dirpath)
    annotated_fpath = _annotate(preprocessed_fpath, work_dirpath, ref_bed_fpath, bedtools)
    _postprocess(preprocessed_fpath, annotated_fpath, bed_params, key_genes_fpath)
    if CLEANUP:
        _remove_files([preprocessed_fpath, annotated_fpath])

if __name__ == '__main__':
    main()
