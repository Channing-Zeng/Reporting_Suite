#!/usr/bin/env python

from __future__ import print_function

from collections import defaultdict, OrderedDict
from genericpath import isdir
from itertools import izip, chain, repeat
import numpy

import sys
import os
from os.path import join, splitext, basename, realpath, expanduser
import subprocess
from datetime import datetime
#downlad hg19.genome
#https://github.com/arq5x/bedtools/tree/master/genomes

#TODO
# check on the input file format
# format result and calculation to the 2 decimal places on the header report    .00
# multi - sample report                                                         header only
#       sample1 sample2
#number 2       3
#bases  10      20

# check if samtools and bedtools exist
# log file
# yaml
# take folder name as a sample name (first column on the report)
# give user an option to select type of the report to run ????
from shutil import rmtree
from src.main import common_main
from src.my_utils import verify_file, critical, step_greetings, intermediate_fname
from src.transaction import file_transaction
from src.utils import splitext_plus


def main(args):
    cnf, options = common_main(
        'targetcov',
        opts=[
            (['--bam'], 'align.bam', {
                'dest': 'bam',
                'help': 'used to generate some annotations by GATK'}),

            (['--capture', '--bed'], 'capture.bed', {
                'dest': 'capture',
                'help': ''}),

            (['--genes', '--genes'], 'genes.bed', {
                'dest': 'genes',
                'help': ''}),

            (['--exons', '--exons'], 'exons.bed', {
                'dest': 'exons',
                'help': ''}),

            (['--padding'], '250', {
                'dest': 'padding',
                'help': '',
                'default': 250}),
        ])

    genes_bed = options.get('genes') or cnf.get('genes') or cnf['genome'].get('genes')
    exons_bed = options.get('exons') or cnf.get('exons') or expanduser(cnf['genome'].get('exons'))
    chr_len_fpath = cnf.get('chr_lengths') or cnf['genome'].get('chr_lengths')
    capture_bed = options.get('capture') or cnf.get('capture')
    bam = options.get('bam') or cnf.get('bam')

    genes_bed = expanduser(genes_bed)
    exons_bed = expanduser(exons_bed)
    chr_len_fpath = expanduser(chr_len_fpath)
    bam = expanduser(bam)
    capture_bed = expanduser(capture_bed)

    if not genes_bed:
        critical('Specify sorted genes bed file in system info or in run info.')
    if not exons_bed:
        critical('Specify sorted exons bed file in system info or in run info.')
    if not chr_len_fpath:
        critical('Specify chromosome lengths for the genome'
                 ' in system info or in run info.')
    if not bam:
        critical('Specify bam file by --bam option or in run_config.')
    if not capture_bed:
        critical('Specify capture file by --capture option or in run_config.')

    print('using genes ' + genes_bed)
    print('using exons ' + exons_bed)
    print('using chr lengths ' + chr_len_fpath)
    print('using bam ' + bam)
    print('using capture panel ' + capture_bed)

    if not verify_file(genes_bed): exit(1)
    if not verify_file(exons_bed): exit(1)
    if not verify_file(chr_len_fpath): exit(1)
    if not verify_file(bam): exit(1)
    if not verify_file(capture_bed): exit(1)

    depth_thresholds = cnf['depth_thresholds']
    padding = options.get('padding', cnf.get('padding', 250))
    output_dir = options.get('output_dir') or cnf.get('output_dir') or os.getcwd()
    assert output_dir
    output_dir = expanduser(output_dir)

    work_dir = join(output_dir, 'work')
    if isdir(work_dir):
        rmtree(work_dir)
    os.makedirs(work_dir)

    run_header_report(output_dir, work_dir, capture_bed, bam, chr_len_fpath, depth_thresholds, padding)

    run_cov_report(output_dir, work_dir, capture_bed, bam, depth_thresholds)

    run_cov_report(output_dir, work_dir, capture_bed, bam, depth_thresholds, genes_bed, exons_bed)


# TODO check on division by 0
def run_header_report(output_dir, work_dir, capture_bed, bam, chr_len_fpath, depth_thresholds, padding=250):
    step_greetings('Target coverage summary report')

    base_name, ext = splitext(basename(bam))
    result_fpath = join(output_dir, base_name + '.' + 'summary.report')

    padded_bed = get_padded_bed_file(capture_bed, chr_len_fpath, padding, work_dir)

    v_mapped_reads = number_of_mapped_reads(bam)
    v_unmapped_reads = number_of_unmapped_reads(bam)
    v_number_of_reads = number_of_reads(bam)
    v_percent_mapped = 100.0 * v_mapped_reads / v_number_of_reads if v_number_of_reads else None
    v_reads_on_target = number_reads_on_target(capture_bed, bam)
    v_percent_on_target = 100.0 * v_reads_on_target / v_mapped_reads if v_mapped_reads else None
    v_reads_on_padded_targ = number_reads_on_target(padded_bed, bam)
    v_percent_on_padded = 100.0 * v_reads_on_padded_targ / v_mapped_reads if v_mapped_reads else None
    # v_aligned_read_bases = number_bases_in_aligned_reads(bam)

    bases_per_depth, v_covd_ref_bases_on_targ, max_depth = \
        get_target_depth_analytics_fast(capture_bed, bam, depth_thresholds)

    # v_percent_read_bases_on_targ = 100.0 * v_read_bases_on_targ / v_aligned_read_bases \
    #     if v_aligned_read_bases else None
    # v_avg_cov_depth = float(v_read_bases_on_targ) / v_covd_ref_bases_on_targ \
    #     if v_covd_ref_bases_on_targ else None

    stats = [format_integer('Number of mapped reads', v_mapped_reads),
             format_integer('Number of unmapped reads', v_unmapped_reads),
             format_integer('Number of reads', v_number_of_reads),
             format_decimal('Percent mapped reads', v_percent_mapped, '%'),
             format_integer('Number of reads on target', v_reads_on_target),
             format_decimal('Percent reads on target', v_percent_on_target, '%'),
             format_decimal('Percent reads on padded target', v_percent_on_padded, '%'),
             # format_integer('Total aligned bases in reads', v_aligned_read_bases),
             # format_integer('Total bases in reads on target', v_read_bases_on_targ),
             # format_decimal('Percent bases in reads on target', v_percent_read_bases_on_targ, '%'),
             format_integer('Bases in targeted reference', v_covd_ref_bases_on_targ),
             format_integer('Bases covered (at least 1x)', bases_per_depth[1]),
             # format_decimal('Average coverage depth', v_avg_cov_depth),
             format_integer('Maximum read depth', max_depth)]

    for depth, num in bases_per_depth.items():
        covd_at = 100.0 * num / v_covd_ref_bases_on_targ if v_covd_ref_bases_on_targ else 0
        stats.append(format_decimal('Target covered at ' + str(depth) + 'x', covd_at, '%'))

    max_len = max(len(l.rsplit(':', 1)[0]) for l in stats)
    with file_transaction(result_fpath) as tx, open(tx, 'w') as out:
        for l in stats:
            text, val = l.rsplit(':', 1)
            spaces = ' ' * (max_len - len(text) + 1)
            out.write(text + spaces + val + '\n')

    print('')
    print('*' * 70)
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print('Done. Report: ' + result_fpath)


def run_cov_report(output_dir, work_dir, capture_bed, bam, depth_threshs, genes_bed=None, exons_bed=None):
    sample_name, _ = splitext(basename(bam))

    if genes_bed and exons_bed:
        step_greetings('Coverage report for exons in capture and genes')
        out_fpath = join(output_dir, sample_name + '.exome_cov.report')
    else:
        step_greetings('Coverage report')
        out_fpath = join(output_dir, sample_name + '.cov.report')

    # bed_sorted_path = gnu_sort(capture_bed, work_dir)

    bed = capture_bed
    if genes_bed:
        bed = intersect_bed(genes_bed, bed, work_dir)
    if exons_bed:
        bed = intersect_bed(exons_bed, bed, work_dir)

    header = (['Sample', 'Chr', 'Start', 'End', 'RegionSize', 'MeanDepth']
              + map(str, depth_threshs))
    max_lengths = map(len, header)

    all_values = []
    with open(bed) as bed_f:
        for line in (l.strip() for l in bed_f if l and l.strip()):
            line_vals = get_report_line_values(bam, sample_name, depth_threshs, line)
            all_values.append(line_vals)
            max_lengths = map(max, izip(max_lengths,
                                        chain(map(len, line_vals), repeat(0))))

    with file_transaction(out_fpath) as tx:
        with open(tx, 'w') as out:
            for h, l in zip(header, max_lengths):
                sys.stdout.write(h + ' ' * (l - len(h) + 2))
                out.write(h + ' ' * (l - len(h) + 2))
            out.write('\n')

            for line_vals in all_values:
                for v, l in zip(line_vals, max_lengths):
                    out.write(v + ' ' * (l - len(v) + 2))
                out.write('\n')
    print('')
    print('*' * 70)
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print('Done. Report: ' + out_fpath)


def _call(cmdline, output_fpath=None):
    print('')
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print(cmdline + (' > ' + output_fpath if output_fpath else ''))
    if output_fpath:
        with file_transaction(output_fpath) as tx:
            subprocess.call(cmdline.split(), stdout=open(tx, 'w'))


def _call_and_open_stdout(cmdline):
    print('')
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print(cmdline)
    return subprocess.Popen(cmdline.split(), stdout=subprocess.PIPE)


def _call_check_output(cmdline, stdout=subprocess.PIPE):
    print('')
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print(cmdline)
    return subprocess.check_output(cmdline.split())


# # TODO to check if input files a re not empty
# def _call_and_write(cmdline, fpath, new_ext):
#     base_name, ext = os.path.splitext(fpath)
#     output_fpath = base_name + '.' + new_ext
#     if not os.path.isfile(output_fpath):
#         _call(cmdline, open(output_fpath, 'w'))
#         #TODO check if we have file
#     return output_fpath


def gnu_sort(bed_path, work_dir):
    cmdline = 'sort -k1,1V -k2,2n -k3,3n {bed_path}'.format(**locals())
    output_fpath = intermediate_fname(work_dir, bed_path, 'sorted')
    _call(cmdline, output_fpath)
    return output_fpath


def samtool_depth_range(bam_path, region):
    cmdline = 'samtools depth -r {region} {bam_path}'.format(**locals())
    return _call_and_open_stdout(cmdline)


def intersect_bed(bed1, bed2, work_dir):
    bed1_fname, _ = splitext_plus(basename(bed1))
    bed2_fname, _ = splitext_plus(basename(bed2))
    output_fpath = join(work_dir, bed1_fname + '__' + bed2_fname + '.bed')
    cmdline = 'intersectBed -a {bed1} -b {bed2} -u'.format(**locals())
    _call(cmdline, output_fpath)
    return output_fpath


def number_of_mapped_reads(bam):
    cmdline = 'samtools view -c -F 4 {bam}'.format(**locals())
    res = _call_check_output(cmdline)
    return int(res)


def number_of_unmapped_reads(bam):
    cmdline = 'samtools view -c -f 4 {bam}'.format(**locals())
    res = _call_check_output(cmdline)
    return int(res)


def number_of_reads(bam):
    cmdline = 'samtools view -c {bam}'.format(**locals())
    res = _call_check_output(cmdline)
    return int(res)


def number_reads_on_target(bed, bam):
    cmdline = 'samtools view -c -F 4 -L {bed} {bam}'.format(**locals())
    res = _call_check_output(cmdline)
    return int(res)


# TODO very slow
def number_bases_in_aligned_reads(bam):
    cmdline = 'samtools depth {bam}'.format(**locals())
    proc = _call_and_open_stdout(cmdline)
    count = 0
    while True:
        coverage_line = proc.stdout.readline()
        if coverage_line:
            values = coverage_line.strip().split('\t')
            count += int(values[2])
    return count


def get_target_depth_analytics_fast(bed, bam, depth_thresholds):
    cmdline = 'coverageBed -abam {bam} -b {bed} -hist'.format(**locals())
    proc = _call_and_open_stdout(cmdline)
    covered_bases = 0
    max_depth = 0

    bases_per_depth = OrderedDict([(depth_thres, 0) for depth_thres in depth_thresholds])

    for line in proc.stdout:
        if line and line.startswith('all'):
            _, depth, bases, reg_size, percent_on_depth = line.strip().split('\t')
            # print(line)
            depth = int(depth)
            bases = int(bases)
            reg_size = int(reg_size)

            max_depth = max(max_depth, depth)
            covered_bases += bases

            for depth_thres in depth_thresholds:
                if depth >= depth_thres:
                    bases_per_depth[depth_thres] += bases

    return bases_per_depth, covered_bases, max_depth


# TODO how to pass the data stream to samtools vs. creating file
def get_padded_bed_file(bed, genome, padding, work_dir):
    cmdline = 'bedtools slop -i {bed} -g {genome} -b {padding}'.format(**locals())
    output_fpath = intermediate_fname(work_dir, bed, 'padded')
    _call(cmdline, output_fpath)
    return output_fpath


def bps_by_depth(depth_vals, depth_thresholds):
    bases_by_min_depth = {depth: 0 for depth in depth_thresholds}

    for depth_value in depth_vals:
        for threshold in depth_thresholds:
            if depth_value >= threshold:
                bases_by_min_depth[threshold] += 1

    return [100.0 * bases_by_min_depth[thres] / len(depth_vals) if depth_vals else 0
            for thres in depth_thresholds]


def get_depth_by_bed_range(bam, region):
    proc = samtool_depth_range(bam, region)
    depths = []
    for coverage_line in proc.stdout:
        if coverage_line.strip():
            values = coverage_line.strip().split('\t')
            if len(values) > 2:
                depths.append(int(values[2]))
            else:
                break
    return depths


def format_integer(name, value, unit=''):
    if value is not None:
        return '{name}: {value:,}{unit}'.format(**locals())
    else:
        return '{name}: -'.format(**locals())


def format_decimal(name, value, unit=''):
    if value is not None:
        return '{name}: {value:.2f}{unit}'.format(**locals())
    else:
        return '{name}: -'.format(**locals())


def mean(ints):
    return float(sum(ints)) / len(ints) if len(ints) > 0 else float('nan')


def get_report_line_values(bam, sample_name, depth_threshs, bed_line):
    chrom, start, end = bed_line.strip().split('\t')[:3]
    region_for_samtools = '{chrom}:{start}-{end}'.format(**locals())
    depths = get_depth_by_bed_range(bam, region_for_samtools)
    region_size = (int(end) - int(start)) + 1

    vals = [sample_name, chrom, start, end, str(region_size)]

    avg_depth_str = '{0:.2f}'.format(mean(depths) if depths else 0.0)
    by_depth = bps_by_depth(depths, depth_threshs)
    by_depth_str = ['{0:.2f}'.format(c) for c in by_depth]
    return vals + [avg_depth_str] + by_depth_str


if __name__ == '__main__':
    main(sys.argv)