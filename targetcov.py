#!/usr/bin/env python

from __future__ import print_function

from collections import defaultdict
from genericpath import isdir
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

_header = ['Sample', 'Chr', 'Start', 'End', 'Exome Num',
           'Direction', 'Length', 'MeanDepth']


def main(args):
    cnf, options = common_main(
        'targetcov',
        opts=[
            (['--bam'], 'align.bam', {
             'dest': 'bam',
             'help': 'used to generate some annotations by GATK'}),

            (['--capture', '--capture'], 'capture.bed', {
             'dest': 'capture',
             'help': ''}),

            (['--padding'], '250', {
             'dest': 'padding',
             'help': '',
             'default': 250}),
        ])

    genes_bed = cnf.get('genes')
    if not genes_bed:
        genes_bed = expanduser(cnf['genome'].get('genes'))
    if not verify_file(genes_bed):
        critical('Specify sorted genes bed file in system info or in run info.')

    exons_bed = cnf.get('exons')
    if not exons_bed:
        exons_bed = expanduser(cnf['genome'].get('exons'))
    if not verify_file(exons_bed):
        critical('Specify sorted exons bed file in system info or in run info.')

    chr_len_fpath = cnf.get('chr_lengths')
    if not chr_len_fpath:
        chr_len_fpath = expanduser(cnf['genome'].get('chr_lengths'))
    if not verify_file(chr_len_fpath):
        critical('Specify chromosome lengths for the genome'
                 ' in system info or in run info.')

    bam = expanduser(options.get('bam', cnf.get('bam')))
    if not bam:
        critical('Specify bam file by --bam option or in run_config.')
    if not verify_file(bam):
        exit()

    capture_bed = expanduser(options.get('capture', cnf.get('capture')))
    if not capture_bed:
        critical('Specify capture file by --capture option or in run_config.')
    if not verify_file(capture_bed):
        exit()

    depth_thresholds = cnf['depth_thresholds']
    padding = options.get('padding', cnf.get('padding', 250))
    output_dir = expanduser(options.get('output_dir', cnf.get('output_dir', os.getcwd())))

    work_dir = join(output_dir, 'work')
    if isdir(work_dir):
        rmtree(work_dir)
    os.makedirs(work_dir)

    run_header_report(output_dir, work_dir, capture_bed, bam, chr_len_fpath, depth_thresholds, padding)

    run_cov_exomes_report(output_dir, work_dir, capture_bed, bam, genes_bed, exons_bed, depth_thresholds)


# TODO check on division by 0
def run_header_report(output_dir, work_dir, capture_bed, bam, chr_len_fpath, depth_thresholds, padding=250):
    step_greetings('Target coverage summary report')

    base_name, ext = splitext(basename(bam))
    result_fpath = join(output_dir, base_name + '.' + 'summary.report')

    with file_transaction(result_fpath) as tx:
        with open(tx, 'w') as out:
            v_mapped_reads = int(number_of_mapped_reads(bam))
            out.write(format_integer('Number of mapped reads', v_mapped_reads))

            v_unmapped_reads = int(number_of_unmapped_reads(bam))
            out.write(format_integer('Number of unmapped reads', v_unmapped_reads))

            v_number_of_reads = int(number_of_reads(bam))
            out.write(format_integer('Number of reads', v_number_of_reads))
            v_percent_mapped = 100.0 * v_mapped_reads / v_number_of_reads if v_number_of_reads else None
            out.write(format_integer('Percent mapped reads', v_percent_mapped))

            v_reads_on_target = int(reads_on_targ(capture_bed, bam))
            out.write(format_integer('Number of reads on target', v_reads_on_target))
            vpercent_on_target = 100.0 * v_reads_on_target / v_mapped_reads if v_mapped_reads else None
            out.write(format_integer('Percent reads on target', vpercent_on_target))

            padded_bed = get_padded_bed_file(capture_bed, chr_len_fpath, padding, work_dir)
            v_padded_reads_on_targ = int(reads_on_targ(padded_bed, bam))
            v_percent_on_padded = 100.0 * v_padded_reads_on_targ / v_mapped_reads if v_mapped_reads else None
            out.write(format_integer('Percent reads on padded target', v_percent_on_padded))

            v_aligned_read_bases = total_aligned_base_reads(bam)
            out.write(format_integer('Total aligned base reads', v_aligned_read_bases))

            bases_per_depth, v_covd_ref_bases_on_targ, v_read_bases_on_targ, max_depth = \
                get_analytics_target_depth(capture_bed, bam, depth_thresholds)
            out.write(format_integer('Total base reads on target', v_read_bases_on_targ))

            v_percent_read_bases_on_targ = 100.0 * v_read_bases_on_targ / v_aligned_read_bases \
                if v_aligned_read_bases else None
            v_avg_cov_depth = v_read_bases_on_targ / v_covd_ref_bases_on_targ if v_covd_ref_bases_on_targ else None
            out.write(format_percent('Percent base reads on target', v_percent_read_bases_on_targ))
            out.write(format_integer('Bases in targeted reference ', v_covd_ref_bases_on_targ))
            # out.write(format_integer('Bases covered (at least 1x) ', target_groups[0][1]))
            out.write(format_percent('Average coverage depth ', v_avg_cov_depth))
            out.write(format_integer('Maximum read depth', max_depth))

            for depth, num in bases_per_depth.items():
                cov_at = 100.0 * num / v_covd_ref_bases_on_targ if v_covd_ref_bases_on_targ else '-'
                out.write(format_integer('Target coverage at ' + str(depth) + 'x', cov_at))

    print('')
    print('*' * 70)
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print('Done. Report: ' + result_fpath)


def run_cov_exomes_report(output_dir, work_dir, capture_bed, bam, genes_bed, exons_bed, depth_thresholds):
    step_greetings('Coverage report for exons in capture and genes')

    base_name, ext = splitext(basename(bam))
    result_fpath = join(output_dir, base_name + '.' + 'exons.report')

    capture_in_genes_bed = intersect_bed(genes_bed, capture_bed, work_dir)
    exome_in_capture_in_gene_bed = intersect_bed(exons_bed, capture_in_genes_bed, work_dir)

    with file_transaction(result_fpath) as tx:
        with open(exome_in_capture_in_gene_bed) as bed, \
             open(tx, 'w') as out:
            out.write('\t'.join(_header) + '\n')
            for bed_line in bed.readlines():
                write_report_line(bam, base_name, depth_thresholds, out, bed_line)

    print('')
    print('*' * 70)
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print('Done. Report: ' + result_fpath)


def _call(cmdline, output_fpath=None):
    print('')
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print(cmdline + (' > ' + output_fpath if output_fpath else ''))
    if output_fpath:
        with file_transaction(output_fpath) as tx:
            subprocess.call(cmdline.split(), stdout=open(tx, 'w'))


def _call_and_open_stdout(cmdline, stdout=subprocess.PIPE):
    print('')
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print(cmdline)
    return subprocess.Popen(cmdline.split(), stdout=stdout)


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
    return _call(cmdline, output_fpath)


def samtool_depth_range(bam_path, region):
    cmdline = 'samtools depth {bam_path} -r {region}'.format(**locals())
    return _call_and_open_stdout(cmdline, stdout=subprocess.PIPE)


def intersect_bed(bed_path_gene, bed_path_capture, work_dir):
    cmdline = 'intersectBed -a {bed_path_gene} -b {bed_path_capture} -u'.format(**locals())
    output_fpath = intermediate_fname(work_dir, bed_path_gene, 'intersect')
    _call(cmdline, output_fpath)
    return output_fpath


def number_of_mapped_reads(bam):
    cmdline = 'samtools view -c -F 4 {bam}'.format(**locals())
    res = _call_check_output(cmdline, stdout=subprocess.PIPE)
    return res


def number_of_unmapped_reads(bam):
    cmdline = 'samtools view -c -f 4 {bam}'.format(**locals())
    res = _call_check_output(cmdline, stdout=subprocess.PIPE)
    return res


def number_of_reads(bam):
    cmdline = 'samtools view -c {bam}'.format(**locals())
    res = _call_check_output(cmdline, stdout=subprocess.PIPE)
    return res


def reads_on_targ(bed, bam):
    cmdline = 'samtools view -c -F 4 {bam} -L {bed}'.format(**locals())
    res = _call_check_output(cmdline, stdout=subprocess.PIPE)
    return res


# TODO very slow
def total_aligned_base_reads(bam):
    cmdline = 'samtools depth {bam}'.format(**locals())
    proc = _call_and_open_stdout(cmdline, stdout=subprocess.PIPE)
    count = 0
    while True:
        coverage_line = proc.stdout.readline()
        if coverage_line != '':
            values = coverage_line.strip().split('\t')
            count += int(values[2])
        else:
            break
    return count


# TODO very slow too
def get_analytics_target_depth(bed, bam, depth_thresholds):
    cmdline = 'samtools depth {bam} -b {bed}'.format(**locals())
    proc = _call_and_open_stdout(cmdline, stdout=subprocess.PIPE)
    covered_bases = 0
    total_reads_on_target = 0
    max_depth = 0

    bases_per_depth = {depth: 0 for depth in depth_thresholds}

    for coverage_line in proc.stdout:
        if coverage_line != '':
            values = coverage_line.strip().split('\t')
            depth_value = int(values[2])
            total_reads_on_target += depth_value
            if max_depth < depth_value:
                max_depth = depth_value
            covered_bases += 1

            for depth in bases_per_depth.keys():
                if depth and depth_value >= bases_per_depth[depth]:
                    bases_per_depth[depth] += 1
        else:
            break

    return bases_per_depth, covered_bases, total_reads_on_target, max_depth


# TODO how to pass the data stream to samtools vs. creating file
def get_padded_bed_file(bed, genome, padding, work_dir):
    cmdline = 'bedtools slop -i {bed} -g {genome} -b {padding}'.format(**locals())
    output_fpath = intermediate_fname(work_dir, bed, 'padded')
    return _call(cmdline, output_fpath)


def count_by_group(depth_values, depth_thresholds):
    bases_by_min_depth = {depth: 0 for depth in depth_thresholds}

    for depth_value in depth_values:
        for threshold in depth_thresholds:
            if depth_value >= threshold:
                bases_by_min_depth[threshold] += 1

    return [float(100 * k[1]) / len(depth_values) for k in bases_by_min_depth]


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


def format_integer(name, value):
    if value:
        return '{0}: {1}\n'.format(name, value)
    else:
        return '{0}: -\n'.format(name)


def format_percent(name, value):
    if value:
        return '{0}: {1:.2f}\n'.format(name, value)
    else:
        return '{0}: -\n'.format(name)


def mean(ints):
    float(sum(ints))/len(ints) if len(ints) > 0 else float('nan')


def write_report_line(bam, base_name, depth_thresholds, out, bed_line):
    values = bed_line.strip().split('\t')
    chr = values[0]
    start = values[1]
    end = values[2]
    region_for_samtools = '{0}:{1}-{2}'.format(chr, start, end)
    depths = get_depth_by_bed_range(bam, region_for_samtools)
    length = (int(end) - int(start)) + 1

    if len(depths) > 0:
        depth_mean = '{0:.2f}'.format(mean(depths))
        counts = count_by_group(depths, depth_thresholds)
        out.write('\t'.join([base_name, bed_line.strip(),
                         length, depth_mean] + map(str, counts)) + '\n')
    else:
        out.write('\t'.join([base_name, bed_line.strip(), length, 'na']))
# '\t'.join([sample_name, line, length, mean, start_end_count] + start_end_count)


def run_cov_report(output_dir, work_dir, bed, bam, depth_thresholds):
    step_greetings('Coverage report')

    base_name, ext = os.path.splitext(bam)
    output_path = base_name + '.' + 'report'

    print('Creating report path: ' + output_path)

    bed_sorted_path = gnu_sort(bed, work_dir)

    with open(bed_sorted_path) as sorted_bed, open(output_path, 'w') as out:
        out.write('\t'.join(_header) + '\n')
        for bed_line in sorted_bed.readlines():
            write_report_line(bam, base_name, depth_thresholds, out, bed_line)


# def run_cov_gene_report(gene_bed, capture_bed, bam, depth_thresholds, work_dir):
#     base_name, ext = os.path.splitext(bam)
#     output_path = base_name + '.' + 'gene.report'
#
#     # a_bed_sorted = gnu_sort(a_bed)
#     # b_bed_sorted = gnu_sort(b_bed)
#     int_bed = intersect_bed(gene_bed, capture_bed, work_dir)
#
#     with open(int_bed) as bed, open(output_path, 'w') as out:
#         out.write('\t'.join(_header) + '\n')
#         for bed_line in bed.readlines():
#             write_report_line(bam, base_name, depth_thresholds, out, bed_line)


if __name__ == '__main__':
    main(sys.argv)