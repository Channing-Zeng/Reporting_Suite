#!/usr/bin/env python

from __future__ import print_function

from collections import defaultdict
import numpy

import sys
import os
from os.path import join, splitext, basename, realpath
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
# take folder name as a sample name ( first column on the report)
# give user an option to select type of the report to run ????
from src.main import common_main

_header = ['Sample', 'Chr', 'Start', 'End', 'Exome Num',
           'Direction', 'Length', 'MeanDepth']


def main(args):
    config, options = common_main(
        'targetcov',
        opts=[
            (['--bam'], 'align.bam', {
             'dest': ['bam'],
             'help': 'used to generate some annotations by GATK'}),

            (['--capture', '--capture'], 'capture.bed', {
             'dest': ['capture'],
             'help': ''}),

            (['--padding'], '250', {
             'dest': ['padding'],
             'help': '',
             'default': 250}),
        ])

    genes_bed = config['genome']['genes']
    exons_bed = config['genome']['exons']
    chr_len_fpath = config['genome']['chr_lengths']
    depth_thresholds = config['depth_thresholds']

    bam = options.bam
    capture = options.capture

    run_header_report(capture, bam, chr_len_fpath, depth_thresholds, padding=250)

    run_cov_exomes_report(capture, bam, genes_bed, exons_bed, depth_thresholds)


def _call(cmdline, stdout=None):
    print('')
    print('*' * 70)
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print(cmdline)
    subprocess.call(cmdline.split(), stdout=stdout)


def _call_stdout(cmdline, stdout=subprocess.PIPE):
    print('')
    print('*' * 70)
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print(cmdline)
    return subprocess.Popen(cmdline.split(), stdout=stdout)


# TODO to check if input files a re not empty
def _call_and_write(cmdline, fpath, new_ext):
    base_name, ext = os.path.splitext(fpath)
    output_fpath = base_name + '.' + new_ext
    if not os.path.isfile(output_fpath):
        _call(cmdline, open(output_fpath, 'w'))
        #TODO check if we have file
    return output_fpath


def gnu_sort(bed_path):
    cmdline = 'sort -k1,1V  -k2,2n -k3,3n {bed_path}'.format(**locals())
    return _call_and_write(cmdline, bed_path, 'sorted.bed')


def samtool_depth_range(bam_path, region):
    cmdline = 'samtools depth {bam_path} -r {region}'.format(**locals())
    return _call_stdout(cmdline, stdout=subprocess.PIPE)


def intersect_bed(bed_path_gene, bed_path_capture):
    cmdline = 'intersectBed -a {bed_path_gene} -b {bed_path_capture} -u'.format(**locals())
    return _call_and_write(cmdline, bed_path_gene, 'intersect.bed')


def number_of_mapped_reads(bam):
    cmdline = 'samtools view -c -F 4 {bam}'
    proc = _call_stdout(cmdline, stdout=subprocess.PIPE)
    return proc.stdout.readline()


def number_of_unmapped_reads(bam):
    cmdline = 'samtools view -c -f 4 %s' % ( bam)
    proc = _call_stdout(cmdline, stdout=subprocess.PIPE)
    return proc.stdout.readline()


def number_of_reads(bam):
    cmdline = 'samtools view -c %s' % bam
    proc = _call_stdout(cmdline, stdout=subprocess.PIPE)
    return proc.stdout.readline()


def reads_on_target(bed, bam):
    cmdline = 'samtools view -c -F 4 %s -L %s' % (bam, bed)
    proc = _call_stdout(cmdline, stdout=subprocess.PIPE)
    return proc.stdout.readline()


# TODO very slow
def total_aligned_base_reads(bam):
    cmdline = 'samtools depth %s ' % bam
    proc = _call_stdout(cmdline, stdout=subprocess.PIPE)
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
def get_analytics_target_depth(bed, bam):
    cmdline = 'samtools depth {bam} -b {bed}'.format(**locals())
    proc = _call_stdout(cmdline, stdout=subprocess.PIPE)
    covered_bases = 0
    total_reads_on_target = 0
    max_depth = 0

    bases_per_depth = {depth: 0 for depth in _depth_thresholds}

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
def get_padded_bed_file(bed, genome, padding):
    cmdline = 'bedtools slop -i {bed} -g {genome} -b {padding}'.format(**locals())
    return _call_and_write(cmdline, bed, 'padded.bed')


def count_by_group(depth_values):
    bases_by_min_depth = {depth: 0 for depth in _depth_thresholds}

    for depth_value in depth_values:
        for threshold in _depth_thresholds:
            if depth_value >= threshold:
                bases_by_min_depth[threshold] += 1

    return [float(100 * k[1]) / len(depth_values) for k in bases_by_min_depth]
# print('%.2f' % a)

def get_depth_by_bed_range(bam, region):
    proc = samtool_depth_range(bam, region)
    depths = []
    for coverage_line in proc.stdout:
        values = coverage_line.strip().split('\t')
        if coverage_line != '' and len(values) > 2:
            depths.append(values[2])
        else:
            break

    return depths


def print_name_value(name, value):
    return '{0}: {1}\n'.format(name, value)


def write_report_line(bam, base_name, out, bed_line):
    values = bed_line.strip().split('\t')
    chr = values[0]
    start = values[1]
    end = values[2]
    region_for_samtools = '{0}:{1}-{2}'.format(chr, start, end)
    depths = get_depth_by_bed_range(bam, region_for_samtools)
    length = (int(end) - int(start)) + 1
    if len(depths) > 0:
        depths_int = map(int, depths)
        depth_mean = '{0:.2f}'.format(numpy.mean(depths_int))
        counts = count_by_group(depths_int)
        print('\t'.join([base_name, bed_line.strip(), length, depth_mean] + map(str, counts)), file=out)
    else:
        out.write(print_line(base_name, bed_line.strip(), length, 'na', []))
# '\t'.join([sample_name, line, length, mean, start_end_count] + start_end_count)


# TODO check on division by 0
def run_header_report(bed, bam, genome, padding, depth_thresholds):
    print('')
    print('*' * 70)
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print('Printing header report')
    base_name, ext = os.path.splitext(bam)
    output_path = base_name + '.' + 'header.report'
    with open(output_path, 'w') as out:
        vnumber_of_mapped_reads = int(number_of_mapped_reads(bam))
        out.write(print_name_value('Number of mapped reads', vnumber_of_mapped_reads))

        vnumber_of_unmapped_reads = int(number_of_unmapped_reads(bam))
        out.write(print_name_value('Number of unmapped reads', vnumber_of_unmapped_reads))

        vnumber_of_reads = int(number_of_reads(bam))
        out.write(print_name_value('Number of reads', vnumber_of_reads))
        out.write(print_name_value('Percent mapped reads', 100 * vnumber_of_mapped_reads / vnumber_of_reads))

        vnumber_reads_on_target = int(reads_on_target(bed, bam))
        out.write(print_name_value('Number of reads on target', vnumber_reads_on_target))

        percent_reads_on_target = 100 * vnumber_reads_on_target / vnumber_of_mapped_reads
        out.write(print_name_value('Percent reads on target', percent_reads_on_target))

        vnumber_padded_reads_on_target = int(reads_on_target(get_padded_bed_file(bed, genome, padding), bam))
        out.write(print_name_value('Percent reads on padded target',
                                   100 * vnumber_padded_reads_on_target / vnumber_of_mapped_reads))

        vtotal_aligned_base_reads = total_aligned_base_reads(bam)
        out.write(print_name_value('Total aligned base reads', vtotal_aligned_base_reads))

        target_groups, cov, vtotal_aligned_base_reads_on_target, max_depth = get_analytics_target_depth(bed, bam)
        out.write(print_name_value('Total base reads on target', vtotal_aligned_base_reads_on_target))

        percent_based_reads_on_target = 100 * vtotal_aligned_base_reads_on_target / vtotal_aligned_base_reads
        out.write(print_name_value('Percent base reads on target', percent_based_reads_on_target))
        out.write(print_name_value('Bases in targeted reference ', cov))
        out.write(print_name_value('Bases covered (at least 1x) ', target_groups[0][1]))
        out.write(print_name_value('Average base coverage depth ', vtotal_aligned_base_reads_on_target / cov))
        out.write(print_name_value('Maximum base read depth', max_depth))

        #Std.Dev base read depth
        [out.write(print_name_value('Target coverage at ' + str(group[0]) + 'x', 100 * group[1] / cov)) for group in
         target_groups]

    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print('Done with header report: ' + output_path)


def run_cov_report(bed, bam):
    print('')
    print('*' * 70)
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    base_name, ext = os.path.splitext(bam)
    output_path = base_name + '.' + 'report'

    print('Creating report path: ' + output_path)

    bed_sorted_path = gnu_sort(bed)

    with open(bed_sorted_path) as sorted_bed, open(output_path, 'w') as out:
        out.write('\t'.join(_header) + '\n')
        for bed_line in sorted_bed.readlines():
            write_report_line(bam, base_name, out, bed_line)


def run_cov_gene_report(gene_bed, capture_bed, bam):
    base_name, ext = os.path.splitext(bam)
    output_path = base_name + '.' + 'gene.report'

    # a_bed_sorted = gnu_sort(a_bed)
    # b_bed_sorted = gnu_sort(b_bed)
    int_bed = intersect_bed(gene_bed, capture_bed)

    with open(int_bed) as bed, open(output_path, 'w') as out:
        out.write('\t'.join(_header) + '\n')
        for bed_line in bed.readlines():
            write_report_line(bam, base_name, out, bed_line)


def run_cov_exomes_report(gene_bed, exome_bed, capture_bed, bam):
    base_name, ext = os.path.splitext(bam)
    output_path = base_name + '.' + 'exomes.report'

    int_gene_bed = intersect_bed(gene_bed, capture_bed)
    int_gene_exone_bed = intersect_bed(exome_bed, int_gene_bed)

    with open(int_gene_exone_bed) as  bed, open(output_path, 'w') as out:
        out.write('\t'.join(_header) + '\n')
        for bed_line in bed.readlines():
            write_report_line(bam, base_name, out, bed_line)

# _gene_exome,bed_path_capture,bam_path)