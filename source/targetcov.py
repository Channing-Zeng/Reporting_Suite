#!/usr/bin/env python

from __future__ import print_function

from collections import OrderedDict
from itertools import izip, chain, repeat

import sys
from os.path import join, splitext, basename
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
import math
from source.my_utils import step_greetings, intermediate_fname, timestamp
from source.transaction import file_transaction
from source.utils import splitext_plus


def log(msg=''):
    print(timestamp() + msg)


def run_header_report(output_dir, work_dir, bed, bam, chr_len_fpath,
                      depth_thresholds, padding,
                      bases_per_depth, avg_depth, std_dev, bases_within_normal_deviation,
                      max_depth, total_bed_size):
    step_greetings('Target coverage summary report')

    base_name, ext = splitext(basename(bam))
    result_fpath = join(output_dir, base_name + '.' + 'summary.report')

    stats = []
    def append_stat(stat):
        stats.append(stat)
        log(stat)

    log('* General coverage statistics *')
    log('Getting number of reads...')
    v_number_of_reads = number_of_reads(bam)
    append_stat(format_integer('Reads', v_number_of_reads))

    log('Getting number of mapped reads...')
    v_mapped_reads = number_of_mapped_reads(bam)
    append_stat(format_integer('Mapped reads', v_mapped_reads))
    append_stat(format_integer('Unmapped reads', v_number_of_reads - v_mapped_reads))

    v_percent_mapped = 100.0 * v_mapped_reads / v_number_of_reads if v_number_of_reads else None
    append_stat(format_decimal('Percentage of mapped reads', v_percent_mapped, '%'))
    log('')

    log('* Target coverage statistics *')
    append_stat(format_integer('Bases in target', total_bed_size))

    v_covered_bases_in_targ = bases_per_depth.items()[0][1]
    append_stat(format_integer('Covered bases in target', v_covered_bases_in_targ))

    v_percent_covered_bases_in_targ = 100.0 * v_covered_bases_in_targ / total_bed_size if total_bed_size else None
    append_stat(format_decimal('Percentage of target covered by at least 1 read', v_percent_covered_bases_in_targ, '%'))
    log('Getting number of mapped reads on target...')

    v_mapped_reads_on_target = number_mapped_reads_on_target(bed, bam)
    append_stat(format_integer('Reads mapped on target', v_mapped_reads_on_target))

    v_percent_mapped_on_target = 100.0 * v_mapped_reads_on_target / v_mapped_reads if v_mapped_reads else None
    append_stat(format_decimal('Percentage of reads mapped on target ', v_percent_mapped_on_target, '%'))

    log('Making bed file for padded regions...')
    padded_bed = get_padded_bed_file(bed, chr_len_fpath, padding, work_dir)

    log('Getting number of mapped reads on padded target...')
    v_reads_on_padded_targ = number_mapped_reads_on_target(padded_bed, bam)
    append_stat(format_integer('Reads mapped on padded target', v_reads_on_padded_targ, '%'))

    v_percent_mapped_on_padded_target = 100.0 * v_reads_on_padded_targ / v_mapped_reads if v_mapped_reads else None
    append_stat(format_decimal('Percentage of reads mapped on padded target', v_percent_mapped_on_padded_target, '%'))

    # v_aligned_read_bases = number_bases_in_aligned_reads(bam)
    # append_stat(format_integer('Total aligned bases in reads', v_aligned_read_bases))

    v_read_bases_on_targ = avg_depth * total_bed_size  # sum of all coverages
    append_stat(format_integer('Read bases mapped on target', v_read_bases_on_targ))

    log('')
    percent_within_normal_deviation = bases_within_normal_deviation / total_bed_size
    append_stat(format_decimal('Average target coverage depth', avg_depth))
    append_stat(format_decimal('Std. dev. of target coverage depth', std_dev))
    append_stat(format_integer('Maximum target coverage depth', max_depth))
    append_stat(format_decimal('Percentage thing 20% of avarage depth in a region',
                               percent_within_normal_deviation, '%'))

    # v_percent_read_bases_on_targ = 100.0 * v_read_bases_on_targ / v_aligned_read_bases \
    #     if v_aligned_read_bases else None
    # format_decimal('Percent bases in reads on target', v_percent_read_bases_on_targ, '%'),

    # format_integer('Bases covered (at least 1x) in target', bases_per_depth[1]),

    # for depth, bases in bases_per_depth.items():
    #     append_stat(format_integer('Bases on target covered at least by ' + str(depth) +
    #                                ' read' + ('s' if depth != 1 else ''), bases))

    for depth, bases in bases_per_depth.items():
        percent = 100.0 * bases / total_bed_size if total_bed_size else 0
        append_stat(format_decimal('Part of target covered at least by ' + str(depth) +
                                   ('x'), percent, '%'))

    max_len = max(len(l.rsplit(':', 1)[0]) for l in stats)
    with file_transaction(result_fpath) as tx, open(tx, 'w') as out:
        for l in stats:
            text, val = l.rsplit(':', 1)
            spaces = ' ' * (max_len - len(text) + 1)
            out.write(text + spaces + val + '\n')

    log('')
    log('Result: ' + result_fpath)
    return result_fpath


def run_cov_report(output_dir, work_dir, bed, bam, depth_threshs,
                   bases_per_depth_per_region):
    sample_name, _ = splitext(basename(bam))

    out_fpath = join(output_dir, sample_name + '.exons.report')

    # bed_sorted_path = gnu_sort(capture_bed, work_dir)

    header, max_lengths = None, None

    all_values = []
    # with open(bed) as bed_f:
    #     for line in (l.strip() for l in bed_f if l and l.strip()):
    #         # line_vals = get_report_line_values(bam, sample_name, depth_threshs, line)

    required_fields_start = ['#Sample', 'Chr', 'Start', 'End', 'Transcript', 'Gene', 'ExonNum', 'Strand']
    required_fields_end = ['Size', 'AvgDepth', 'StdDev', 'Within20%'] + map(str, depth_threshs)

    for i, (region, (bp_per_depths, avg_depth, std_dev, bases_within_normal_deviation)) \
            in enumerate(bases_per_depth_per_region.items()):
        # for depth, num in bases_per_depth.items():
        #     covd_at = 100.0 * num / v_covd_ref_bases_on_targ if v_covd_ref_bases_on_targ else 0
        #     line_vals.append(covd_at + '%'))
        region_tokens = region.split()  # Chr, Start, End, RegionSize, F1, F2,...
        region_size = int(region_tokens[-1])
        region_tokens = region_tokens[:-1]
        if i == 0:
            header = (required_fields_start +
                      (['-'] * (len(region_tokens) - len(required_fields_start) + 1)) +
                      required_fields_end)
            max_lengths = map(len, header)

        # avg_depth_str = '{0:.2f}'.format(mean(depths) if depths else 0.0)

        line_tokens = [sample_name] + region_tokens
        line_tokens += ['-'] * (len(header) - len(line_tokens) - len(required_fields_end))
        line_tokens += [str(region_size)]
        line_tokens += ['{0:.2f}'.format(avg_depth) if avg_depth else '0']
        line_tokens += ['{0:.2f}'.format(std_dev) if std_dev else '0']

        percent_within_normal_dev = 100.0 * bases_within_normal_deviation / region_size if region_size else '-'
        line_tokens += ['{0:.2f}%'.format(percent_within_normal_dev) if percent_within_normal_dev else '0']

        for depth_thres, bases in bp_per_depths.items():
            if int(region_size) == 0:
                percent_str = '-'
            else:
                percent = 100.0 * float(bases) / region_size
                percent_str = '{0:.2f}%'.format(percent) if percent != 0 else '0'
            line_tokens.append(percent_str)

        all_values.append(line_tokens)
        max_lengths = map(max, izip(max_lengths, chain(map(len, line_tokens), repeat(0))))

    with file_transaction(out_fpath) as tx:
        with open(tx, 'w') as out:
            for h, l in zip(header, max_lengths):
                out.write(h + ' ' * (l - len(h) + 2))
            out.write('\n')

            for line_tokens in all_values:
                for v, l in zip(line_tokens, max_lengths):
                    out.write(v + ' ' * (l - len(v) + 2))
                out.write('\n')

    log('')
    log('Result: ' + out_fpath)
    return out_fpath


def _call(cmdline, output_fpath=None):
    log('  $ ' + cmdline + (' > ' + output_fpath if output_fpath else ''))
    if output_fpath:
        with file_transaction(output_fpath) as tx:
            subprocess.call(cmdline.split(), stdout=open(tx, 'w'))


def _call_and_open_stdout(cmdline):
    log('  $ ' + cmdline)
    return subprocess.Popen(cmdline.split(), stdout=subprocess.PIPE)


def _call_check_output(cmdline, stdout=subprocess.PIPE):
    log('  $ ' + cmdline)
    return subprocess.check_output(cmdline.split())


def samtool_depth_range(bam_path, region):
    cmdline = 'samtools depth -r {region} {bam_path}'.format(**locals())
    return _call_and_open_stdout(cmdline)


def mapped_bases_in_bed_using_samtoolsdepth(bam, bed):
    total = 0
    for st, end in (l.split()[1:3] for l in open(bed).readlines() if l.strip()):
        total += len([l for l in samtool_depth_range(bam, 'chrM:' + str(int(st) + 1) + '-' + end).stdout
                   if l.strip() and not l.startswith('#')])
    return total


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


# def total_bed_length(bed_fpath):
#     cmdline = 'cat {bed} | awk -F"\t" ' \
#               '"BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}"'.format(**locals())
#     return int(_call_check_output(cmdline))


def intersect_bed(bed1, bed2, work_dir):
    bed1_fname, _ = splitext_plus(basename(bed1))
    bed2_fname, _ = splitext_plus(basename(bed2))
    output_fpath = join(work_dir, bed1_fname + '__' + bed2_fname + '.bed')
    cmdline = 'intersectBed -a {bed1} -b {bed2} -u'.format(**locals())
    _call(cmdline, output_fpath)
    return output_fpath


# TODO very slow :(
def number_of_mapped_reads(bam):
    cmdline = 'samtools view -c -F 4 {bam}'.format(**locals())
    res = _call_check_output(cmdline)
    return int(res)


# TODO very slow :(
def number_of_unmapped_reads(bam):
    cmdline = 'samtools view -c -f 4 {bam}'.format(**locals())
    res = _call_check_output(cmdline)
    return int(res)


# TODO very slow :(
def number_of_reads(bam):
    cmdline = 'samtools view -c {bam}'.format(**locals())
    res = _call_check_output(cmdline)
    return int(res)


# TODO very slow :(
def number_mapped_reads_on_target(bed, bam):
    cmdline = 'samtools view -c -F 4 -L {bed} {bam}'.format(**locals())
    res = _call_check_output(cmdline)
    return int(res)


# TODO very slow :(
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
    max_depth = 0
    total_size = 0

    bases_per_depth_per_region = OrderedDict()

    _prev_region_bases_by_depths = []
    _prev_region_line = None
    _prev_region_size = 0

    regions_number = 0

    for line in proc.stdout:
        if not line.strip() or line.startswith('#'):
            continue

        def set_avg_depth_for_prev_region():
            if _prev_region_size and _prev_region_line:
                depth_sum = sum(depth * bases
                                for depth, bases
                                in _prev_region_bases_by_depths)
                avg_depth = depth_sum / _prev_region_size

                sum_of_sq_var = sum((depth - avg_depth)**2 * bases
                                    for depth, bases
                                    in _prev_region_bases_by_depths)
                std_dev = math.sqrt(sum_of_sq_var / _prev_region_size)

                bases_within_normal_deviation = sum(
                    bases for depth, bases
                    in _prev_region_bases_by_depths
                    if math.fabs(avg_depth - depth) <= 0.2 * avg_depth)

                bases_per_depth_per_region[_prev_region_line][1] = avg_depth
                bases_per_depth_per_region[_prev_region_line][2] = std_dev
                bases_per_depth_per_region[_prev_region_line][3] = bases_within_normal_deviation

        tokens = line.strip().split('\t')
        depth, bases_for_depth, region_size = map(int, tokens[-4:-1])
        region_tokens = [str(region_size)]

        if line.startswith('all'):
            set_avg_depth_for_prev_region()
            region_tokens = [tokens[0]] + region_tokens
            max_depth = max(max_depth, depth)
            total_size += bases_for_depth

            # for depth_thres in depth_thresholds:
            #     if depth >= depth_thres:
            #         bases_per_depth_all[depth_thres] += bases_for_depth
        else:
            region_tokens = tokens[:-4] + region_tokens

        region_line = '\t'.join(region_tokens)
        if region_line not in bases_per_depth_per_region:
            set_avg_depth_for_prev_region()
            regions_number += 1

            bases_per_depth_per_region[region_line] = \
                [OrderedDict([(depth_thres, 0.0) for depth_thres
                              in depth_thresholds]),
                 None,
                 None,
                 None]
            _prev_region_bases_by_depths = []
            _prev_region_line = region_line
            _prev_region_size = region_size

        _prev_region_bases_by_depths.append((depth, bases_for_depth))
        # _prev_region_sum_of_sq_var += bases_for_depth * depth -

        for depth_thres in depth_thresholds:
            if depth >= depth_thres:
                bases_per_depth_per_region[region_line][0][depth_thres] += float(bases_for_depth)

        if regions_number > 0 and regions_number % 100000 == 0:
            log('processed %i regions' % regions_number)

    if regions_number % 100000 != 0:
        log('processed %i regions' % regions_number)

    return bases_per_depth_per_region, max_depth, total_size


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


# def get_depth_by_bed_range(bam, region):
#     proc = samtool_depth_range(bam, region)
#     depths = []
#     for coverage_line in proc.stdout:
#         if coverage_line.strip():
#             values = coverage_line.strip().split('\t')
#             if len(values) > 2:
#                 depths.append(int(values[2]))
#             else:
#                 break
#     return depths


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


# def get_report_line_values(bam, sample_name, depth_threshs, bed_line):
#     chrom, start, end = bed_line.strip().split('\t')[:3]
#     region_for_samtools = '{chrom}:{start}-{end}'.format(**locals())
#     depths = get_depth_by_bed_range(bam, region_for_samtools)
#     region_size = (int(end) - int(start)) + 1
#
#     vals = [sample_name, chrom, start, end, str(region_size)]
#
#     avg_depth_str = '{0:.2f}'.format(mean(depths) if depths else 0.0)
#     by_depth = bps_by_depth(depths, depth_threshs)
#     by_depth_str = ['{0:.2f}'.format(c) for c in by_depth]
#     return vals + [avg_depth_str] + by_depth_str