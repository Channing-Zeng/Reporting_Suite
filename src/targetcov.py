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
from src.my_utils import step_greetings, intermediate_fname, timestamp
from src.transaction import file_transaction
from src.utils import splitext_plus


def run_header_report(output_dir, work_dir, capture_bed, bam, chr_len_fpath, depth_thresholds, padding=250):
    step_greetings('Target coverage summary report')

    base_name, ext = splitext(basename(bam))
    result_fpath = join(output_dir, base_name + '.' + 'summary.report')

    padded_bed = get_padded_bed_file(capture_bed, chr_len_fpath, padding, work_dir)

    # v_mapped_reads = number_of_mapped_reads(bam)
    # v_unmapped_reads = number_of_unmapped_reads(bam)
    # v_number_of_reads = number_of_reads(bam)
    # v_percent_mapped = 100.0 * v_mapped_reads / v_number_of_reads if v_number_of_reads else None
    # v_reads_on_target = number_reads_on_target(capture_bed, bam)
    # v_percent_on_target = 100.0 * v_reads_on_target / v_mapped_reads if v_mapped_reads else None
    # v_reads_on_padded_targ = number_reads_on_target(padded_bed, bam)
    # v_percent_on_padded = 100.0 * v_reads_on_padded_targ / v_mapped_reads if v_mapped_reads else None
    # v_aligned_read_bases = number_bases_in_aligned_reads(bam)

    bases_per_depth, v_covd_ref_bases_on_targ, max_depth, \
        bases_per_depth_per_region = \
            get_target_depth_analytics_fast(capture_bed, bam, depth_thresholds)

    # v_percent_read_bases_on_targ = 100.0 * v_read_bases_on_targ / v_aligned_read_bases \
    #     if v_aligned_read_bases else None
    # v_avg_cov_depth = float(v_read_bases_on_targ) / v_covd_ref_bases_on_targ \
    #     if v_covd_ref_bases_on_targ else None

    # stats = [format_integer('Number of mapped reads', v_mapped_reads),
    #          format_integer('Number of unmapped reads', v_unmapped_reads),
    #          format_integer('Number of reads', v_number_of_reads),
    #          format_decimal('Percent mapped reads', v_percent_mapped, '%'),
    #          format_integer('Number of reads on target', v_reads_on_target),
    #          format_decimal('Percent reads on target', v_percent_on_target, '%'),
    #          format_decimal('Percent reads on padded target', v_percent_on_padded, '%'),
    #          # format_integer('Total aligned bases in reads', v_aligned_read_bases),
    #          # format_integer('Total bases in reads on target', v_read_bases_on_targ),
    #          # format_decimal('Percent bases in reads on target', v_percent_read_bases_on_targ, '%'),
    #          format_integer('Bases in targeted reference', v_covd_ref_bases_on_targ),
    #          format_integer('Bases covered (at least 1x)', bases_per_depth[1]),
    #          # format_decimal('Average coverage depth', v_avg_cov_depth),
    #          format_integer('Maximum read depth', max_depth)]

    # for depth, num in bases_per_depth.items():
    #     covd_at = 100.0 * num / v_covd_ref_bases_on_targ if v_covd_ref_bases_on_targ else 0
    #     stats.append(format_decimal('Target covered at ' + str(depth) + 'x', covd_at, '%'))

    # max_len = max(len(l.rsplit(':', 1)[0]) for l in stats)
    # with file_transaction(result_fpath) as tx, open(tx, 'w') as out:
    #     for l in stats:
    #         text, val = l.rsplit(':', 1)
    #         spaces = ' ' * (max_len - len(text) + 1)
    #         out.write(text + spaces + val + '\n')

    print('')
    print('*' * 70)
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print('Done. Report: ' + result_fpath)

    return bases_per_depth_per_region


def run_cov_report(output_dir, work_dir, capture_bed, bam, depth_threshs,
                   bases_per_depth_per_region, genes_bed=None, exons_bed=None):
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

    header, max_lengths = None, None

    all_values = []
    # with open(bed) as bed_f:
    #     for line in (l.strip() for l in bed_f if l and l.strip()):
    #         # line_vals = get_report_line_values(bam, sample_name, depth_threshs, line)

    required_fields_start = ['#Sample', 'Chr', 'Start', 'End']
    required_fields_end = ['Size', 'AvgDepth'] + map(str, depth_threshs)

    for i, (region, (bp_per_depths, avg_depth)) in enumerate(bases_per_depth_per_region.items()):
        # for depth, num in bases_per_depth.items():
        #     covd_at = 100.0 * num / v_covd_ref_bases_on_targ if v_covd_ref_bases_on_targ else 0
        #     line_vals.append(covd_at + '%'))
        region_tokens = region.split()  # Chr, Start, End, RegionSize, F1, F2,...
        region_size = region_tokens[-1]
        region_tokens = region_tokens[:-1]
        if i == 0:
            header = (required_fields_start +
                      ([''] * (len(region_tokens) - len(required_fields_start))) +
                      required_fields_end)
            max_lengths = map(len, header)

        # avg_depth_str = '{0:.2f}'.format(mean(depths) if depths else 0.0)

        line_tokens = [sample_name] + region_tokens
        line_tokens += ['-'] * (len(header) - len(line_tokens) - len(required_fields_end))
        line_tokens += [region_size]
        line_tokens += ['{0:.2f}'.format(avg_depth)]

        for depth_thres, bases in bp_per_depths.items():
            if int(region_size) == 0:
                percent_str = '-'
            else:
                percent = 100.0 * float(bases) / int(region_size)
                percent_str = '{0:.2f}%'.format(percent) if percent != 0 else '0'
            line_tokens.append(percent_str)

        all_values.append(line_tokens)
        max_lengths = map(max, izip(max_lengths, chain(map(len, line_tokens), repeat(0))))


    with file_transaction(out_fpath) as tx:
        with open(tx, 'w') as out:
            for h, l in zip(header, max_lengths):
                sys.stdout.write(h + ' ' * (l - len(h) + 2))
                out.write(h + ' ' * (l - len(h) + 2))
            out.write('\n')

            for line_tokens in all_values:
                for v, l in zip(line_tokens, max_lengths):
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
#


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


def samtool_depth_range(bam_path, region):
    cmdline = 'samtools depth -r {region} {bam_path}'.format(**locals())
    return _call_and_open_stdout(cmdline)


#
# def samtool_depth_range_fast(bam_path, region):
#     cmdline = 'coverageBed -abam {bam} -b {bed} -hist'.format(**locals())
#     return _call_and_open_stdout(cmdline)

# chr17   62006585        62006684        NM_000626_cds_0_0_chr17_62006586_r      0       -       10      4       99      0.0404040
# chr17   62006585        62006684        NM_000626_cds_0_0_chr17_62006586_r      0       -       11      14      99      0.1414141
# chr17   62006585        62006684        NM_000626_cds_0_0_chr17_62006586_r      0       -       12      2       99      0.0202020
# chr17   62006585        62006684        NM_000626_cds_0_0_chr17_62006586_r      0       -       13      2       99      0.0202020
# chr17   62006585        62006684        NM_000626_cds_0_0_chr17_62006586_r      0       -       14      5       99      0.0505050
# chr17   62006585        62006684        NM_000626_cds_0_0_chr17_62006586_r      0       -       15      17      99      0.1717172
# chr17   62006585        62006684        NM_000626_cds_0_0_chr17_62006586_r      0       -       16      8       99      0.0808081
# chr17   62006585        62006684        NM_000626_cds_0_0_chr17_62006586_r      0       -       17      19      99      0.1919192
# chr17   62006585        62006684        NM_000626_cds_0_0_chr17_62006586_r      0       -       18      25      99      0.2525252
# chr17   62006585        62006684        NM_000626_cds_0_0_chr17_62006586_r      0       -       19      3       99      0.0303030
# chr17   62006793        62006835        NM_000626_cds_1_0_chr17_62006794_r      0       -       5       12      42      0.2857143
# chr17   62006793        62006835        NM_000626_cds_1_0_chr17_62006794_r      0       -       6       30      42      0.7142857
# chr17   62007129        62007248        NM_000626_cds_2_0_chr17_62007130_r      0       -       7       15      119     0.1260504
# chr17   62007129        62007248        NM_000626_cds_2_0_chr17_62007130_r      0       -       8       5       119     0.0420168
# chr17   62007129        62007248        NM_000626_cds_2_0_chr17_62007130_r      0       -       9       5       119     0.0420168
# chr17   62007129        62007248        NM_000626_cds_2_0_chr17_62007130_r      0       -       10      44      119     0.3697479
# chr17   62007129        62007248        NM_000626_cds_2_0_chr17_62007130_r      0       -       11      11      119     0.0924370
# chr17   62007129        62007248        NM_000626_cds_2_0_chr17_62007130_r      0       -       12      26      119     0.2184874
# chr17   62007129        62007248        NM_000626_cds_2_0_chr17_62007130_r      0       -       13      6       119     0.0504202
# chr17   62007129        62007248        NM_000626_cds_2_0_chr17_62007130_r      0       -       14      7       119     0.0588235
# chr17   62007433        62007745        NM_000626_cds_3_0_chr17_62007434_r      0       -       3       10      312     0.0320513
# chr17   62007433        62007745        NM_000626_cds_3_0_chr17_62007434_r      0       -       4       10      312     0.0320513

def get_target_depth_analytics_fast(bed, bam, depth_thresholds):
    cmdline = 'coverageBed -abam {bam} -b {bed} -hist'.format(**locals())
    proc = _call_and_open_stdout(cmdline)
    covered_bases = 0
    max_depth = 0

    bases_per_depth_all = OrderedDict([(depth_thres, 0) for depth_thres in depth_thresholds])
    percent_per_depth_all = OrderedDict([(depth_thres, 0.0) for depth_thres in depth_thresholds])

    bases_per_depth_per_region = OrderedDict()
    # percent_per_depth_per_region = OrderedDict()

    _prev_region_line = None
    _prev_region_size = 0
    _prev_region_depth_sum = 0

    regions_number = 0

    for line in proc.stdout:
        if not line.strip() or line.startswith('#'):
            continue

        def set_avg_depth_for_prev_region():
            if _prev_region_size and _prev_region_line:
                avg_depth = _prev_region_depth_sum / _prev_region_size
                bases_per_depth_per_region[_prev_region_line][1] = avg_depth

        tokens = line.strip().split('\t')
        depth, bases_for_depth, region_size = map(int, tokens[-4:-1])
        region_tokens = [str(region_size)]

        if line.startswith('all'):
            set_avg_depth_for_prev_region()

            region_tokens = [tokens[0]] + region_tokens

            max_depth = max(max_depth, depth)
            covered_bases += bases_for_depth

            for depth_thres in depth_thresholds:
                if depth >= depth_thres:
                    bases_per_depth_all[depth_thres] += bases_for_depth
        else:
            region_tokens = tokens[:-4] + region_tokens

        region_line = '\t'.join(region_tokens)
        if region_line not in bases_per_depth_per_region:
            set_avg_depth_for_prev_region()

            regions_number += 1

            bases_per_depth_per_region[region_line] = [bases_per_depth_all.copy(), None]
            _prev_region_line = region_line
            _prev_region_depth_sum = 0
            _prev_region_size = region_size

        _prev_region_depth_sum += depth * bases_for_depth

        for depth_thres in depth_thresholds:
            if depth >= depth_thres:
                bases_per_depth_per_region[region_line][0][depth_thres] += float(bases_for_depth)

    print(timestamp() + ' Processed ' + str(regions_number) + ' regions.')

    return bases_per_depth_all, covered_bases, max_depth, \
           bases_per_depth_per_region


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