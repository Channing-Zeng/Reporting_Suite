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
from source.utils import step_greetings, intermediate_fname, timestamp, get_tool_cmdline
from source.transaction import file_transaction
from source.bcbio_utils import splitext_plus


def log(msg=''):
    print(timestamp() + msg)


def run_header_report(cnf, result_fpath, output_dir, work_dir,
                      bed, bam, chr_len_fpath,
                      depth_thresholds, padding,
                      all_region, max_depth, total_bed_size):
    stats = []

    def append_stat(stat):
        stats.append(stat)
        log(stat)

    log('* General coverage statistics *')
    log('Getting number of reads...')
    v_number_of_reads = number_of_reads(cnf, bam)
    append_stat(format_integer('Reads', v_number_of_reads))

    log('Getting number of mapped reads...')
    v_mapped_reads = number_of_mapped_reads(cnf, bam)
    append_stat(format_integer('Mapped reads', v_mapped_reads))
    append_stat(format_integer('Unmapped reads', v_number_of_reads - v_mapped_reads))

    v_percent_mapped = 100.0 * v_mapped_reads / v_number_of_reads if v_number_of_reads else None
    append_stat(format_decimal('Percentage of mapped reads', v_percent_mapped, '%'))
    log('')

    log('* Target coverage statistics *')
    append_stat(format_integer('Bases in target', total_bed_size))

    v_covered_bases_in_targ = all_region.bases_by_depth.items()[0][1]
    append_stat(format_integer('Covered bases in target', v_covered_bases_in_targ))

    v_percent_covered_bases_in_targ = 100.0 * v_covered_bases_in_targ / total_bed_size if total_bed_size else None
    append_stat(format_decimal('Percentage of target covered by at least 1 read', v_percent_covered_bases_in_targ, '%'))
    log('Getting number of mapped reads on target...')

    v_mapped_reads_on_target = number_mapped_reads_on_target(cnf, bed, bam)
    append_stat(format_integer('Reads mapped on target', v_mapped_reads_on_target))

    v_percent_mapped_on_target = 100.0 * v_mapped_reads_on_target / v_mapped_reads if v_mapped_reads else None
    append_stat(format_decimal('Percentage of reads mapped on target ', v_percent_mapped_on_target, '%'))

    log('Making bed file for padded regions...')
    padded_bed = get_padded_bed_file(cnf, bed, chr_len_fpath, padding, work_dir)

    log('Getting number of mapped reads on padded target...')
    v_reads_on_padded_targ = number_mapped_reads_on_target(cnf, padded_bed, bam)
    append_stat(format_integer('Reads mapped on padded target', v_reads_on_padded_targ, '%'))

    v_percent_mapped_on_padded_target = 100.0 * v_reads_on_padded_targ / v_mapped_reads if v_mapped_reads else None
    append_stat(format_decimal('Percentage of reads mapped on padded target', v_percent_mapped_on_padded_target, '%'))

    # v_aligned_read_bases = number_bases_in_aligned_reads(bam)
    # append_stat(format_integer('Total aligned bases in reads', v_aligned_read_bases))

    v_read_bases_on_targ = all_region.avg_depth * total_bed_size  # sum of all coverages
    append_stat(format_integer('Read bases mapped on target', v_read_bases_on_targ))

    log('')
    append_stat(format_decimal('Average target coverage depth', all_region.avg_depth))
    append_stat(format_decimal('Std. dev. of target coverage depth', all_region.std_dev))
    append_stat(format_integer('Maximum target coverage depth', max_depth))
    append_stat(format_decimal('Percentage thing 20% of avarage depth in a region',
                               all_region.percent_within_normal, '%'))

    # v_percent_read_bases_on_targ = 100.0 * v_read_bases_on_targ / v_aligned_read_bases \
    #     if v_aligned_read_bases else None
    # format_decimal('Percent bases in reads on target', v_percent_read_bases_on_targ, '%'),

    # format_integer('Bases covered (at least 1x) in target', bases_per_depth[1]),

    # for depth, bases in bases_per_depth.items():
    #     append_stat(format_integer('Bases on target covered at least by ' + str(depth) +
    #                                ' read' + ('s' if depth != 1 else ''), bases))

    for depth, bases in all_region.bases_by_depth.items():
        percent = 100.0 * bases / total_bed_size if total_bed_size else 0
        append_stat(format_decimal('Part of target covered at least by ' + str(depth) +
                                   'x', percent, '%'))

    max_len = max(len(l.rsplit(':', 1)[0]) for l in stats)
    with file_transaction(result_fpath) as tx, open(tx, 'w') as out:
        for l in stats:
            text, val = l.rsplit(':', 1)
            spaces = ' ' * (max_len - len(text) + 1)
            out.write(text + spaces + val + '\n')

    log('')
    log('Result: ' + result_fpath)
    return result_fpath


def run_amplicons_cov_report(cnf, report_fpath, depth_threshs,
                             bases_per_depth_per_amplicon):
    return run_cov_report(cnf, report_fpath, depth_threshs,
                          bases_per_depth_per_amplicon)


def run_exons_cov_report(cnf, report_fpath, depth_threshs,
                         bases_per_depth_per_exon_and_gene):
    return run_cov_report(cnf, report_fpath, depth_threshs,
                          bases_per_depth_per_exon_and_gene,
                          extra_fields=['Transcript', 'Gene'])


def add_genes_cov_analytics(bases_per_depth_per_exon):
    bases_per_depth_per_exon_and_gene = []

    current_genename = None
    current_gene = None
    current_geneexons = []

    def sum_up_gene(gene, gene_exons):
        gene.end = gene_exons[-1].end
        gene.avg_depth = 0
        gene.percent_within_normal = 0
        gene.std_dev = 0
        bases_per_depth_per_exon_and_gene.append(current_gene)

    for exon in bases_per_depth_per_exon:
        if len(exon.extra_fields) < 2:
           sys.exit('no gene info in exons record: ' + str(exon))
        exon_genename = exon.extra_fields[1]

        if exon_genename == current_genename:
            current_geneexons.append(exon)

        else:
            if current_gene and current_geneexons:
                sum_up_gene(current_gene, current_geneexons)

            current_geneexons.append(exon)

            current_genename = exon_genename
            current_geneexons = [exon]
            current_gene = Region(
                sample=exon.sample, chrom=exon.chrom,
                start=exon.start, feature='Gene',
                extra_fields=[current_genename, '.'],
                depth_thresholds=exon.depth_thresholds)

        bases_per_depth_per_exon_and_gene.append(exon)

    if current_gene and current_geneexons:
        sum_up_gene(current_gene, current_geneexons)

    return bases_per_depth_per_exon_and_gene


def run_cov_report(cnf, report_fpath, depth_threshs,
                   bases_per_depth_per_region, extra_fields=list()):
    first_fields = ['SAMPLE', 'Chr', 'Start', 'End', 'Feature']
    last_fields = ['Size', 'Mean Depth', 'Standard Dev.', 'Within 20% of Mean']
    header = first_fields + extra_fields + last_fields
    header += ['{}x'.format(thres) for thres in depth_threshs]
    max_lengths = map(len, header)

    all_values = []

    for i, region in enumerate(bases_per_depth_per_region):
        line_tokens = map(str, [region.sample,
                                region.chrom,
                                region.start,
                                region.end,
                                region.feature])
        line_tokens += region.extra_fields[:len(extra_fields)]
        line_tokens += [str(region.get_size())]
        line_tokens += ['{0:.2f}'.format(region.avg_depth)]
        line_tokens += ['{0:.2f}'.format(region.std_dev)]
        line_tokens += ['{0:.2f}%'.format(region.percent_within_normal)
                        if region.percent_within_normal is not None else '-']

        for depth_thres, bases in region.bases_by_depth.items():
            if int(region.get_size()) == 0:
                percent_str = '-'
            else:
                percent = 100.0 * bases / region.get_size()
                percent_str = '{0:.2f}%'.format(percent)
            line_tokens.append(percent_str)

        all_values.append(line_tokens)
        max_lengths = map(max, izip(max_lengths, chain(map(len, line_tokens), repeat(0))))

    with file_transaction(report_fpath) as tx:
        with open(tx, 'w') as out:
            for h, l in zip(header, max_lengths):
                out.write(h + ' ' * (l - len(h) + 2))
            out.write('\n')

            for line_tokens in all_values:
                for v, l in zip(line_tokens, max_lengths):
                    out.write(v + ' ' * (l - len(v) + 2))
                out.write('\n')

    log('')
    log('Result: ' + report_fpath)
    return report_fpath


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


# def samtool_depth_range(cnf, bam_path, region):
#     bedtools = get_tool_cmdline(cnf, 'samtools')
#     cmdline = '{samtools} depth -r {region} {bam_path}'.format(**locals())
#     return _call_and_open_stdout(cmdline)
#
#
# def mapped_bases_in_bed_using_samtoolsdepth(cnf, bam, bed):
#     total = 0
#     for st, end in (l.split()[1:3] for l in open(bed).readlines() if l.strip()):
#         total += len([l for l in samtool_depth_range(cnf, bam, 'chrM:' + str(int(st) + 1) + '-' + end).stdout
#                    if l.strip() and not l.startswith('#')])
#     return total


# # TODO to check if input files a re not empty
# def _call_and_write(cmdline, fpath, new_ext):
#     base_name, ext = os.path.splitext(fpath)
#     output_fpath = base_name + '.' + new_ext
#     if not os.path.isfile(output_fpath):
#         _call(cmdline, open(output_fpath, 'w'))
#         #TODO check if we have file
#     return output_fpath


def gnu_sort(cnf, bed_path, work_dir):
    sort = get_tool_cmdline(cnf, 'sort')
    cmdline = '{sort} -k1,1V -k2,2n -k3,3n {bed_path}'.format(**locals())
    output_fpath = intermediate_fname(work_dir, bed_path, 'sorted')
    _call(cmdline, output_fpath)
    return output_fpath


# def total_bed_length(bed_fpath):
#     cmdline = 'cat {bed} | awk -F"\t" ' \
#               '"BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}"'.format(**locals())
#     return int(_call_check_output(cmdline))


def intersect_bed(cnf, bed1, bed2, work_dir):
    bed1_fname, _ = splitext_plus(basename(bed1))
    bed2_fname, _ = splitext_plus(basename(bed2))
    output_fpath = join(work_dir, bed1_fname + '__' + bed2_fname + '.bed')
    bedtools = get_tool_cmdline(cnf, 'bedtools')
    cmdline = '{bedtools} intersect -a {bed1} -b {bed2} -u'.format(**locals())
    _call(cmdline, output_fpath)
    return output_fpath


# TODO very slow :(
def number_of_mapped_reads(cnf, bam):
    samtools = get_tool_cmdline(cnf, 'samtools')
    cmdline = '{samtools} view -c -F 4 {bam}'.format(**locals())
    res = _call_check_output(cmdline)
    return int(res)


# TODO very slow :(
def number_of_unmapped_reads(cnf, bam):
    samtools = get_tool_cmdline(cnf, 'samtools')
    cmdline = '{samtools} view -c -f 4 {bam}'.format(**locals())
    res = _call_check_output(cmdline)
    return int(res)


# TODO very slow :(
def number_of_reads(cnf, bam):
    samtools = get_tool_cmdline(cnf, 'samtools')
    cmdline = '{samtools} view -c {bam}'.format(**locals())
    res = _call_check_output(cmdline)
    return int(res)


# TODO very slow :(
def number_mapped_reads_on_target(cnf, bed, bam):
    samtools = get_tool_cmdline(cnf, 'samtools')
    cmdline = '{samtools} view -c -F 4 -L {bed} {bam}'.format(**locals())
    res = _call_check_output(cmdline)
    return int(res)


# TODO very slow :(
def number_bases_in_aligned_reads(cnf, bam):
    samtools = get_tool_cmdline(cnf, 'samtools')
    cmdline = '{samtools} depth {bam}'.format(**locals())
    proc = _call_and_open_stdout(cmdline)
    count = 0
    while True:
        coverage_line = proc.stdout.readline()
        if coverage_line:
            values = coverage_line.strip().split('\t')
            count += int(values[2])
    return count


#Region = namedtuple('Region', 'line '
#                              'size '
#                              'bases_by_depth '
#                              'avg_depth '
#                              'std_dev '
#                              'percent_within_normal')
#

def _sum_up_region(region):
    depth_sum = sum(
        depth * bases
        for depth, bases
        in region.bases_by_depth.items())
    region.avg_depth = float(depth_sum) / region.get_size()

    sum_of_sq_var = sum(
        (depth - region.avg_depth)**2 * bases
        for depth, bases
        in region.bases_by_depth.items())
    region.std_dev = math.sqrt(float(sum_of_sq_var) / region.get_size())

    bases_within_normal = sum(
        bases
        for depth, bases
        in region.bases_by_depth.items()
        if math.fabs(region.avg_depth - depth) <= 0.2 * region.avg_depth)

    region.percent_within_normal = 100.0 * bases_within_normal / region.get_size() if region.get_size else None
    return region


class Region():
    def __init__(self, sample=None, chrom=None,
                 start=None, end=None, size=None, extra_fields=list(),
                 feature=None, depth_thresholds=None):
        self.sample = sample
        self.chrom = chrom
        self.start = start
        self.end = end
        self.size = size
        self.feature = feature
        self.extra_fields = extra_fields

        self.depth_thresholds = depth_thresholds or []
        self.bases_by_depth = OrderedDict([(depth_thres, 0.0)
                                           for depth_thres in depth_thresholds])
        self.avg_depth = None
        self.std_dev = None
        self.percent_within_normal = None

    def get_size(self):
        if self.size:
            return self.size
        if self.start is None or self.end is None:
            return None
        return abs(self.end - self.start)

    def key(self):
        return hash((self.sample, self.chrom, self.start, self.end))

    def __str__(self):
        ts = [self.sample, self.chrom, self.start, self.end, self.feature] + self.extra_fields
        return '\t'.join(map(str, ts))


def get_target_depth_analytics(cnf, sample, bed, bam, feature, depth_thresholds):
    total_regions_count = 0
    total_bed_size = 0
    max_depth = 0

    bedtools = get_tool_cmdline(cnf, 'bedtools')
    cmdline = '{bedtools} coverage -abam {bam} -b {bed} -hist'.format(**locals())

    bases_per_depth_per_region = []
    cur_region = None

    for next_line in _call_and_open_stdout(cmdline).stdout:
        if not next_line.strip() or next_line.startswith('#'):
            continue

        tokens = next_line.strip().split()
        chrom = tokens[0]
        start, end = None, None
        depth, bases, region_size = map(int, tokens[-4:-1])

        region_tokens = ()

        if next_line.startswith('all'):
            max_depth = max(max_depth, depth)
            total_bed_size += bases

            region_tokens, extra_tokens = (sample, chrom, None, None), []
        else:
            start, end = map(int, tokens[1:3])
            region_tokens, extra_tokens = (sample, chrom, start, end), tokens[3:-4]

        if cur_region:
            k = cur_region.key()
            h = hash(region_tokens)
            # print (k)
            # print (h)
        if not cur_region or hash(region_tokens) != cur_region.key():
            total_regions_count += 1

            if cur_region:
                _sum_up_region(cur_region)

            cur_region = Region(sample=sample, chrom=chrom,
                start=start, end=end, size=region_size,
                feature=feature, extra_fields=extra_tokens,
                depth_thresholds=depth_thresholds)
            bases_per_depth_per_region.append(cur_region)

        for depth_thres in depth_thresholds:
            if depth >= depth_thres:
                cur_region.bases_by_depth[depth_thres] += float(bases)

        if total_regions_count > 0 and total_regions_count % 100000 == 0:
            log('processed %i regions' % total_regions_count)

    if total_regions_count % 100000 != 0:
        log('processed %i regions' % total_regions_count)

    if cur_region:  # "all" region
        _sum_up_region(cur_region)

    return bases_per_depth_per_region, max_depth, total_bed_size


# TODO how to pass the data stream to samtools vs. creating file
def get_padded_bed_file(cnf, bed, genome, padding, work_dir):
    bedtools = get_tool_cmdline(cnf, 'bedtools')
    cmdline = '{bedtools} slop -i {bed} -g {genome} -b {padding}'.format(**locals())
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