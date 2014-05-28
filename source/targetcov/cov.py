#!/usr/bin/env python

from __future__ import print_function

from collections import defaultdict
from itertools import izip, chain, repeat

import sys
from os.path import join, basename
import subprocess

#downlad hg19.genome
#https://github.com/arq5x/bedtools/tree/master/genomes

#TODO
# log file
# yaml
# take folder name as a sample name (first column on the report)
from source.targetcov.Region import Region
from source.utils import intermediate_fname, timestamp, get_tool_cmdline
from source.transaction import file_transaction
from source.bcbio_utils import splitext_plus


def log(msg=''):
    print(timestamp() + msg)


def run_header_report(cnf, result_fpath, output_dir, work_dir,
                      bed, bam, chr_len_fpath,
                      depth_thresholds, padding,
                      combined_region, max_depth, total_bed_size):
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

    bases_within_threshs, avg_depth, std_dev, percent_within_normal = combined_region.sum_up(depth_thresholds)

    v_covered_bases_in_targ = bases_within_threshs.items()[0][1]
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
    append_stat(format_integer('Reads mapped on padded target', v_reads_on_padded_targ))

    v_percent_mapped_on_padded_target = 100.0 * v_reads_on_padded_targ / v_mapped_reads if v_mapped_reads else None
    append_stat(format_decimal('Percentage of reads mapped on padded target', v_percent_mapped_on_padded_target, '%'))

    # v_aligned_read_bases = number_bases_in_aligned_reads(bam)
    # append_stat(format_integer('Total aligned bases in reads', v_aligned_read_bases))

    v_read_bases_on_targ = avg_depth * total_bed_size  # sum of all coverages
    append_stat(format_integer('Read bases mapped on target', v_read_bases_on_targ))

    log('')
    append_stat(format_decimal('Average target coverage depth', avg_depth))
    append_stat(format_decimal('Std. dev. of target coverage depth', std_dev))
    append_stat(format_integer('Maximum target coverage depth', max_depth))
    append_stat(format_decimal('Percentage of target within 20% of mean depth',
                               percent_within_normal, '%'))

    # v_percent_read_bases_on_targ = 100.0 * v_read_bases_on_targ / v_aligned_read_bases \
    #     if v_aligned_read_bases else None
    # format_decimal('Percent bases in reads on target', v_percent_read_bases_on_targ, '%'),

    # format_integer('Bases covered (at least 1x) in target', bases_per_depth[1]),

    # for depth, bases in bases_per_depth.items():
    #     append_stat(format_integer('Bases on target covered at least by ' + str(depth) +
    #                                ' read' + ('s' if depth != 1 else ''), bases))

    for depth, bases in bases_within_threshs.items():
        percent = 100.0 * bases / total_bed_size if total_bed_size else 0
        append_stat(format_decimal('Part of target covered at least by ' + str(depth) +
                                   'x', percent, '%'))

    max_len = max(len(l.rsplit(':', 1)[0]) for l in stats)
    with file_transaction(result_fpath) as tx, open(tx, 'w') as out:
        for l in stats:
            text, val = l.rsplit(':', 1)
            # spaces = ' ' * (max_len - len(text) + 1)
            spaces = '\t'
            out.write(text + spaces + val + '\n')

    log('')
    log('Result: ' + result_fpath)
    return result_fpath


def run_amplicons_cov_report(cnf, report_fpath, sample_name, depth_threshs, regions):
    for region in regions:
        region.feature = 'Amplicon'
        region.sample = sample_name

    return run_cov_report(cnf, report_fpath, depth_threshs, regions)


def run_exons_cov_report(cnf, report_fpath, sample_name, depth_threshs, regions):
    extra_fields = ['Gene']

    for region in regions:
        region.feature = 'Exon'
        region.sample = sample_name
        region.extra_fields = region.extra_fields[:2]

    exons_and_genes = add_genes_cov_analytics(regions, gene_pos=0, exon_num_pos=1)

    return run_cov_report(cnf, report_fpath, depth_threshs,
                          exons_and_genes, extra_fields=extra_fields)


def add_genes_cov_analytics(exons, gene_pos=0, exon_num_pos=1):
    exons_and_genes = []

    genes_by_name = dict()

    def update_gene(gene, exon):
        gene.end = max(gene.end, exon.end)
        gene.size += exon.get_size()
        for depth, bases in exon.bases_by_depth.items():
            gene.add_bases_for_depth(depth, bases)

    for exon in exons:
        if len(exon.extra_fields) <= gene_pos:
            sys.exit('no gene info in exons record: ' + str(exon))

        gene_name = exon.extra_fields[gene_pos]
        gene = genes_by_name.get(gene_name)
        if gene is None:
            extra_fields = ['-'] * len(exon.extra_fields)
            extra_fields[gene_pos] = gene_name
            gene = Region(
                sample=exon.sample, chrom=exon.chrom,
                start=exon.start, end=0, size=0, feature='Gene',
                extra_fields=extra_fields)
            genes_by_name[gene_name] = gene
            exons_and_genes.append(gene)
        update_gene(gene, exon)

        exon.extra_fields[0] = exon.extra_fields[exon_num_pos]
        exons_and_genes.append(exon)

    return exons_and_genes


def run_cov_report(cnf, report_fpath, depth_threshs, regions, extra_fields=list()):
    first_fields = ['SAMPLE', 'Chr', 'Start', 'End', 'Feature']
    last_fields = ['Size', 'Mean Depth', 'Standard Dev.', 'Within 20% of Mean']
    header_fields = first_fields + extra_fields + last_fields
    header_fields += ['{}x'.format(thres) for thres in depth_threshs]
    max_lengths = map(len, header_fields)

    all_values = []

    for i, region in enumerate(regions):
        bases_within_threshs, avg_depth, std_dev, percent_within_normal = region.sum_up(depth_threshs)

        line_fields = map(str, [region.sample,
                                region.chrom,
                                region.start,
                                region.end,
                                region.feature])
        line_fields += region.extra_fields[:len(extra_fields)]
        line_fields += [str(region.get_size())]
        line_fields += ['{0:.2f}'.format(avg_depth)]
        line_fields += ['{0:.2f}'.format(std_dev)]
        line_fields += ['{0:.2f}%'.format(percent_within_normal)
                        if percent_within_normal is not None else '-']

        for depth_thres, bases in bases_within_threshs.items():
            if int(region.get_size()) == 0:
                percent_str = '-'
            else:
                percent = 100.0 * bases / region.get_size()
                percent_str = '{0:.2f}%'.format(percent)
            line_fields.append(percent_str)

        all_values.append(line_fields)
        max_lengths = map(max, izip(max_lengths, chain(map(len, line_fields), repeat(0))))

    with file_transaction(report_fpath) as tx:
        with open(tx, 'w') as out, \
                open(join(cnf['work_dir'], basename(report_fpath)), 'w') as nice_out:
            out.write('\t'.join(header_fields) + '\n')

            for line_fields in all_values:
                out.write('\t'.join(line_fields) + '\n')

            for h, l in zip(header_fields, max_lengths):
                nice_out.write(h + ' ' * (l - len(h) + 2))
            nice_out.write('\n')

            for line_tokens in all_values:
                for v, l in zip(line_tokens, max_lengths):
                    nice_out.write(v + ' ' * (l - len(v) + 2))
                nice_out.write('\n')
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


# def gnu_sort(cnf, bed_path, work_dir):
#     sort = get_tool_cmdline(cnf, 'sort')
#     cmdline = '{sort} -k1,1V -k2,2n -k3,3n {bed_path}'.format(**locals())
#     output_fpath = intermediate_fname(work_dir, bed_path, 'sorted')
#     _call(cmdline, output_fpath)
#     return output_fpath


def sort_bed(cnf, bed_fpath):
    bedtools = get_tool_cmdline(cnf, 'bedtools')
    cmdline = '{bedtools} sort -i {bed_fpath}'.format(**locals())
    output_fpath = intermediate_fname(cnf['work_dir'], bed_fpath, 'sorted')
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


def bedcoverage_hist_stats(cnf, bed, bam):
    regions, max_depth, total_bed_size = [], 0, 0

    bedtools = get_tool_cmdline(cnf, 'bedtools')
    cmdline = '{bedtools} coverage -abam {bam} -b {bed} -hist'.format(**locals())

    _total_regions_count = 0

    for next_line in _call_and_open_stdout(cmdline).stdout:
        if not next_line.strip() or next_line.startswith('#'):
            continue

        line_tokens = next_line.strip().split()
        chrom = line_tokens[0]
        start, end = None, None
        depth, bases, region_size = map(int, line_tokens[-4:-1])

        if next_line.startswith('all'):
            max_depth = max(max_depth, depth)
            total_bed_size += bases
            extra_tokens = []
        else:
            start, end = map(int, line_tokens[1:3])
            extra_tokens = line_tokens[3:-4]

        line_region_key_tokens = (None, chrom, start, end)

        if regions == [] or hash(line_region_key_tokens) != regions[-1].key():
            _total_regions_count += 1

            region = Region(sample=None, chrom=chrom,
                            start=start, end=end, size=region_size,
                            extra_fields=extra_tokens)
            regions.append(region)

        regions[-1].add_bases_for_depth(depth, bases)

        if _total_regions_count > 0 and _total_regions_count % 100000 == 0:
            log('processed %i regions' % _total_regions_count)

    if _total_regions_count % 100000 != 0:
        log('processed %i regions' % _total_regions_count)

    return regions[:-1], regions[-1], max_depth, total_bed_size


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


def format_integer(name, value, unit=''):
    value = int(value)
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


