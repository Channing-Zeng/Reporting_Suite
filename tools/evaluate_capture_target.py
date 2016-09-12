#!/usr/bin/env python
# noinspection PyUnresolvedReferences
from itertools import izip

import math

import bcbio_postproc

from collections import defaultdict
import os
import sys
from os.path import join, basename, splitext
from optparse import OptionParser

from source.bcbio.bcbio_structure import BCBioStructure, process_post_bcbio_args
from source.calling_process import call
from source.logger import info, critical, warn
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug, set_up_log
from source.file_utils import safe_mkdir, adjust_path, verify_file, splitext_plus
from source.targetcov.bam_and_bed_utils import sort_bed
from source.targetcov.flag_regions import _intersect_with_tricky_regions, tricky_regions_fnames_d
from source.targetcov.summarize_targetcov import get_val, get_float_val
from source.tools_from_cnf import get_system_path
from tools.prepare_data_for_exac import calculate_coverage_use_grid, get_exac_dir, add_project_to_exac


def main():
    info(' '.join(sys.argv))
    info()
    description = 'This script evaluate capture target.'

    parser = OptionParser(description=description)
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser)

    parser.add_option('--log-dir', dest='log_dir', default='-')
    parser.add_option('--bed', '--capture', '--amplicons', dest='bed', help='BED file to overlap.')
    parser.add_option('--tricky-regions', dest='tricky_regions', action='store_true', default=False,
                      help='Use high GC, low GC, low complexity regions to overlap.')
    parser.add_option('--min-percent', dest='min_percent', help='Minimal percent of region which has low coverage.')
    parser.add_option('--min-ratio', dest='min_ratio', help='Minimal percent of samples which share the same feature.')
    parser.add_option('--min-depth', dest='min_depth', help='Coverage threshold.')
    parser.add_option('--metadata', dest='metadata', help='Samples type for each project (plasma, cell_line, ffpe).')
    parser.add_option('-o', dest='output_dir', help='Output directory.')

    cnf, bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths, tags, is_wgs_in_bcbio, is_rnaseq \
        = process_post_bcbio_args(parser)

    cnf.min_percent = 1 - float(cnf.min_percent)
    cnf.min_ratio = float(cnf.min_ratio)
    if cnf.min_depth:
        cnf.min_depth = int(cnf.min_depth)
    if cnf.metadata:
        cov_thresholds = {'plasma': 100, 'ffpe': 10, 'cell_line': 10}
        cnf.min_depths = [cov_thresholds[type] for type in cnf.metadata.split(',')]

    if len(bcbio_project_dirpaths) < 1:
        critical('Usage: ' + __file__ + ' project_bcbio_path [project_bcbio_path] [-o output_dir]')

    info()
    info('*' * 70)
    add_to_exac = True

    if not cnf.project_name:
        add_to_exac = False
        cnf.project_name = 'CaptureTargetEvaluation'
        warn('Please specify --project-name if you want to export data to ExAC browser')

    if cnf.output_dir is None:
        cnf.output_dir = join(os.getcwd(), cnf.project_name)
    safe_mkdir(cnf.output_dir)

    cnf.log_dir = join(cnf.output_dir, 'log')
    info('log_dirpath: ' + cnf.log_dir)
    safe_mkdir(cnf.log_dir)
    set_up_log(cnf, 'evaluate_capture_target', cnf.project_name, cnf.output_dir)

    bcbio_structures = []
    for bcbio_project_dirpath, bcbio_cnf, final_dirpath in zip(
            bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths):
        bs = BCBioStructure(cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath)
        bcbio_structures.append(bs)

    cnf.work_dir = adjust_path(join(cnf.output_dir, 'work'))
    safe_mkdir(cnf.work_dir)

    info('')
    info('*' * 70)
    evaluate_capture(cnf, bcbio_structures)
    if add_to_exac:
        add_project_to_exac(cnf)


def evaluate_capture(cnf, bcbio_structures):
    samples = [s for bs in bcbio_structures for s in bs.samples]
    min_samples = math.ceil(cnf.min_ratio * len(samples))

    regions = check_regions_depth(cnf, bcbio_structures, min_samples)
    if cnf.bed or cnf.tricky_regions:
        regions = intersect_regions(cnf, bcbio_structures, regions, min_samples)

    regions_fpath = join(cnf.output_dir, 'filtered_regions.txt')
    with open(regions_fpath, 'w') as out:
        if not cnf.min_depth:
            out.write('## Coverage threshold Nx is 10x for cell line and 100x for plasma\n')
        else:
            out.write('## Coverage threshold Nx is ' + str(cnf.min_depth) + 'x\n')
        out.write('\t'.join(['#Chr', 'Start', 'End', 'Size', 'Gene', 'Depth<Nx', 'SamplesSharingSameFeature', 'Annotation']) + '\n')
        for region in sorted(regions, key=lambda x: (x[0], int(x[1]))):
            out.write('\t'.join([val for val in region]) + '\n')

    info()
    info(str(len(regions)) + ' regions were saved into ' + regions_fpath)


def intersect_regions(cnf, bcbio_structures, all_regions, min_samples):
    all_regions_bed_fpath = join(cnf.output_dir, 'all_regions.bed')

    with open(all_regions_bed_fpath, 'w') as out:
        if not cnf.min_depth:
            out.write('## Coverage threshold Nx is 10x for cell line and 100x for plasma\n')
        else:
            out.write('## Coverage threshold Nx is ' + str(cnf.min_depth) + 'x\n')
        out.write('\t'.join(['#Chr', 'Start', 'End', 'Size', 'Gene', 'Depth<Nx', 'SamplesSharingSameFeature']) + '\n')
        for region in all_regions:
            out.write('\t'.join([str(val) for val in region]) + '\n')

    regions_overlaps = defaultdict(lambda: defaultdict(list))
    regions = []
    if cnf.tricky_regions:
        intersection_fpath = _intersect_with_tricky_regions(cnf, all_regions_bed_fpath, 'samples')
    else:
        bed_fpath = cnf.bed
        intersection_fpath = join(cnf.work_dir, splitext(basename(all_regions_bed_fpath))[0] + '_bed.intersect')
        bedtools = get_system_path(cnf, 'bedtools')
        if not cnf.reuse_intermediate or not verify_file(intersection_fpath, silent=True, is_critical=False):
            cmdline = '{bedtools} intersect -header -a {all_regions_bed_fpath} -b {bed_fpath} -wo'.format(**locals())
            res = call(cnf, cmdline, output_fpath=intersection_fpath, max_number_of_tries=1, exit_on_error=False)
            if not res:
                return None

    output_dirpath = join(get_exac_dir(cnf), 'coverage', cnf.project_name)
    samples = [s for bs in bcbio_structures for s in bs.samples]
    calculate_coverage_use_grid(cnf, samples, output_dirpath)
    with open(intersection_fpath) as f:
        for l in f:
            l = l.strip()
            if not l or l.startswith('#'):
                continue
            fs = l.split('\t')
            chrom, start, end, size, symbol, pct_depth, num_samples = fs[:7]
            overlap_bps = int(fs[-1])
            r = (chrom, start, end, size, symbol, pct_depth, num_samples)
            if cnf.tricky_regions:
                filename = tricky_regions_fnames_d[splitext_plus(basename(fs[7]))[0]]
                regions_overlaps[r][filename].append(overlap_bps)
            else:
                regions_overlaps[r][basename(cnf.bed)].append(overlap_bps)
    for r in all_regions:
        if r in regions_overlaps:
            overlaps = ''
            chrom, start, end, size, symbol, pct_depth, num_samples = r
            overlaps_txt = ', '.join(
                fname + ': %.0f' % (sum(regions_overlaps[r][fname]) / float(size) * 100) + '%'
                for fname in regions_overlaps[r])
            r = list(r)
            r.append(overlaps_txt)
        regions.append(r)
    os.remove(intersection_fpath)
    return regions


def check_regions_depth(cnf, bcbio_structures, min_samples):
    regions = defaultdict(list)
    samples = [s for bs in bcbio_structures for s in bs.samples]
    for i, bs in enumerate(bcbio_structures):
        for s in samples:
            tsv_fpath = s.targetcov_detailed_tsv
            if not verify_file(tsv_fpath, is_critical=False):
                continue
            with open(tsv_fpath) as f_inp:
                for l in f_inp:
                    if l.startswith('#'):
                        def filter_digits(s):
                            return ''.join(c for c in s if c.isdigit())
                        fs = l.split('\t')
                        if len(fs) > 13:
                            depth_thresholds = [int(filter_digits(d)) for d in fs[13:]]
                        continue

                    fs = l.split('\t')  # Chr	Start	End	Size	Gene	Strand	Feature	Biotype	TranscriptID    Min depth	Ave depth	Std dev	W/n 20% of ave depth	1x	5x	10x	25x	50x	100x	500x	1000x	5000x	10000x	50000x
                    chrom, start, end, size, symbol, strand, feature, biotype, transcript_id, min_depth, ave_depth, std_dev, wn20pcnt = fs[:13]
                    pcnt_val_by_thresh = fs[13:]

                    symbol = get_val(symbol)
                    chrom = get_val(chrom)

                    if start == '.' or end == '.':
                        continue

                    cov_by_threshs = dict((t, get_float_val(f)) for t, f in izip(depth_thresholds, pcnt_val_by_thresh))

                    if feature in ['Capture']:
                        region = (chrom, start, end, size, symbol)
                        if not cnf.min_depth and not cnf.min_depths[i]:  # no filtering
                            regions[region] += 1
                            continue
                        min_depth = cnf.min_depth or cnf.min_depths[i]
                        if min_depth not in cov_by_threshs:
                            warn()
                            continue
                        if cov_by_threshs[min_depth] < cnf.min_percent:
                            regions[region].append(1 - cov_by_threshs[min_depth])
        filtered_regions = []

        for r, depths in regions.iteritems():
            num_samples = len(depths)
            if num_samples >= min_samples:
                percent_samples = int(num_samples * 100.0 / len(samples))
                str_num_samples = '{num_samples} ({percent_samples}%)'.format(**locals())
                r = list(r)
                pct_depth = str(int(sum(depths) * 100.0 / num_samples)) + '%'
                r.append(pct_depth)
                r.append(str_num_samples)
                filtered_regions.append(tuple(r))
        return filtered_regions


if __name__ == '__main__':
    main()
