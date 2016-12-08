#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc

import math
import os
from os.path import join, basename, splitext
import sys
import shutil
from collections import defaultdict
from itertools import izip
from optparse import OptionParser

from source.bcbio.bcbio_structure import BCBioStructure, process_post_bcbio_args
from source.calling_process import call
from source.logger import info, critical, warn, err
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug, set_up_log
from source.file_utils import safe_mkdir, adjust_path, verify_file, splitext_plus, add_suffix
from source.targetcov.bam_and_bed_utils import sort_bed
from source.targetcov.flag_regions import _intersect_with_tricky_regions, tricky_regions_fnames_d
from source.targetcov.summarize_targetcov import get_val, get_float_val
from source.tools_from_cnf import get_system_path
from source.utils import is_us
from source.variants.vcf_processing import bgzip_and_tabix

from tools.prepare_data_for_exac import calculate_coverage_use_grid, get_exac_dir, add_project_to_exac

from ngs_reporting.combine_reports import get_uniq_sample_key


def main():
    info(' '.join(sys.argv))
    info()
    description = 'This script evaluate capture target.'

    parser = OptionParser(description=description)
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser)

    parser.add_option('--log-dir', dest='log_dir')
    parser.add_option('--exac-only-filtering', dest='prepare_for_exac', action='store_true', default=False, help='Export filtered regions to ExAC browser.')
    parser.add_option('--exac', dest='add_to_exac', action='store_true', default=False, help='Export coverage data to ExAC browser.')
    parser.add_option('--bed', '--capture', '--amplicons', dest='bed', help='BED file to overlap.')
    parser.add_option('--tricky-regions', dest='tricky_regions', action='store_true', default=False,
                      help='Use high GC, low GC, low complexity regions to overlap.')
    parser.add_option('--min-percent', dest='min_percent', default='0.5', help='Minimal percent of region which has low coverage.')
    parser.add_option('--min-ratio', dest='min_ratio', default='0.5', help='Minimal percent of samples which share the same feature.')
    parser.add_option('--min-depth', dest='min_depth', help='Coverage threshold.')
    parser.add_option('--metadata', dest='metadata', help='Samples type for each project '
                      '(plasma, cell_line, ffpe, deepseq, exome, wgs).')
    parser.add_option('-o', dest='output_dir', help='Output directory.')

    cnf, bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths, tags, is_wgs_in_bcbio, is_rnaseq \
        = process_post_bcbio_args(parser)

    if not cnf.project_name:
        cnf.add_to_exac = False
        cnf.project_name = 'CaptureTargetEvaluation'

    if cnf.prepare_for_exac:
        cnf.output_dir = join(get_exac_dir(cnf), 'coverage', cnf.project_name)
    elif cnf.output_dir is None:
        cnf.output_dir = join(os.getcwd(), cnf.project_name)

    cnf.output_dir = safe_mkdir(adjust_path(cnf.output_dir))
    cnf.work_dir = safe_mkdir(join(cnf.output_dir, 'work'))
    cnf.log_dir = safe_mkdir(join(cnf.work_dir), 'log')

    cnf.min_percent = 1 - float(cnf.min_percent)
    cnf.min_ratio = float(cnf.min_ratio)
    if cnf.min_depth:
        cnf.min_depth = int(cnf.min_depth)
    if cnf.metadata:
        cov_thresholds = {'deepseq': 250, 'plasma': 100, 'exome': 20, 'ffpe': 10, 'cell_line': 10, 'wgs': 10}
        cnf.min_depths = [cov_thresholds[type] for type in cnf.metadata.split(',')]

    if len(bcbio_project_dirpaths) < 1:
        critical('Usage: ' + __file__ + ' project_bcbio_path [project_bcbio_path] [-o output_dir]')

    info()
    info('*' * 70)

    safe_mkdir(cnf.output_dir)

    if cnf.log_dir:
        info('log_dirpath: ' + cnf.log_dir)
        safe_mkdir(cnf.log_dir)
        set_up_log(cnf, 'evaluate_capture_target', cnf.project_name, cnf.output_dir)

    bcbio_structures = []
    for bcbio_project_dirpath, bcbio_cnf, final_dirpath in zip(
            bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths):
        bs = BCBioStructure(cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath)
        bcbio_structures.append(bs)

    cnf.work_dir = cnf.work_dir or adjust_path(join(cnf.output_dir, 'work'))
    safe_mkdir(cnf.work_dir)

    info('')
    info('*' * 70)
    regions_fpath = evaluate_capture(cnf, bcbio_structures)
    if cnf.add_to_exac:
        if not is_us():
            err('Exposing to ExAC browser is available only on US server')
            return
        output_dirpath = join(get_exac_dir(cnf), 'coverage', cnf.project_name)
        safe_mkdir(output_dirpath)
        if regions_fpath and regions_fpath != join(output_dirpath, basename(regions_fpath)):
            shutil.copy(regions_fpath, join(output_dirpath, basename(regions_fpath)))
            shutil.copy(regions_fpath + '.tbi', join(output_dirpath, basename(regions_fpath + '.tbi')))
        samples = []
        sample_names = [s.name for bs in bcbio_structures for s in bs.samples]
        for bs in bcbio_structures:
            for sample in bs.samples:
                sample.name = get_uniq_sample_key(bs.project_name, sample, sample_names)
                samples.append(sample)
        calculate_coverage_use_grid(cnf, samples, output_dirpath)
        add_project_to_exac(cnf)
    else:
        info('Use --exac if you want to export data to ExAC browser')
    info('Done.')


def evaluate_capture(cnf, bcbio_structures):
    samples = [s for bs in bcbio_structures for s in bs.samples]
    min_samples = math.ceil(cnf.min_ratio * len(samples))

    info('Filtering regions by depth')
    regions = check_regions_depth(cnf, bcbio_structures, min_samples)
    if not regions:
        err('No regions were filtered.')
        return None
    if cnf.bed or cnf.tricky_regions:
        regions = intersect_regions(cnf, bcbio_structures, regions, min_samples)

    regions_fname = 'filtered_regions.txt'
    regions_fpath = join(cnf.output_dir, add_suffix(regions_fname, str(cnf.min_depth)) if cnf.min_depth else regions_fname)
    with open(regions_fpath, 'w') as out:
        out.write('## Minimal percent of region with low coverage: ' + str((1 - cnf.min_percent) * 100) + '%\n')
        out.write('## Minimal percent of samples that share the same feature: ' + str(cnf.min_ratio * 100) + '%\n')
        if not cnf.min_depth:
            out.write('## Coverage threshold Nx is 10x for cell line and 100x for plasma\n')
        else:
            out.write('## Coverage threshold Nx is ' + str(cnf.min_depth) + 'x\n')
        out.write('\t'.join(['#Chr', 'Start', 'End', 'Size', 'Gene', 'Depth<Nx', 'SamplesSharingSameFeature', 'Annotation']) + '\n')
        for region in sorted(regions, key=lambda x: (x[0], int(x[1]))):
            out.write('\t'.join([str(val) for val in region]) + '\n')

    info()
    info(str(len(regions)) + ' regions were saved into ' + regions_fpath)
    bgzip_and_tabix(cnf, regions_fpath, tabix_parameters='-p bed')
    return regions_fpath


def intersect_regions(cnf, bcbio_structures, all_regions, min_samples):
    all_regions_fname = 'all_regions.bed'
    all_regions_bed_fpath = join(cnf.output_dir, add_suffix(all_regions_fname, str(cnf.min_depth)) if cnf.min_depth else all_regions_fname)

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
                filename = tricky_regions_fnames_d[basename(fs[7]).split('.')[0]]
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
        else:
            r = list(r)
            r.append('')
        regions.append(r)
    os.remove(intersection_fpath)
    return regions


def check_regions_depth(cnf, bcbio_structures, min_samples):
    regions = defaultdict(list)
    filtered_regions = []
    total_samples_count = len([s for bs in bcbio_structures for s in bs.samples])
    for i, bs in enumerate(bcbio_structures):
        depth_threshold = cnf.min_depth or cnf.min_depths[i]
        for s in bs.samples:
            tsv_fpath = s.targetcov_detailed_tsv
            if not verify_file(tsv_fpath, is_critical=False):
                continue

            chrom_col = None
            start_col = None
            end_col = None
            size_col = None
            gene_col = None
            strand_col = None
            feature_col = None
            biotype_col = None
            tx_col = None
            ave_depth_col = None
            med_depth_col = None
            stddev_depth_col = None
            wn20pcnt_col = None
            depth_thresholds = None
            with open(tsv_fpath) as f_inp:
                for l in f_inp:
                    if l.startswith('#'):
                        def filter_digits(s):
                            return ''.join(c for c in s if c.isdigit())
                        fs = l.strip('\n').split('\t')
                        chrom_col = fs.index('#Chr')
                        start_col = fs.index('Start')
                        end_col = fs.index('End')
                        gene_col = fs.index('Gene')
                        strand_col = fs.index('Strand')
                        feature_col = fs.index('Feature')
                        biotype_col = fs.index('Biotype')
                        tx_col = fs.index('Transcript')
                        ave_depth_col = fs.index('Avg depth')
                        med_depth_col = fs.index('Median depth')
                        stddev_depth_col = fs.index('Std dev')
                        wn20pcnt_col = fs.index('W/n 20% of median')
                        depth_thresholds = [int(filter_digits(d)) for d in fs if d.endswith('x') and filter_digits(d)]
                        if depth_threshold not in depth_thresholds:
                            err('Depth ' + str(depth_threshold) + 'x is not used in ' + tsv_fpath)
                            break
                        continue
                    fs = l.strip('\n').split('\t')  # Chr	Start	End	Size	Gene	Strand	Feature	Biotype	TranscriptID    Min depth	Ave depth	Std dev	W/n 20% of ave depth	1x	5x	10x	25x	50x	100x	500x	1000x	5000x	10000x	50000x
                    chrom = fs[chrom_col]
                    start = fs[start_col]
                    end = fs[end_col]
                    size = fs[size_col]
                    symbol = fs[gene_col]
                    strand = fs[strand_col]
                    feature = fs[feature_col]
                    biotype = fs[biotype_col]
                    transcript_id = fs[tx_col]
                    med_depth = fs[med_depth_col]
                    stddev_depth = fs[stddev_depth_col]
                    wn20pcnt = fs[wn20pcnt_col]
                    pcnt_val_by_thresh = fs[wn20pcnt_col + 1:]

                    symbol = get_val(symbol)
                    chrom = get_val(chrom)
                    if start == '.' or end == '.': continue
                    if not depth_thresholds:
                        critical('No depth_thresholds header in ' + tsv_fpath)
                    cov_by_threshs = dict((t, get_float_val(f)) for t, f in izip(depth_thresholds, pcnt_val_by_thresh))

                    if feature in ['Capture']:
                        region = (chrom, start, end, size, symbol)
                        if not depth_threshold:  # no filtering
                            regions[region].append(1)
                        elif cov_by_threshs[depth_threshold] < cnf.min_percent:
                            regions[region].append(1 - cov_by_threshs[depth_threshold])

    for r, depths in regions.iteritems():
        num_samples = len(depths)
        if num_samples >= min_samples:
            percent_samples = int(num_samples * 100.0 / total_samples_count)
            str_num_samples = '{num_samples} ({percent_samples}%)'.format(**locals())
            r = list(r)
            pct_depth = str(int(sum(depths) * 100.0 / num_samples)) + '%'
            r.append(pct_depth)
            r.append(str_num_samples)
            filtered_regions.append(tuple(r))
    return filtered_regions


if __name__ == '__main__':
    main()
