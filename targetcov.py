#!/usr/bin/env python

from __future__ import print_function

from genericpath import isdir

import sys
import os
from os.path import join, expanduser
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
from source.main import common_main
from source.my_utils import verify_file, critical, step_greetings
from source.targetcov import run_cov_report, run_header_report, get_target_depth_analytics_fast, intersect_bed, log


if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    exit('Python 2, versions 2.7 and higher is supported (you are running ' +
         '.'.join(map(str, sys.version_info[:3])) + ')\n')


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

            (['--only-summary'], '', {
                'dest': 'only_summary',
                'help': '',
                'action': 'store_true',
                'default': False}),

            (['--only-regions'], '', {
                'dest': 'only_regions',
                'help': '',
                'action': 'store_true',
                'default': False}),
        ])

    genes_bed = options.get('genes') or cnf.get('genes') or cnf['genome'].get('genes')
    exons_bed = options.get('exons') or cnf.get('exons') or expanduser(cnf['genome'].get('exons'))
    chr_len_fpath = cnf.get('chr_lengths') or cnf['genome'].get('chr_lengths')
    capture_bed = options.get('capture') or cnf.get('capture')
    bam = options.get('bam') or cnf.get('bam')

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

    genes_bed = expanduser(genes_bed)
    exons_bed = expanduser(exons_bed)
    chr_len_fpath = expanduser(chr_len_fpath)
    bam = expanduser(bam)
    capture_bed = expanduser(capture_bed)

    if not verify_file(genes_bed): exit(1)
    if not verify_file(exons_bed): exit(1)
    if not verify_file(chr_len_fpath): exit(1)
    if not verify_file(bam): exit(1)
    if not verify_file(capture_bed): exit(1)

    depth_thresholds = cnf['depth_thresholds']
    padding = options.get('padding', cnf.get('padding', 250))
    output_dir = options.get('output_dir') or cnf.get('output_dir') or os.getcwd()
    print('writing to output dir ' + output_dir)
    output_dir = expanduser(output_dir)

    work_dir = join(output_dir, 'work')
    if isdir(work_dir):
        rmtree(work_dir)
    os.makedirs(work_dir)

    print('')

    #########################################
    log('Calculation of coverage statistics for the regions in the input BED file...')
    bases_per_depth_per_region, max_depth, total_bed_size = \
        get_target_depth_analytics_fast(capture_bed, bam, depth_thresholds)

    _all_reg_line, (bases_per_depth_all, all_avg_depth, all_std_dev, bases_within_normal_deviation) = \
        bases_per_depth_per_region.items()[-1]
    print(_all_reg_line)

    header_report_fpath = None
    if not options.get('only_regions'):
        header_report_fpath = run_header_report(output_dir, work_dir, capture_bed, bam, chr_len_fpath,
            depth_thresholds, padding,
            bases_per_depth_all, all_avg_depth, all_std_dev, bases_within_normal_deviation,
            max_depth, total_bed_size)

    cov_report_fpath = None
    if not options.get('only_summary'):
        step_greetings('Coverage report for regions')

        bed = capture_bed
        if genes_bed:
            log('Getting the gene regions that intersect with our capture panel.')
            bed = intersect_bed(genes_bed, capture_bed, work_dir)
        if exons_bed:
            log('Getting the exons of the genes.')
            bed = intersect_bed(exons_bed, genes_bed, work_dir)

        log('Calculation of coverage statistics for exons of the genes ovelapping with the input regions...')
        bases_per_depth_per_region, max_depth, _ = \
            get_target_depth_analytics_fast(bed, bam, depth_thresholds)

        cov_report_fpath = run_cov_report(output_dir, work_dir, bed, bam, depth_thresholds,
            bases_per_depth_per_region)

    print('')
    print('*' * 70)
    if header_report_fpath:
        log('Summary report: ' + header_report_fpath)
    if cov_report_fpath:
        log('Region coverage report: ' + cov_report_fpath)


if __name__ == '__main__':
    main(sys.argv)