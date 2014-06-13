#!/usr/bin/env python

import sys
from source.logger import critical
from source.targetcov.cov import run_target_cov
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from source.runner import run_one
from source.main import read_opts_and_cnfs, check_system_resources, load_genome_resources, check_inputs
from source.utils import info


REPORT_TYPES = 'summary,genes'


def main(args):
    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--bam'], 'FILE.bam', {
                'dest': 'bam',
                'help': 'used to generate some annotations by GATK'}),

            (['--bed', '--capture', '--amplicons'], 'FILE.bed', {
                'dest': 'bed',
                'help': 'capture panel/amplicons'}),

            (['--padding'], 'N', {
                'dest': 'padding',
                'help': 'input regions will be extended by this value in both directions',
                'default': 250}),

            (['--reports'], REPORT_TYPES, {
                'dest': 'reports',
                'help': 'default: --reports ' + REPORT_TYPES,
                'default': REPORT_TYPES}),

            # (['--depth-thresholds'], 'A,B,C', {
            #     'dest': 'depth_thresholds',
            #     'help': 'A,B,C..',
            #     'default': ','.join(map(int, [5, 10, 25, 50, 100, 500, 1000, 5000,
            #                                   10000, 50000, 100000, 500000]))}),
        ],
        key_for_sample_name='bam')

    check_system_resources(
        cnf,
        required=['samtools', 'bedtools'],
        optional=[])

    load_genome_resources(
        cnf,
        required=['chr_lengths', 'genes', 'exons'],
        optional=[])

    check_inputs(cnf,
        required_keys=['bam', 'bed'],
        file_keys=['bam', 'bed'])

    if 'coverage_reports' not in cnf:
        critical('No coverage_reports section in the report, cannot make coverage reports.')

    info('Using alignement ' + cnf['bam'])
    info('Using amplicons/capture panel ' + cnf['bed'])

    run_one(cnf, process_one, finalize_one)


def process_one(cnf):
    bed_fpath = cnf['bed']
    bam_fpath = cnf['bam']

    return run_target_cov(cnf, bam_fpath, bed_fpath)


def finalize_one(cnf, summary_report_fpath, gene_report_fpath):
    if summary_report_fpath:
        info('Summary report: ' + summary_report_fpath)
    if gene_report_fpath:
        info('Exons coverage report: ' + gene_report_fpath)


if __name__ == '__main__':
    main(sys.argv)