#!/usr/bin/env python

import sys
from source.config import Defaults
from source.logger import critical
from source.targetcov.cov import run_target_cov
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from source.runner import run_one
from source.main import read_opts_and_cnfs, check_system_resources, load_genome_resources, check_inputs
from source.utils import info


def main(args):
    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--bam'], dict(
                dest='bam',
                help='path to the BAM file')
             ),
            (['--bed', '--capture', '--amplicons'], dict(
                dest='bed',
                help='capture panel/amplicons')
             ),
            (['--padding'], dict(
                dest='padding',
                help='integer indicating the number of bases to extend each target region up and down-stream',
                type='int',
                default=Defaults.coverage_reports['padding'])
             ),
            (['--reports'], dict(
                dest='reports',
                metavar=Defaults.coverage_reports['report_types'],
                help='Comma-separated report names',
                default=Defaults.coverage_reports['report_types'])
             ),
            (['--depth-thresholds'], dict(
                dest='depth_thresholds',
                help='A,B,C..',
                metavar='A,B,C',
                default=','.join(map(str, Defaults.coverage_reports['depth_thresholds'])))
             ),
        ],
        required_keys=['bam', 'bed'],
        file_keys=['bam', 'bed'],
        key_for_sample_name='bam')

    check_system_resources(
        cnf,
        required=['samtools', 'bedtools'],
        optional=[])

    load_genome_resources(
        cnf,
        required=['chr_lengths', 'exons'],
        optional=['genes'])

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