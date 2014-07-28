#!/usr/bin/env python
import sys
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from os.path import abspath, dirname, realpath, pardir, join
from site import addsitedir
source_dir = abspath(dirname(realpath(__file__)))
addsitedir(join(source_dir, 'ext_modules'))

import shutil
from source.main import read_opts_and_cnfs, check_system_resources, load_genome_resources
from source.config import Defaults
from source.targetcov.cov import run_target_cov
from source.runner import run_one
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
                help='integer indicating the number of bases to extend each target region up and down-stream. '
                     'Default is ' + str(Defaults.coverage_reports['padding']),
                type='int')
             ),
            (['--reports'], dict(
                dest='report_types',
                metavar=Defaults.coverage_reports['report_types'],
                help='Comma-separated report names.')
             ),
            (['--depth-thresholds'], dict(
                dest='depth_thresholds',
                metavar='A,B,C',
                help='Default: ' + ','.join(map(str,
                      Defaults.coverage_reports['depth_thresholds']))),
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
        required=['seq', 'exons'],
        optional=['chr_lengths', 'genes'])

    if cnf['report_types']:
        cnf['coverage_reports']['report_types'] = cnf['report_types']
        del cnf['report_types']

    info('Using alignement ' + cnf['bam'])
    info('Using amplicons/capture panel ' + cnf['bed'])

    run_one(cnf, process_one, finalize_one)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


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