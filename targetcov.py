#!/usr/bin/env python

import sys
from source.targetcov.cov import run_target_cov
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from source.runner import run_one
from source.main import read_opts_and_cnfs, check_system_resources, load_genome_resources
from source.utils import info


REPORT_TYPES = 'summary,genes'


def main(args):
    required_keys = ['bam', 'bed']
    optional_keys = []

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
        required_keys=required_keys,
        optional_keys=optional_keys)

    check_system_resources(cnf, ['samtools', 'bedtools'])
    load_genome_resources(cnf, ['chr_lengths', 'genes', 'exons'])

    run_one(cnf, required_keys, optional_keys, process_one, finalize_one)


def process_one(cnf, *inputs):
    return run_target_cov(cnf, *inputs)


def finalize_one(cnf, summary_report_fpath, gene_report_fpath):
    if summary_report_fpath:
        info('Summary report: ' + summary_report_fpath)
    if gene_report_fpath:
        info('Exons coverage report: ' + gene_report_fpath)


if __name__ == '__main__':
    main(sys.argv)