#!/usr/bin/env python

import __common

import sys
import shutil
from source.bcbio_structure import Sample
from source.main import read_opts_and_cnfs, check_system_resources, check_genome_resources
from source.config import defaults
from source.targetcov.cov import make_targetseq_reports
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
                     'Default is ' + str(defaults['coverage_reports']['padding']),
                type='int')
             ),
            (['--depth-thresholds'], dict(
                dest='depth_thresholds',
                metavar='A,B,C',
                help='Default: ' + ','.join(map(str,
                      defaults['coverage_reports']['depth_thresholds']))),
             ),
        ],
        required_keys=['bam', 'bed'],
        file_keys=['bam', 'bed'],
        key_for_sample_name='bam')

    check_system_resources(
        cnf,
        required=['samtools', 'bedtools'],
        optional=[])

    check_genome_resources(cnf)

    info('Using alignement ' + cnf['bam'])
    info('Using amplicons/capture panel ' + cnf['bed'])

    run_one(cnf, process_one, finalize_one)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


def process_one(cnf):
    sample = Sample(cnf.name, bam=cnf.bam, bed=cnf.bed)
    return make_targetseq_reports(cnf, sample)  # cnf.vcfs_by_callername


def finalize_one(cnf, summary_report_txt_path, gene_report_fpath):
    msg = ['TargetSeq reprots finished for ' + cnf.name + ':']

    if summary_report_txt_path:
        msg.append('Summary TXT:  ' + summary_report_txt_path)
        # msg.append('Summary HTML: ' + summary_report_html_path)
        info('Summary report: ' + summary_report_txt_path)
    if gene_report_fpath:
        msg.append('Per-region report: ' + gene_report_fpath)
        info('Per-region report:')
        info('  ' + gene_report_fpath)
    # if abnormal_regions_reports:
    #     msg.append('Abnormal region reports: ')
    #     info('Abnormal region reports:')
    #     for rep in abnormal_regions_reports:
    #         msg.append('  ' + rep)
    #         info('  ' + rep)

    # send_email('\n'.join(msg))


if __name__ == '__main__':
    main(sys.argv)