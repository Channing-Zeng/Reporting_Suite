#!/usr/bin/env python

import __common

import sys
import shutil
from source.bcbio_structure import BCBioSample
from source.main import read_opts_and_cnfs, check_system_resources, check_genome_resources
from source.config import defaults
from source.targetcov.cov import make_targetseq_reports
from source.runner import run_one
from source.utils import info
from source.file_utils import adjust_path


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
            (['--exons', '--exome'], dict(
                dest='exons',
                help='exome (default is in system_config)')
             ),
            (['--genes'], dict(
                dest='genes',
                help='custom list of genes')
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

    if cnf.exons:
        exons_bed_fpath = adjust_path(cnf.exons)
    else:
        exons_bed_fpath = adjust_path(cnf.genome.exons)
    info('Exons: ' + exons_bed_fpath)

    if cnf.genes:
        genes_fpath = adjust_path(cnf.genes)
        info('Custom genes list: ' + exons_bed_fpath)
    else:
        genes_fpath = None

    info('Using alignement ' + cnf['bam'])
    info('Using amplicons/capture panel ' + cnf['bed'])

    run_one(cnf, process_one, finalize_one, multiple_samples=False, exons_bed_fpath=exons_bed_fpath, genes_fpath=genes_fpath)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


def process_one(cnf, exons_bed_fpath, genes_fpath):
    sample = BCBioSample(cnf.name, bam=cnf.bam, bed=cnf.bed)
    return make_targetseq_reports(cnf, sample, exons_bed_fpath, genes_fpath)  # cnf.vcfs_by_callername


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