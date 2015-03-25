#!/usr/bin/env python

import __check_python_version

import sys
import shutil
from source import BaseSample
from source.main import read_opts_and_cnfs
from source.config import defaults
from source.prepare_args_and_cnf import check_genome_resources, check_system_resources
from source.targetcov.cov import make_targetseq_reports
from source.runner import run_one
from source.utils import info
from source.file_utils import adjust_path


def main(args):
    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--bam'], dict(
                dest='bam',
                help='a path to the BAM file to study')
             ),
            (['--bed', '--capture', '--amplicons'], dict(
                dest='bed',
                help='a BED file for capture panel or amplicons')
             ),
            (['--reannotate'], dict(
                dest='reannotate',
                help='re-annotate BED file with gene names',
                action='store_true',
                default=False)
             ),
            (['--count-dups'], dict(
                dest='count_dups',
                help='count duplicates when calculating coverage metrics',
                action='store_true',
                default=False)
             ),
            (['--exons', '--exome'], dict(
                dest='exons',
                help='a BED file with real CDS regions (default Ensembl is in system_config)')
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
        ],
        required_keys=['bam', 'bed'],
        file_keys=['bam', 'bed'],
        key_for_sample_name='bam')

    check_system_resources(
        cnf,
        required=['samtools', 'bedtools'],
        optional=[])

    check_genome_resources(cnf)

    exons_bed_fpath = adjust_path(cnf.exons) if cnf.exons else adjust_path(cnf.genome.exons)
    info('Exons: ' + exons_bed_fpath)

    if cnf.genes:
        genes_fpath = adjust_path(cnf.genes)
        info('Custom genes list: ' + genes_fpath)
    else:
        genes_fpath = None

    info('Using alignement ' + cnf['bam'])

    bed_fpath = cnf.bed or cnf.genome.az_exome or exons_bed_fpath
    info('Using amplicons/capture panel ' + bed_fpath)

    run_one(cnf, process_one, finalize_one, multiple_samples=False, output_dir=cnf.output_dir, exons_bed_fpath=exons_bed_fpath, genes_fpath=genes_fpath)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


class Sample(BaseSample):
    def __init__(self, name, output_dir, **kwargs):
        BaseSample.__init__(self, name, output_dir, path_base=output_dir, **kwargs)


def process_one(cnf, output_dir, exons_bed_fpath, genes_fpath):
    sample = Sample(cnf.name, output_dir, bam=cnf.bam, bed=cnf.bed)
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