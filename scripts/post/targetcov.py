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
from source.targetcov.flag_regions import generate_flagged_regions_report
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
            (['--dedup'], dict(
                dest='dedup',
                help='count duplicates when calculating coverage metrics',
                action='store_true',
                default=False)
             ),
            (['-e', '--extended'], dict(
                dest='extended',
                help='extended - flagged regions and missed variants',
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
        required_keys=['bam'],
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

    bed_fpath = cnf.bed
    if bed_fpath:
        info('Using amplicons/capture panel ' + bed_fpath)
    elif exons_bed_fpath:
        info('WGS, taking CDS as target')

    run_one(cnf, process_one, finalize_one, multiple_samples=False, output_dir=cnf.output_dir, exons_bed_fpath=exons_bed_fpath, genes_fpath=genes_fpath)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


class Sample(BaseSample):
    def __init__(self, name, output_dir, **kwargs):
        BaseSample.__init__(self, name, output_dir, path_base=output_dir, **kwargs)


def process_one(cnf, output_dir, exons_bed_fpath, genes_fpath):
    sample = Sample(cnf.name, output_dir, bam=cnf.bam, bed=cnf.bed)

    avg_depth, gene_by_name, reports = make_targetseq_reports(cnf, sample, exons_bed_fpath, genes_fpath)  # cnf.vcfs_by_callername

    if cnf.extended:
        info('Generating flagged regions report...')
        flagged_report = generate_flagged_regions_report(cnf, cnf.output_dir, sample,
            avg_depth, gene_by_name, cnf.coverage_reports.depth_thresholds)
        reports.append(flagged_report)

    return reports


def finalize_one(cnf, *args):
    summary_report, gene_report = args[:2]

    if summary_report.txt_fpath:
        info('Summary report: ' + summary_report.txt_fpath)

    if gene_report.txt_fpath:
        info('All regions: ' + gene_report.txt_fpath + ' (' + str(len(gene_report.regions)) + ' regions)')

    if len(args) >= 3:
        selected_regions_report = args[2]
        if selected_regions_report.txt_fpath:
            info('Selected regions: ' + selected_regions_report.txt_fpath +
                 ' (' + str(len(selected_regions_report.regions)) + ' regions)')


if __name__ == '__main__':
    main(sys.argv)