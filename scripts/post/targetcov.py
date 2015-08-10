#!/usr/bin/env python
import os
from os.path import isfile, join

import __check_python_version

import sys
import shutil
from source import BaseSample
from source.logger import critical
from source.logger import err
from source.main import read_opts_and_cnfs
from source.config import defaults
from source.prepare_args_and_cnf import check_genome_resources, check_system_resources
from source.targetcov.bam_and_bed_utils import index_bam, prepare_beds
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
            (['--exons', '--exome'], dict(
                dest='exons',
                help='a BED file with real CDS regions (default Ensembl is in system_config)')
             ),
            (['--exons-no-genes'], dict(
                dest='exons_no_genes',
                help='a BED file with real CDS regions, w/o Gene records (default Ensembl is in system_config)')
             ),
            (['--original-bed'], dict(
                dest='original_target_bed',
                help='original bed file path (just for reporting)')
             ),
            (['--original-exons'], dict(
                dest='original_exons_bed',
                help='original exons genes bed file path (just for reporting)')
             ),
            (['--reannotate'], dict(
                dest='reannotate',
                help='re-annotate BED file with gene names',
                action='store_true',
                default=False)
             ),
            (['--no-prep-bed'], dict(
                dest='no_prep_bed',
                help='do not fix input beds and exons',
                action='store_true',
                default=False)
             ),
            (['-e', '--extended'], dict(
                dest='extended',
                help='extended - flagged regions and missed variants',
                action='store_true',
                default=False)
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

    exons_bed = adjust_path(cnf.exons) if cnf.exons else adjust_path(cnf.genome.exons)
    info('Exons: ' + exons_bed)

    if cnf.genes:
        genes_fpath = adjust_path(cnf.genes)
        info('Custom genes list: ' + genes_fpath)
    else:
        genes_fpath = None

    info('Using alignement ' + cnf['bam'])

    bed_fpath = cnf.bed
    if bed_fpath:
        info('Using amplicons/capture panel ' + bed_fpath)
    elif exons_bed:
        info('WGS, taking CDS as target')

    run_one(cnf, process_one, finalize_one, multiple_samples=False, output_dir=cnf.output_dir,
        exons_bed=exons_bed, exons_no_genes_bed=cnf.exons_no_genes, genes_fpath=genes_fpath)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


class Sample(BaseSample):
    def __init__(self, name, output_dir, **kwargs):
        BaseSample.__init__(self, name, output_dir, path_base=output_dir, **kwargs)


def process_one(cnf, output_dir, exons_bed, exons_no_genes_bed, genes_fpath):
    sample = Sample(cnf.sample, output_dir, bam=cnf.bam, bed=cnf.bed)

    bam_fpath = cnf.bam
    target_bed = cnf.bed
    cnf.original_exons_bed = cnf.original_exons_bed or exons_no_genes_bed or exons_bed
    cnf.original_target_bed = cnf.original_target_bed or cnf.bed

    if not cnf.no_prep_bed:
        bam_fpath, exons_bed, exons_no_genes_bed, target_bed, _ = \
            prep_files(cnf, sample.name, sample.bam, sample.bed, exons_bed)

    avg_depth, gene_by_name, reports = make_targetseq_reports(
        cnf, output_dir, sample, bam_fpath, exons_bed, exons_no_genes_bed, target_bed, genes_fpath)

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


def prep_files(cnf, sample_name, bam_fpath, bed_fpath, exons_bed):
    bam_fpath = bam_fpath
    if not bam_fpath:
        critical(sample_name + ': BAM file is required.')
    if not isfile(bam_fpath + '.bai'):
        info('Indexing bam ' + bam_fpath)
        index_bam(cnf, bam_fpath)

    target_bed = bed_fpath
    # if not sample.bed:
    #     info(sample.name + ': BED file was not provided. Using Exons as default: ' + exons_bed)

    if not exons_bed:
        err('Error: no exons specified for the genome in system config.')

    exons_bed, exons_no_genes_bed, target_bed, seq2c_bed = prepare_beds(cnf, exons_bed, target_bed)

    return bam_fpath, exons_bed, exons_no_genes_bed, target_bed, seq2c_bed


if __name__ == '__main__':
    main(sys.argv)