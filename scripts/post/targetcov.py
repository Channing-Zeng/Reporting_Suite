#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc


import os
from os.path import isfile, join, basename, splitext, dirname
import sys
from traceback import format_exc
import shutil
from source import BaseSample, TargQC_Sample
from source.calling_process import call
from source.logger import critical
from source.logger import err
from source.main import read_opts_and_cnfs
from source.config import defaults
from source.prepare_args_and_cnf import check_genome_resources, check_system_resources
from source.targetcov.bam_and_bed_utils import index_bam, prepare_beds, extract_gene_names_and_filter_exons
from source.targetcov.cov import make_targetseq_reports
from source.runner import run_one
from source.targetcov.flag_regions import generate_flagged_regions_report
from source.tools_from_cnf import get_system_path
from source.utils import info
from source.file_utils import adjust_path, safe_mkdir


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
                dest='prep_bed',
                help='do not fix input beds and exons',
                action='store_false',
                default=True)
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
            # (['--pcr'], dict(
            #     dest='pcr',
            #     action='store_true',
            #     help='deduplication was not perfome, thud do not try to dedup')
            #  ),
        ],
        required_keys=['bam'],
        file_keys=['bam', 'bed'],
        key_for_sample_name='bam')

    if cnf.padding:
        cnf.coverage_reports.padding = cnf.padding

    check_system_resources(
        cnf,
        required=['samtools', 'bedtools'],
        optional=[])

    check_genome_resources(cnf)

    exons_bed = adjust_path(cnf.exons) if cnf.exons else adjust_path(cnf.genome.exons)
    if exons_bed:
        info('Exons: ' + exons_bed)
    else:
        info('no exons found')

    if cnf.genes:
        genes_fpath = adjust_path(cnf.genes)
        info('Custom genes list: ' + genes_fpath)
    else:
        genes_fpath = None

    info('Using alignment ' + cnf['bam'])

    bed_fpath = cnf.bed
    if bed_fpath:
        info('Using amplicons/capture panel ' + bed_fpath)
    elif exons_bed:
        info('WGS, taking CDS as target')

    run_one(cnf, process_one, finalize_one, output_dir=cnf.output_dir,
        exons_bed=exons_bed, exons_no_genes_bed=cnf.exons_no_genes, genes_fpath=genes_fpath)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


def picard_ins_size_hist(cnf, sample, bam_fpath, output_dir):
    picard = get_system_path(cnf, 'java', 'picard')
    if picard:
        safe_mkdir(dirname(sample.picard_ins_size_hist_txt_fpath))
        safe_mkdir(dirname(sample.picard_ins_size_hist_pdf_fpath))
        info()
        info('Picard ins size hist for "' + basename(bam_fpath) + '"')
        cmdline = '{picard} CollectInsertSizeMetrics' \
                  ' I={bam_fpath}' \
                  ' O={sample.picard_ins_size_hist_txt_fpath}' \
                  ' H={sample.picard_ins_size_hist_pdf_fpath}' \
                  ' VALIDATION_STRINGENCY=LENIENT'

        cmdline = cmdline.format(**locals())
        call(cnf, cmdline, output_fpath=sample.picard_ins_size_hist_txt_fpath,
             stdout_to_outputfile=False, exit_on_error=False)


def process_one(cnf, output_dir, exons_bed, exons_no_genes_bed, genes_fpath):
    sample = TargQC_Sample(cnf.sample, output_dir, bam=cnf.bam, bed=cnf.bed)

    bam_fpath = cnf.bam
    target_bed = cnf.bed
    cnf.original_exons_bed = cnf.original_exons_bed or exons_no_genes_bed or exons_bed
    cnf.original_target_bed = cnf.original_target_bed or cnf.bed

    bam_fpath = bam_fpath
    if not bam_fpath:
        critical(sample.name + ': BAM file is required.')
    if not isfile(bam_fpath + '.bai'):
        info('Indexing bam ' + bam_fpath)
        index_bam(cnf, bam_fpath)

    gene_names_list = None
    if cnf.prep_bed:
        info('Preparing the BED file.')
        exons_bed, exons_no_genes_bed, target_bed, seq2c_bed = prepare_beds(cnf, exons_bed, target_bed)

        gene_names_set, gene_names_list, target_bed, exons_bed, exons_no_genes_bed = \
            extract_gene_names_and_filter_exons(cnf, target_bed, exons_bed, exons_no_genes_bed, genes_fpath)
    else:
        info('The BED file is ready, skipping preparing.')
        with open(genes_fpath) as f:
            gene_names_list = [g.strip() for g in f.read().split('\n')]

    avg_depth, gene_by_name, reports = make_targetseq_reports(
        cnf, output_dir, sample, bam_fpath, exons_bed, exons_no_genes_bed, target_bed, gene_names_list)

    picard_ins_size_hist(cnf, sample, bam_fpath, output_dir)

    if cnf.extended:
        try:
            info('Generating flagged regions report...')
            flagged_report = generate_flagged_regions_report(cnf, cnf.output_dir, sample, avg_depth, gene_by_name)
        except:
            err(format_exc())

    return reports


def finalize_one(cnf, *args):
    summary_report, gene_report = args[:2]

    if summary_report.txt_fpath:
        info('Summary report: ' + summary_report.txt_fpath)

    if gene_report:
        if gene_report.txt_fpath:
            info('All regions: ' + gene_report.txt_fpath + ' (' + str(len(gene_report.rows)) + ' regions)')

    if len(args) > 2:
        # key_genes_report = args[2]
        # if key_genes_report:
        #     info('Key genes: ' + key_genes_report.tsv_fpath)
        selected_regions_report = args[2]
        if selected_regions_report.txt_fpath:
            info('Flagged regions: ' + selected_regions_report.txt_fpath +
                 ' (' + str(len(selected_regions_report.rows)) + ' regions)')


if __name__ == '__main__':
    main(sys.argv)