#!/usr/bin/env python
# noinspection PyUnresolvedReferences

import bcbio_postproc


import os
from os.path import isfile, join, basename, splitext, dirname
import sys
from traceback import format_exc
import shutil
from optparse import SUPPRESS_HELP

from source import BaseSample, TargQC_Sample
from source.calling_process import call
from source.fastqc.fastq_utils import downsample
from source.logger import critical
from source.logger import err
from source.main import read_opts_and_cnfs
from source.config import defaults
from source.prepare_args_and_cnf import check_genome_resources, check_system_resources
from source.targetcov.bam_and_bed_utils import index_bam, prepare_beds, extract_gene_names_and_filter_exons, verify_bam
from source.targetcov.cov import make_targqc_reports
from source.runner import run_one
from source.targetcov.flag_regions import generate_flagged_regions_report
from source.tools_from_cnf import get_system_path
from source.utils import info
from source.file_utils import adjust_path, safe_mkdir, verify_file, remove_quotes


def main(args):
    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--bam'], dict(
                dest='bam',
                help='a path to the BAM file to study')
             ),
            (['-1'], dict(
                dest='l_fpath')
             ),
            (['-2'], dict(
                dest='r_fpath')
             ),
            (['--bed', '--capture', '--amplicons'], dict(
                dest='bed',
                help='a BED file for capture panel or amplicons')
             ),
            (['--exons', '--exome', '--features'], dict(
                dest='features',
                help='a BED file with real CDS/Exon/Gene/Transcript regions with annotations (default "features" is in system_config)')
             ),
            (['--exons-no-genes', '--features-no-genes'], dict(
                dest='features_no_genes',
                help='a BED file with real CDS/Exon regions with annotations, w/o Gene/Transcript records (default "features" is in system_config)')
             ),
            (['--original-bed'], dict(
                dest='original_target_bed',
                help=SUPPRESS_HELP)
             ),
            (['--original-exons', '--original-features'], dict(
                dest='original_features_bed',
                help='original features genes bed file path (just for reporting)')
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
            (['--no-dedup'], dict(
                dest='no_dedup',
                action='store_true',
                help=SUPPRESS_HELP)
             ),
            (['--downsample-to'], dict(
                dest='downsample_to',
                type='int',
                help=SUPPRESS_HELP)
             ),
        ],
        file_keys=['bam', 'l_fpath', 'r_fpath', 'bed'],
        key_for_sample_name='bam')

    if cnf.padding:
        cnf.coverage_reports.padding = cnf.padding

    check_system_resources(
        cnf,
        required=['bedtools'],
        optional=[])

    check_genome_resources(cnf)

    features_bed = adjust_path(cnf.features) if cnf.features else adjust_path(cnf.genome.features)
    if features_bed:
        info('Features: ' + features_bed)
        features_bed = verify_file(features_bed)
    else:
        info('No features BED found')

    if cnf.bed:
        cnf.bed = verify_file(cnf.bed, is_critical=True)
        info('Using amplicons/capture panel ' + cnf.bed)
    elif features_bed:
        info('WGS, taking CDS as target')

    cnf.bam = verify_bam(cnf.bam, is_critical=True)

    reports = process_one(cnf, cnf.output_dir, cnf.bam,
        features_bed=features_bed, features_no_genes_bed=cnf.features_no_genes)
    summary_report, gene_report = reports[:2]

    info('')
    info('*' * 70)
    if summary_report.txt_fpath:
        info('Summary report: ' + summary_report.txt_fpath)
    if gene_report:
        if gene_report.txt_fpath:
            info('All regions: ' + gene_report.txt_fpath + ' (' + str(len(gene_report.rows)) + ' regions)')

    if len(reports) > 2:
        selected_regions_report = reports[2]
        if selected_regions_report.txt_fpath:
            info('Flagged regions: ' + selected_regions_report.txt_fpath +
                 ' (' + str(len(selected_regions_report.rows)) + ' regions)')

    for fpaths in reports:
        if fpaths:
            ok = True
            info('Checking expected results...')
            if not isinstance(fpaths, list):
                fpaths = [fpaths]
            for fpath in fpaths:
                if isinstance(fpath, basestring):
                    if not verify_file(fpath):
                        ok = False
            if ok:
                info('The results are good.')

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


def proc_fastq(cnf, sample, l_fpath, r_fpath):
    if cnf.downsample_to:
        info('Downsampling the reads to ' + str(cnf.downsample_to))
        l_fpath, r_fpath = downsample(cnf, sample.nname, l_fpath, r_fpath, cnf.downsample_to, output_dir=cnf.work_dir, suffix='subset')

    sambamba = get_system_path(cnf, 'sambamba')
    bwa = get_system_path(cnf, 'bwa')
    seqtk = get_system_path(cnf, 'seqtk')
    bammarkduplicates = get_system_path(cnf, 'bammarkduplicates')
    if not (sambamba and bwa and seqtk and bammarkduplicates):
        critical('sambamba, BWA, seqtk and bammarkduplicates are required to align BAM')
    info()
    info('Alignming reads to the reference')
    bam_fpath = align(cnf, sample, l_fpath, r_fpath,
        sambamba, bwa, seqtk, bammarkduplicates,
        cnf.genome.bwa, cnf.is_pcr)
    bam_fpath = verify_bam(bam_fpath)
    if not bam_fpath:
        critical('Sample ' + sample + ' was not aligned successfully.')
    return bam_fpath


def picard_ins_size_hist(cnf, sample, bam_fpath, output_dir):
    picard = get_system_path(cnf, 'java', 'picard')
    if picard:
        safe_mkdir(dirname(sample.picard_ins_size_hist_txt_fpath))
        safe_mkdir(dirname(sample.picard_ins_size_hist_pdf_fpath))
        info('Picard ins size hist for "' + basename(bam_fpath) + '"')
        cmdline = '{picard} CollectInsertSizeMetrics' \
                  ' I={bam_fpath}' \
                  ' O={sample.picard_ins_size_hist_txt_fpath}' \
                  ' H={sample.picard_ins_size_hist_pdf_fpath}' \
                  ' VALIDATION_STRINGENCY=LENIENT'

        cmdline = cmdline.format(**locals())
        call(cnf, cmdline, output_fpath=sample.picard_ins_size_hist_txt_fpath,
             stdout_to_outputfile=False, exit_on_error=False)


def process_one(cnf, output_dir, bam_fpath, features_bed, features_no_genes_bed):
    sample = TargQC_Sample(cnf.sample, output_dir, bed=cnf.bed, bam=cnf.bam)
    sample.l_fpath = cnf.l_fpath
    sample.r_fpath = cnf.r_fpath

    # if not sample.bam and sample.l_fpath and sample.r_fpath:
    #     sample.bam = proc_fastq(cnf, sample, verify_file(cnf.l_fpath), verify_file(cnf.r_fpath))

    info('Using alignment ' + sample.bam)

    if not bam_fpath:
        critical(sample.name + ': BAM file is required.')

    target_bed = verify_file(cnf.bed, is_critical=True) if cnf.bed else None
    bam_fpath = verify_file(sample.bam, is_critical=True)
    index_bam(cnf, bam_fpath)

    gene_keys_list = None
    if cnf.prep_bed is not False:
        info('Preparing the BED file.')
        features_bed, features_no_genes_bed, target_bed, seq2c_bed = prepare_beds(cnf, features_bed, target_bed)

        gene_keys_set, gene_keys_list, target_bed, features_bed, features_no_genes_bed = \
            extract_gene_names_and_filter_exons(cnf, target_bed, features_bed, features_no_genes_bed)
    else:
        info('The BED file is ready, skipping preparing.')
        gene_keys_set, gene_keys_list, _, _, _ = \
            extract_gene_names_and_filter_exons(cnf, target_bed, features_bed, features_no_genes_bed)

    picard_ins_size_hist(cnf, sample, bam_fpath, output_dir)

    avg_depth, gene_by_name_and_chrom, reports = make_targqc_reports(
        cnf, output_dir, sample, bam_fpath, features_bed, features_no_genes_bed, target_bed, gene_keys_list)

    # #if cnf.extended:
    # try:
    #     info('Generating flagged regions report...')
    #     flagged_report = generate_flagged_regions_report(cnf, cnf.output_dir, sample, avg_depth, gene_by_name_and_chrom)
    #     if not flagged_report:
    #         err('Flagged regions report was not generated')
    #         err()
    # except:
    #     err(format_exc())

    return reports


if __name__ == '__main__':
    main(sys.argv)