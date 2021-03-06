#!/usr/bin/env python
import bcbio_postproc  # checking for python version and adding site dirs inside

import sys
import os
from os.path import relpath, join, exists, abspath, pardir, basename
from optparse import OptionParser

from source import logger
from source.config import Config, defaults
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug, check_genome_resources, determine_sys_cnf, \
    determine_run_cnf
from source.logger import info, err, warn, critical, send_email
from source.file_utils import verify_dir, safe_mkdir, adjust_path, verify_file, adjust_system_path, remove_quotes, \
    file_exists, isfile
from source.main import set_up_dirs
from source.targetcov.submit_jobs import run_targqc
from source.targetcov.bam_and_bed_utils import verify_bam, verify_bed


def main():
    info(' '.join(sys.argv))
    info()

    description = 'This script generates target QC reports for each BAM provided as an input.'
    parser = OptionParser(description=description)
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser, threads=1)
    parser.add_option('--work-dir', dest='work_dir', metavar='DIR')
    parser.add_option('--log-dir', dest='log_dir')
    parser.add_option('--only-summary', dest='only_summary', action='store_true')
    parser.add_option('-o', dest='output_dir', metavar='DIR', default=join(os.getcwd(), 'targetqc'))
    parser.add_option('--reannotate', dest='reannotate', action='store_true', default=False, help='re-annotate BED file with gene names')
    parser.add_option('--dedup', dest='dedup', action='store_true', default=False, help='count duplicates in coverage metrics')
    parser.add_option('--bed', dest='bed', help='BED file to run targetSeq and Seq2C analysis on.')
    parser.add_option('--exons', '--exome', '--features', dest='features', help='Annotated CDS/Exon/Gene/Transcripts BED file to make targetSeq exon/amplicon regions reports.')

    (opts, args) = parser.parse_args()
    logger.is_debug = opts.debug

    if len(args) == 0:
        critical('No BAMs provided to input.')
    bam_fpaths = list(set([abspath(a) for a in args]))

    bad_bam_fpaths = []
    for fpath in bam_fpaths:
        if not verify_bam(fpath):
            bad_bam_fpaths.append(fpath)
    if bad_bam_fpaths:
        critical('BAM files cannot be found, empty or not BAMs:' + ', '.join(bad_bam_fpaths))

    run_cnf = determine_run_cnf(opts, is_wgs=not opts.__dict__.get('bed'))
    cnf = Config(opts.__dict__, determine_sys_cnf(opts), run_cnf)

    if not cnf.project_name:
        cnf.project_name = basename(cnf.output_dir)
    info('Project name: ' + cnf.project_name)

    cnf.proc_name = 'TargQC'
    set_up_dirs(cnf)
    # cnf.name = 'TargQC_' + cnf.project_name

    check_genome_resources(cnf)

    verify_bed(cnf.bed, is_critical=True)
    bed_fpath = adjust_path(cnf.bed)
    info('Using amplicons/capture panel ' + bed_fpath)

    features_bed_fpath = adjust_path(cnf.features) if cnf.features else adjust_path(cnf.genome.features)
    info('Features: ' + features_bed_fpath)

    genes_fpath = None
    if cnf.genes:
        genes_fpath = adjust_path(cnf.genes)
        info('Custom genes list: ' + genes_fpath)

    if not cnf.only_summary:
        cnf.qsub_runner = adjust_system_path(cnf.qsub_runner)
        if not cnf.qsub_runner: critical('Error: qsub-runner is not provided is sys-config.')
        verify_file(cnf.qsub_runner, is_critical=True)

    info('*' * 70)
    info()

    targqc_html_fpath = run_targqc(cnf, cnf.output_dir, bam_fpaths, bed_fpath, features_bed_fpath, genes_fpath)
    if targqc_html_fpath:
        send_email(cnf, 'TargQC report for ' + cnf.project_name + ':\n  ' + targqc_html_fpath)


if __name__ == '__main__':
    main()