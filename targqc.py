#!/usr/bin/env python
import __check_python_version  # checking for python version and adding site dirs inside

import sys
import os
from os.path import relpath, join, exists, abspath, pardir, basename, splitext
from optparse import OptionParser

from source import logger
import source
from source.config import Config, defaults
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_reuse_marker_genome, check_genome_resources, determine_run_cnf, \
    determine_sys_cnf
from source.logger import info, err, warn, critical, send_email
from source.file_utils import verify_dir, safe_mkdir, adjust_path, verify_file, adjust_system_path, remove_quotes, \
    file_exists, isfile
from source.main import set_up_dirs
from source.targetcov.bam_and_bed_utils import prepare_beds, extract_gene_names_and_filter_exons
from source.targetcov.submit_jobs import run_targqc
from source.ngscat.bed_file import verify_bam, verify_bed
from source.targetcov.summarize_targetcov import get_bed_targqc_inputs


def proc_args(argv):
    info(' '.join(sys.argv))
    info()

    description = 'This script generates target QC reports for each BAM provided as an input.'
    parser = OptionParser(description=description)
    add_cnf_t_reuse_prjname_reuse_marker_genome(parser)
    parser.add_option('--log-dir', dest='log_dir')
    parser.add_option('--only-summary', dest='only_summary', action='store_true')
    parser.add_option('-o', dest='output_dir', metavar='DIR', default=join(os.getcwd(), 'targetqc'))
    parser.add_option('--reannotate', dest='reannotate', action='store_true', default=False, help='re-annotate BED file with gene names')
    parser.add_option('--dedup', dest='dedup', action='store_true', default=False, help='count duplicates in coverage metrics')
    # parser.add_option('--no-dedup', dest='dedup', action='store_false', default=False, help='not counting duplicates in coverage metrics')
    parser.add_option('-e', '--extended', dest='extended', action='store_true', default=False, help='count missed variants')
    parser.add_option('--deep-seq', dest='deep_seq', action='store_true', default=False, help='deep targeted sequencing')
    parser.add_option('--no-qualimap', dest='qualimap', action='store_false', default=True, help='do not run qualimap')
    parser.add_option('--bed', dest='bed', help='BED file to run targetSeq and Seq2C analysis on.')
    parser.add_option('--exons', '--exome', dest='exons', help='Exons BED file to make targetSeq exon/amplicon regions reports.')

    (opts, args) = parser.parse_args()

    if len(args) == 0:
        critical('No BAMs provided to input.')

    sample_names, bam_fpaths = read_samples(args)

    run_cnf = determine_run_cnf(opts, is_wgs=not opts.__dict__.get('bed'), is_targeteq=opts.deep_seq)
    cnf = Config(opts.__dict__, determine_sys_cnf(opts), run_cnf)

    if not cnf.project_name:
        cnf.project_name = basename(cnf.output_dir)
    info('Project name: ' + cnf.project_name)

    cnf.proc_name = 'TargQC'
    set_up_dirs(cnf)
    # cnf.name = 'TargQC_' + cnf.project_name

    samples = [
        source.TargQC_Sample(s_name, join(cnf.output_dir, s_name), bam=bam_fpath)
            for s_name, bam_fpath in zip(sample_names, bam_fpaths)]
    samples.sort(key=lambda _s: _s.key_to_sort())

    check_genome_resources(cnf)

    target_bed, exons_bed, genes_fpath = get_bed_targqc_inputs(cnf, cnf.bed)
    exons_no_genes_bed = None
    if not target_bed:
        info('No bed is specified, using exons instead: ' + exons_bed)

    if not cnf.only_summary:
        cnf.qsub_runner = adjust_system_path(cnf.qsub_runner)
        if not cnf.qsub_runner: critical('Error: qsub-runner is not provided is sys-config.')
        verify_file(cnf.qsub_runner, is_critical=True)

    return cnf, samples, target_bed, exons_bed, genes_fpath


def main():
    cnf, samples, target_bed, exons_bed, genes_fpath = proc_args(sys.argv)

    targqc_html_fpath = run_targqc(cnf, samples, basename(__file__), target_bed, exons_bed, genes_fpath)

    # if targqc_html_fpath:
    #     send_email('TargQC report for ' + cnf.project_name + ':\n  ' + targqc_html_fpath)


def read_samples(args):
    bam_fpaths = []
    sample_names = []
    bad_bam_fpaths = []

    for arg in args or [os.getcwd()]:
        # /ngs/oncology/Analysis/bioscience/Bio_0038_KudosCellLinesExomes/Bio_0038_150521_D00443_0159_AHK2KTADXX/bcbio,Kudos159 /ngs/oncology/Analysis/bioscience/Bio_0038_KudosCellLinesExomes/Bio_0038_150521_D00443_0160_BHKWMNADXX/bcbio,Kudos160
        bam_fpath = verify_bam(arg.split(',')[0])
        if not verify_bam(bam_fpath):
            bad_bam_fpaths.append(bam_fpath)
        bam_fpaths.append(bam_fpath)
        if len(arg.split(',')) > 1:
            sample_names.append(arg.split(',')[1])
        else:
            sample_names.append(basename(splitext(bam_fpath)[0]))
    if bad_bam_fpaths:
        critical('BAM files cannot be found, empty or not BAMs:' + ', '.join(bad_bam_fpaths))

    return sample_names, bam_fpaths


if __name__ == '__main__':
    main()