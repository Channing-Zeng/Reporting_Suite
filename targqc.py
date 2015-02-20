#!/usr/bin/env python

import sub_scripts.__check_python_version  # checking for python version and adding site dirs inside

import sys
import os
from os.path import relpath, join, exists, abspath, pardir, basename
from optparse import OptionParser

from source import logger
from source.config import Config, defaults
from source.prepare_args_and_cnf import add_post_bcbio_args, check_genome_resources
from source.logger import info, err, warn, critical
from source.file_utils import verify_dir, safe_mkdir, adjust_path, verify_file, adjust_system_path, remove_quotes, \
    file_exists, isfile
from source.main import determine_cnf_files, set_up_dirs
from source.standalone_targqc.submit_jobs import run_targqc
from source.ngscat.bed_file import verify_bam, verify_bed


def main():
    info(' '.join(sys.argv))
    info()

    description = 'This script generates target QC reports for each BAM provided as an input.'
    parser = OptionParser(description=description)
    add_post_bcbio_args(parser)
    parser.add_option('--work-dir', dest='work_dir', metavar='DIR')
    parser.add_option('--log-dir', dest='log_dir')
    parser.add_option('--only-summary', dest='only_summary', action='store_true')
    parser.add_option('-o', dest='output_dir', metavar='DIR', default=join(os.getcwd(), 'targetqc'))
    parser.add_option('--reannotate', dest='reannotate', action='store_true', default=False, help='re-annotate BED file with gene names')

    (opts, args) = parser.parse_args()

    if len(args) == 0:
        critical('No BAMs provided to input.')
    bam_fpaths = [abspath(a) for a in args]
    if any(not verify_bam(fpath) for fpath in bam_fpaths):
        sys.exit(1)

    determine_cnf_files(opts)
    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)

    if not cnf.project_name:
        cnf.project_name = basename(cnf.output_dir)
    info('Project name: ' + cnf.project_name)

    set_up_dirs(cnf)
    cnf.name = 'TargQC_' + cnf.project_name

    check_genome_resources(cnf)

    bed_fpath = cnf.bed
    info('Using amplicons/capture panel ' + abspath(bed_fpath))

    exons_bed_fpath = adjust_path(cnf.exons) if cnf.exons else adjust_path(cnf.genome.exons)
    info('Exons: ' + exons_bed_fpath)

    genes_fpath = None
    if cnf.genes:
        genes_fpath = adjust_path(cnf.genes)
        info('Custom genes list: ' + genes_fpath)

    if not cnf.only_summary:
        cnf.qsub_runner = adjust_system_path(cnf.qsub_runner)
        if not cnf.qsub_runner: critical('Error: qsub-runner is not provided is sys-config.')
        if not verify_file(cnf.qsub_runner): sys.exit(1)

        if not verify_bed(bed_fpath):
            sys.exit(1)

    info('*' * 70)
    info()

    run_targqc(cnf, bam_fpaths, basename(__file__), bed_fpath, exons_bed_fpath, genes_fpath)


if __name__ == '__main__':
    main()