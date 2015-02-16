#!/usr/bin/env python

import sub_scripts.__check_python_version  # checking for python version and adding site dirs inside

import sys
import os
from os.path import relpath, join, exists, abspath, pardir, basename
from optparse import OptionParser

from source.config import Config, defaults
from source.prepare_args_and_cnf import add_post_bcbio_args
from source.logger import info, err, warn, critical
from source.file_utils import verify_dir, safe_mkdir, adjust_path, verify_file, adjust_system_path
from source.main import check_genome_resources, determine_cnf_files
from source.standalone_targqc.submit_jobs import run
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
    parser.add_option('-o', dest='output_dir', metavar='DIR')
    parser.add_option('--reannotate', dest='reannotate', action='store_true', default=False, help='re-annotate BED file with gene names')

    (opts, args) = parser.parse_args()

    if len(args) == 0:
        critical('No BAMs provided to input.')
    bam_fpaths = [abspath(a) for a in args]
    if any(not verify_bam(fpath) for fpath in bam_fpaths):
        sys.exit(1)

    determine_cnf_files(opts)
    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)

    output_dir = adjust_path(cnf.output_dir or join(os.getcwd(), 'targetqc'))
    if not verify_dir(join(output_dir, pardir)): sys.exit(1)
    safe_mkdir(output_dir)
    info('Output to ' + output_dir)
    cnf.output_dir = output_dir

    if not cnf.work_dir: cnf.work_dir = join(cnf.output_dir, 'work')
    safe_mkdir(cnf.work_dir)
    cnf.log_dir = join(cnf.work_dir, 'log')
    safe_mkdir(cnf.log_dir)

    if not cnf.project_name: cnf.project_name = basename(cnf.output_dir)

    check_genome_resources(cnf)

    bed_fpath = cnf.bed

    if not cnf.only_summary:
        cnf.qsub_runner = adjust_system_path(cnf.qsub_runner)
        if not cnf.qsub_runner: critical('Error: qsub-runner is not provided is sys-config.')
        if not verify_file(cnf.qsub_runner): sys.exit(1)

        if not verify_bed(bed_fpath):
            sys.exit(1)

    info('*' * 70)
    info()

    run(cnf, bed_fpath, bam_fpaths, basename(__file__))


if __name__ == '__main__':
    main()