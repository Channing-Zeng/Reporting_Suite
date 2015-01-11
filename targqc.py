#!/usr/bin/env python

import sub_scripts.__common  # checking for python version and adding site dirs inside

import sys
import os
from os.path import relpath, join, exists, abspath, pardir, basename
from optparse import OptionParser

from source.config import Config, defaults
from source.targetcov.summarize_targqc import summary_reports
from source.prepare_args_and_cnf import add_post_bcbio_args, detect_sys_cnf
from source.logger import info, err, warn, critical
from source.file_utils import verify_dir, safe_mkdir, adjust_path, verify_file, adjust_system_path
from source.main import check_genome_resources
from source.standalone_targqc.submit_jobs import run
from source.ngscat.bed_file import verify_bam, verify_bed


def main():
    info(' '.join(sys.argv))
    info()

    description = 'This script generates target QC reports for each BAM provided as an input.'
    parser = OptionParser(description=description)
    add_post_bcbio_args(parser)
    parser.add_option('--work-dir', dest='work_dir', metavar='DIR')
    parser.add_option('--genome', dest='genome', default='hg19')
    parser.add_option('--only-summary', dest='only_summary', action='store_true')
    parser.add_option('-o', dest='output_dir', metavar='DIR')

    (opts, args) = parser.parse_args()

    if len(args) == 0:
        critical('No BAMs provided to input.')
    bam_fpaths = [abspath(a) for a in args]
    if any(not verify_bam(fpath) for fpath in bam_fpaths):
        sys.exit(1)

    opts.sys_cnf = adjust_path(opts.sys_cnf) if opts.sys_cnf else detect_sys_cnf(opts)
    if not verify_file(opts.sys_cnf): sys.exit(1)
    info('Using ' + opts.sys_cnf)

    opts.run_cnf = adjust_path(opts.run_cnf) if opts.run_cnf else defaults['run_cnf']
    if not verify_file(opts.run_cnf): sys.exit(1)
    info('Using ' + opts.run_cnf)

    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)

    output_dir = adjust_path(cnf.output_dir or join(os.getcwd(), 'targetqc'))
    if not verify_dir(join(output_dir, pardir)): sys.exit(1)
    safe_mkdir(output_dir)
    info('Output to ' + output_dir)
    cnf.output_dir = output_dir

    cnf.qsub_runner = adjust_system_path(cnf.qsub_runner)
    if not cnf.qsub_runner: critical('Error: qsub-runner is not provided is sys-config.')
    if not verify_file(cnf.qsub_runner): sys.exit(1)

    check_genome_resources(cnf)

    if not cnf.work_dir: cnf.work_dir = join(cnf.output_dir, 'work')
    safe_mkdir(cnf.work_dir)
    cnf.log_dir = join(cnf.work_dir, 'log')
    safe_mkdir(cnf.log_dir)

    bed_fpath = cnf.bed
    if not verify_bed(bed_fpath):
        sys.exit(1)

    if not cnf.project_name: cnf.project_name = basename(cnf.output_dir)

    info('*' * 70)
    info()

    run(cnf, bed_fpath, bam_fpaths, basename(__file__))


if __name__ == '__main__':
    main()