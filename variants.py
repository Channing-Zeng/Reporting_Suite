#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc

import sys
import os
from os.path import relpath, join, exists, abspath, pardir, basename, splitext
from optparse import OptionParser, SUPPRESS_HELP

from source import logger
import source
from source.config import Config, defaults
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug, check_genome_resources, determine_run_cnf, \
    determine_sys_cnf
from source import logger
from source.logger import info, err, warn, critical, send_email
from source.file_utils import verify_dir, safe_mkdir, adjust_path, verify_file, adjust_system_path, remove_quotes, \
    file_exists, isfile, splitext_plus
from source.main import set_up_dirs
from source.variants.variants import run_variants


def proc_args(argv):
    info(' '.join(sys.argv))
    info()

    description = 'This script generates target QC reports for each BAM provided as an input.'
    parser = OptionParser(description=description)
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser)
    parser.add_option('--log-dir', dest='log_dir')
    parser.add_option('--is-wgs', dest='is_wgs', action='store_true', default=False, help='whole genome sequencing')
    parser.add_option('--no-check', dest='no_check', action='store_true', default=False, help=SUPPRESS_HELP)
    parser.add_option('--is-deep-seq', dest='is_deep_seq', action='store_true', default=False, help='deep targeted sequencing')
    parser.add_option('--only-summary', dest='only_summary', action='store_true')
    parser.add_option('-o', dest='output_dir', metavar='DIR', default=join(os.getcwd(), 'targetqc'))

    (opts, args) = parser.parse_args()
    logger.is_debug = opts.debug

    if len(args) == 0:
        critical('No BAMs provided to input.')

    sample_names, vcf_fpaths = read_samples(args)

    run_cnf = determine_run_cnf(opts, is_targeteq=opts.is_deep_seq, is_wgs=opts.is_wgs)
    cnf = Config(opts.__dict__, determine_sys_cnf(opts), run_cnf)

    if not cnf.project_name:
        cnf.project_name = basename(cnf.output_dir)
    info('Project name: ' + cnf.project_name)

    cnf.proc_name = 'Variants'
    set_up_dirs(cnf)
    # cnf.name = 'TargQC_' + cnf.project_name
    info(' '.join(sys.argv))

    samples = [
        source.VarSample(s_name, join(cnf.output_dir, s_name), vcf=vcf_fpath)
            for s_name, vcf_fpath in zip(sample_names, vcf_fpaths)]
    samples.sort(key=lambda _s: _s.key_to_sort())

    check_genome_resources(cnf)

    if not cnf.only_summary:
        cnf.qsub_runner = adjust_system_path(cnf.qsub_runner)
        if not cnf.qsub_runner: critical('Error: qsub-runner is not provided is sys-config.')
        verify_file(cnf.qsub_runner, is_critical=True)

    return cnf, samples


def main():
    cnf, samples = proc_args(sys.argv)

    html_fpath = run_variants(cnf, samples, basename(__file__))

    # if targqc_html_fpath:
    #     send_email('TargQC report for ' + cnf.project_name + ':\n  ' + targqc_html_fpath)


def read_samples(args):
    vcf_fpaths = []
    sample_names = []
    bad_vcf_fpaths = []

    for arg in args or [os.getcwd()]:
        vcf_fpath = verify_file(arg.split(',')[0])
        if not verify_file(vcf_fpath):
            bad_vcf_fpaths.append(vcf_fpath)
        vcf_fpaths.append(vcf_fpath)
        if len(arg.split(',')) > 1:
            sample_names.append(arg.split(',')[1])
        else:
            sample_names.append(basename(splitext_plus(vcf_fpath)[0]))
    if bad_vcf_fpaths:
        critical('VCF files cannot be found, empty or not VCFs:' + ', '.join(bad_vcf_fpaths))

    return sample_names, vcf_fpaths


if __name__ == '__main__':
    main()
