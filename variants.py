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
from source.variants.vcf_processing import verify_vcf


def proc_args(argv):
    info(' '.join(sys.argv))
    info()

    description = 'This script generates target QC reports for each BAM provided as an input.'
    parser = OptionParser(description=description)
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser)
    parser.add_option('--log-dir', dest='log_dir')
    parser.add_option('--is-wgs', dest='is_wgs', action='store_true', default=False, help='whole genome sequencing')
    parser.add_option('--is-deep-seq', dest='is_deep_seq', action='store_true', default=False, help='deep targeted sequencing')
    parser.add_option('--only-summary', dest='only_summary', action='store_true')
    parser.add_option('-o', dest='output_dir', metavar='DIR', default=join(os.getcwd(), 'targetqc'))
    parser.add_option('-c', '--caller', dest='caller')
    parser.add_option('--qc-caption', dest='qc_caption', help=SUPPRESS_HELP)
    parser.add_option('--no-tsv', dest='tsv', action='store_false', default=True, help=SUPPRESS_HELP)

    (opts, args) = parser.parse_args()
    logger.is_debug = opts.debug

    if len(args) == 0:
        critical('No vcf files provided to input.')

    run_cnf = determine_run_cnf(opts, is_targetseq=opts.is_deep_seq, is_wgs=opts.is_wgs)
    cnf = Config(opts.__dict__, determine_sys_cnf(opts), run_cnf)

    sample_names, vcf_fpaths = read_samples(args, cnf.caller)
    info()

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

    run_variants(cnf, samples, basename(__file__))


def read_samples(args, caller_name=None):
    vcf_fpaths = []
    sample_names = []
    bad_vcf_fpaths = []

    info('Reading samples...')

    if len(args) == 1:
        first_fpath = args[0]
        if not first_fpath.endswith('.vcf') and not first_fpath.endswith('.vcf.gz'):  # TODO: check ##fileformat=VCF ?
            info('First argument file name does not look like VCF, assuming TSV with files names')

            with open(first_fpath) as f:
                for i, l in enumerate(f):
                    fs = l.strip().split('\t')
                    if len(fs) != 2:
                        critical('Line ' + str(i) + ' has only ' + str(len(fs)) + ' fields. Expecting 2 (sample and vcf_fpath)')
                    sn, vcf_fpath = fs
                    if not verify_file(vcf_fpath):
                        bad_vcf_fpaths.append(vcf_fpath)
                    sample_names.append(sn)
                    vcf_fpaths.append(adjust_path(vcf_fpath))

            if bad_vcf_fpaths:
                critical('VCF files cannot be found, empty or not VCFs:' + ', '.join(bad_vcf_fpaths))
            info('Done reading ' + str(len(sample_names)) + ' samples')
            return sample_names, vcf_fpaths

    for arg in args or [os.getcwd()]:
        vcf_fpath = verify_vcf(arg.split(',')[0])
        if not verify_file(vcf_fpath):
            bad_vcf_fpaths.append(vcf_fpath)
        vcf_fpaths.append(vcf_fpath)
        if len(arg.split(',')) > 1:
            sample_names.append(arg.split(',')[1])
        else:
            sn = basename(splitext_plus(vcf_fpath)[0])
            if caller_name and sn.endswith('-' + caller_name):
                sn = sn[:-len(caller_name) - 1]
            sample_names.append(sn)
            info('  ' + sn)
    if bad_vcf_fpaths:
        critical('VCF files cannot be found, empty or not VCFs:' + ', '.join(bad_vcf_fpaths))
    info('Done reading ' + str(len(sample_names)) + ' samples')

    # TODO: read sample names from VCF
    # def get_main_sample(self, main_sample_index=None):
    #     if len(self._sample_indexes) == 0:
    #         return None
    #     if main_sample_index is not None:
    #         return self.samples[main_sample_index]
    #     try:
    #         sample_index = [sname.lower() for sname in self._sample_indexes] \
    #                         .index(self.sample_name_from_file.lower())
    #     except ValueError:
    #         return self.samples[0]
    #     else:
    #         return self.samples[sample_index]

    return sample_names, vcf_fpaths


if __name__ == '__main__':
    main()
