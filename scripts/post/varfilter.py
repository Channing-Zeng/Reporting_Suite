#!/usr/bin/env python
# noinspection PyUnresolvedReferences

import bcbio_postproc


import gzip
import os
import sys
import shutil
from optparse import SUPPRESS_HELP
from os.path import dirname, join, basename

import source
from source import VarSample
from source.file_utils import iterate_file, open_gzipsafe, safe_mkdir, add_suffix
from source.main import read_opts_and_cnfs
from source.prepare_args_and_cnf import check_genome_resources, check_system_resources
from source.variants import qc
from source.variants.filtering import run_vcf2txt, run_vardict2mut, write_vcfs, write_vcf
from source.variants.vcf_processing import remove_rejected, get_sample_column_index, bgzip_and_tabix
from source.runner import run_one
from source.variants.anno import run_annotators, finialize_annotate_file
from source.utils import info
from source.logger import err, warn, critical


def main(args):
    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--vcf', '--var'], dict(
                dest='vcf',
                help='variants to filter')
             ),
            (['--min-freq'], dict(
                dest='min_freq')
             ),
            (['--qc'], dict(
                dest='qc',
                action='store_true',
                default=None,
                help=SUPPRESS_HELP)
             ),
        ],
        required_keys=['vcf'],
        file_keys=['vcf'],
        key_for_sample_name='vcf',
        proc_name=source.varfilter_name + '_post')

    check_system_resources(cnf, required=['perl'])
    check_genome_resources(cnf)

    safe_mkdir(dirname(cnf.output_file))

    vcf2txt_res_fpath = run_vcf2txt(cnf, {cnf.sample: cnf.vcf}, cnf.output_file, cnf.min_freq)
    if not vcf2txt_res_fpath:
        critical('vcf2txt run returned non-0')
        return None
    info('Saved variants to ' + vcf2txt_res_fpath)

    mut_fpath = run_vardict2mut(cnf, vcf2txt_res_fpath,
        add_suffix(vcf2txt_res_fpath, source.mut_pass_suffix))
    info('Saved passed mutations to ' + mut_fpath)

    var_s = source.VarSample(cnf.sample, cnf.output_dir)
    var_s.anno_vcf_fpath = cnf.vcf
    var_s.varfilter_dirpath = var_s.dirpath
    var_s.filt_vcf_fpath = join(cnf.output_dir, add_suffix(basename(cnf.vcf), 'filt'))
    var_s.pass_filt_vcf_fpath = add_suffix(var_s.filt_tsv_fpath, 'pass')
    var_s.varfilter_result = vcf2txt_res_fpath
    var_s.varfilter_pass_result = add_suffix(vcf2txt_res_fpath, source.mut_pass_suffix)

    write_vcf(cnf, var_s, cnf.output_dir, cnf.caller, vcf2txt_res_fpath, mut_fpath)

    if cnf.qc:
        report = qc.make_report(cnf, var_s.filt_vcf_fpath, var_s)
        qc_dirpath = join(cnf.output_dir, 'qc')
        safe_mkdir(qc_dirpath)
        qc.save_report(cnf, report, var_s, cnf.caller, qc_dirpath, source.varqc_name)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


if __name__ == '__main__':
    main(sys.argv[1:])