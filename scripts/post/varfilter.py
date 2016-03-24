#!/usr/bin/env python
# noinspection PyUnresolvedReferences

import bcbio_postproc


import gzip
import os
import sys
import shutil
from optparse import SUPPRESS_HELP
from os.path import dirname, join, basename, splitext

import source
from source import VarSample
from source.file_utils import iterate_file, open_gzipsafe, safe_mkdir, add_suffix, file_transaction, verify_file
from source.main import read_opts_and_cnfs
from source.prepare_args_and_cnf import check_genome_resources, check_system_resources
from source.tools_from_cnf import get_script_cmdline
from source.variants import qc
from source.variants.filtering import run_vcf2txt, run_vardict2mut, write_vcfs, write_vcf, index_vcf
from source.variants.vcf_processing import remove_rejected, get_sample_column_index, bgzip_and_tabix, verify_vcf
from source.runner import run_one
from source.variants.anno import run_annotators, finialize_annotate_file
from source.utils import info, is_local
from source.logger import err, warn, critical


def main(args):
    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--vcf', '--var'], dict(
                dest='vcf',
                help='variants to filter')
             ),
            (['--vcf2txt'], dict(
                dest='vcf2txt',
                help='variants in vcf2txt to filter')
             ),
            (['--qc'], dict(
                dest='qc',
                action='store_true',
                default=None,
                help=SUPPRESS_HELP)
             ),
            (['--no-tsv'], dict(
                dest='tsv',
                action='store_false',
                default=True,
                help=SUPPRESS_HELP)
             ),
            (['--freq', '--min-freq'], dict(
                dest='min_freq',
                type='float',
                help=SUPPRESS_HELP)
             ),
        ],
        required_keys=['vcf'],
        file_keys=['vcf'],
        key_for_sample_name='vcf',
        proc_name=source.varfilter_name + '_post')

    check_system_resources(cnf, required=['perl'])
    check_genome_resources(cnf)

    if not cnf.output_file:
        cnf.output_file = join(cnf.output_dir, (cnf.caller or 'variants') + '.txt')

    safe_mkdir(dirname(cnf.output_file))
    safe_mkdir(cnf.output_dir)

    if cnf.vcf.endswith('.vcf.gz') or cnf.vcf.endswith('.vcf'):
        verify_vcf(cnf.vcf, is_critical=True)

    if not cnf.vcf2txt:
        vcf2txt_res_fpath = run_vcf2txt(cnf, {cnf.sample: cnf.vcf}, cnf.output_file, cnf.min_freq)
        if not vcf2txt_res_fpath:
            critical('vcf2txt run returned non-0')
        info('Saved vcf2txt output to ' + vcf2txt_res_fpath)
    else:
        cnf.vcf2txt = verify_file(cnf.vcf2txt, is_critical=True)
        info('Input is vcf2txt output, grepping by sample name ' + cnf.sample)
        vcf2txt_res_fpath = cnf.output_file
        with file_transaction(cnf.work_dir, vcf2txt_res_fpath) as tx:
            with open(cnf.vcf2txt) as f, open(tx, 'w') as out:
                for i, l in enumerate(f):
                    if l.strip():
                        if i == 0:
                            out.write(l)
                        else:
                            if l.split('\t')[0] == cnf.sample:
                                out.write(l)
        info('Using vcf2txt from ' + vcf2txt_res_fpath)

    if is_local():
        vardict2mut_pl = get_script_cmdline(cnf, 'perl', join('VarDict', 'vardict2mut.pl'))
        info('Running vardict2mut perl')
        res = run_vardict2mut(cnf, vcf2txt_res_fpath,
            add_suffix(vcf2txt_res_fpath, source.mut_pass_suffix + '_perl'),
            vardict2mut_executable=vardict2mut_pl)
        if not res:
            critical('vardict2mut.pl run returned non-0')

    mut_fpath = run_vardict2mut(cnf, vcf2txt_res_fpath, add_suffix(vcf2txt_res_fpath, source.mut_pass_suffix))
    if not mut_fpath:
        err('vardict2mut failed')
    else:
        info('Saved passed mutations to ' + mut_fpath)

        var_s = source.VarSample(cnf.sample, cnf.output_dir)
        var_s.anno_vcf_fpath = cnf.vcf
        var_s.varfilter_dirpath = var_s.dirpath

        ungz_anno_vcf_fpath = var_s.anno_vcf_fpath if not var_s.anno_vcf_fpath.endswith('.gz') else splitext(var_s.anno_vcf_fpath)[0]
        ungz_filt_vcf_fpath = join(cnf.output_dir, add_suffix(basename(ungz_anno_vcf_fpath), 'filt'))
        var_s.filt_vcf_fpath = ungz_filt_vcf_fpath + '.gz'
        var_s.pass_filt_vcf_fpath = add_suffix(ungz_filt_vcf_fpath, 'pass')

        var_s.variants_fpath = vcf2txt_res_fpath
        var_s.variants_pass_fpath = add_suffix(vcf2txt_res_fpath, source.mut_pass_suffix)

        filt_vcf = write_vcf(cnf, var_s, cnf.output_dir, cnf.caller, vcf2txt_res_fpath, mut_fpath)
        index_vcf(cnf, var_s.name, var_s.pass_filt_vcf_fpath, filt_vcf, cnf.caller)

        if cnf.qc:
            report = qc.make_report(cnf, var_s.pass_filt_vcf_fpath, var_s)
            qc_dirpath = join(cnf.output_dir, 'qc')
            safe_mkdir(qc_dirpath)
            qc.save_report(cnf, report, var_s, cnf.caller, qc_dirpath, source.varqc_after_name)
            info('Saved QC to ' + qc_dirpath + ' (' + report.html_fpath + ')')
            info('-' * 70)
            info()

        if not cnf['keep_intermediate']:
            shutil.rmtree(cnf['work_dir'])

        info()
        info('*' * 70)
        info('Done filtering ' + var_s.name)


if __name__ == '__main__':
    main(sys.argv[1:])