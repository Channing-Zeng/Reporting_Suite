#!/usr/bin/env python
# noinspection PyUnresolvedReferences
from os.path import dirname

import bcbio_postproc


import gzip
import os
import sys
import shutil
from optparse import SUPPRESS_HELP

import source
from source import VarSample
from source.file_utils import iterate_file, open_gzipsafe, safe_mkdir
from source.main import read_opts_and_cnfs
from source.prepare_args_and_cnf import check_genome_resources, check_system_resources
from source.variants.filtering import prep_vcf, run_vcf2txt
from source.variants.vcf_processing import remove_rejected, get_sample_column_index, bgzip_and_tabix
from source.runner import run_one
from source.variants.anno import run_annotators, finialize_annotate_file
from source.utils import info
from source.logger import err, warn


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
        ],
        required_keys=['vcf'],
        file_keys=['vcf'],
        key_for_sample_name='vcf',
        proc_name=source.varfilter_name + '_post',
        with_output_file=True)

    check_system_resources(cnf, required=['perl'])
    check_genome_resources(cnf)

    info('Preparing VCFs for vcf2txt (filtering PASS only)')
    prep_vcf_fpath = prep_vcf(cnf, cnf.vcf, cnf.sample, cnf.caller)
    info('Prepared VCF: ' + prep_vcf_fpath)

    safe_mkdir(dirname(cnf.output_file))
    res = run_vcf2txt(cnf, {cnf.sample: prep_vcf_fpath}, cnf.output_file, cnf.min_freq)
    if not res:
        err('vcf2txt run returned non-0')
        return None

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


if __name__ == '__main__':
    main(sys.argv[1:])