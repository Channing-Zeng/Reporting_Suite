#!/usr/bin/env python
# noinspection PyUnresolvedReferences
from os.path import splitext

import bcbio_postproc


import gzip
import os
import sys
import shutil
from optparse import SUPPRESS_HELP

import source
from source import VarSample
from source.calling_process import call
from source.file_utils import iterate_file, open_gzipsafe, intermediate_fname
from source.main import read_opts_and_cnfs
from source.prepare_args_and_cnf import check_genome_resources, check_system_resources
from source.tools_from_cnf import get_system_path
from source.variants.vcf_processing import remove_rejected, get_sample_column_index, verify_vcf
from source.runner import run_one
from source.variants.anno import run_annotators, finialize_annotate_file
from source.utils import info
from source.logger import err, warn


def main(args):
    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--vcf', '--var'], dict(
                dest='vcf',
                help='variants to annotate')
             ),
            (['--bam'], dict(
                dest='bam',
                help='(outdated) used to generate some annotations by GATK')
             ),
            (['--match-normal-sample-name'], dict(
                dest='match_normal_normal_name')
             ),
            (['--clinical-reporting'], dict(
                dest='clinical_reporting',
                help='used to generate some annotations by GATK',
                action='store_true',
                default=None)
             ),
            (['--qc'], dict(
                dest='qc',
                action='store_true',
                default=None,
                help=SUPPRESS_HELP)
             ),
            (['--transcripts'], dict(
                dest='transcripts_fpath')
             ),
        ],
        required_keys=['vcf'],
        file_keys=['bam', 'vcf'],
        key_for_sample_name='vcf',
        proc_name=source.varannotate_name)

    check_system_resources(cnf,
        required=['java', 'perl', 'snpeff', 'snpsift'],
        optional=['transcripts_fpath'])

    check_genome_resources(cnf)

    # info('Using variants ' + cnf['vcf'])
    # info('Using alignement ' + cnf['bam'])

    run_one(cnf, process_one, finalize_one)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


def process_one(cnf):
    info('process_one')
    sample = VarSample(cnf.sample, cnf.output_dir, vcf=cnf.vcf, bam=cnf.bam, genome=cnf.genome)

    # this method will also gunzip the vcf file
    # sample.vcf = fix_chromosome_names(cnf, sample.vcf)
    # if cnf.vcf.endswith('.gz'):
    #     vcf_fpath = intermediate_fname(cnf, splitext(sample.vcf)[0], None)
    #     info('Ungzipping ' + sample.vcf + ', writing to ' + vcf_fpath)
    #     gunzip = get_system_path(cnf, 'gunzip', is_critical=True)
    #     cmdl = '{gunzip} {sample.vcf} --to-stdout'.format(**locals())
    #     call(cnf, cmdl, output_fpath=vcf_fpath)
    #     verify_vcf(vcf_fpath)
    #     sample.vcf = vcf_fpath

    ungz_pass_vcf_fpath = remove_rejected(cnf, cnf.vcf, intermediate_fname(cnf, cnf.vcf, ''))
    info()

    # if sample.vcf is None:
    #     err('No variants left for ' + cnf.vcf + ': all rejected and removed.')
    #     return None, None, None

    # # In mutect, running paired analysis on a single sample could lead
    # # to a "none" sample column. Removing that column.
    # info('get_sample_column_index')
    # none_idx = get_sample_column_index(sample.vcf, 'none', suppress_warn=True)
    # if none_idx is not None:
    #     info('Removing the "none" column.')
    #     def fn(line, i):
    #         if line and not line.startswith('##'):
    #             ts = line.split('\t')
    #             del ts[9 + none_idx]
    #             return '\t'.join(ts) + '\n'
    #         return line
    #     sample.vcf = iterate_file(cnf, sample.vcf, fn, suffix='none_col')

    # Replacing so the main sample goes first (if it is not already)
    # main_idx = get_sample_column_index(sample.vcf, sample.name)
    # if main_idx:
    #     info('Moving the main sample column (' + sample.name + ') to the first place.')
    #     def fn(line, i):
    #         if line and not line.startswith('##'):
    #             ts = line.split('\t')
    #             main_sample_field = ts[9 + main_idx]
    #             del ts[9 + main_idx]
    #             ts = ts[:9] + [main_sample_field] + ts[9:]
    #             return '\t'.join(ts) + '\n'
    #         return line
    #     sample.vcf = iterate_file(cnf, sample.vcf, fn, suffix='main_col')

    anno_vcf_fpath = run_annotators(cnf, ungz_pass_vcf_fpath, sample.bam)

    return finialize_annotate_file(cnf, anno_vcf_fpath, sample, cnf.caller)


def finalize_one(cnf, anno_vcf_fpath):
    msg = ['Annoatation finished for ' + cnf.sample + ':']
    if anno_vcf_fpath:
        msg.append('VCF: ' + anno_vcf_fpath)
        info('Saved final VCF to ' + anno_vcf_fpath)


def finalize_all(cnf, samples, results):
    for (sample_name, cnf), (vcf, tsv, maf) in zip(samples.items(), results):
        if vcf or tsv:
            info(sample_name + ':')
        if vcf:
            info('  ' + vcf)
        if tsv:
            info('  ' + tsv)
        if maf:
            info('  ' + maf)


if __name__ == '__main__':
    main(sys.argv[1:])