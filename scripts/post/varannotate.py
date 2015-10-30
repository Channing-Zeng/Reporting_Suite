#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc


import gzip
import os
import sys
import shutil

import source
from source import VarSample
from source.file_utils import iterate_file, open_gzipsafe
from source.main import read_opts_and_cnfs
from source.prepare_args_and_cnf import check_genome_resources, check_system_resources
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

    if cnf.get('filter_reject'):
        info('Filtering rejected')
        sample.vcf = remove_rejected(cnf, sample.vcf)
        if sample.vcf is None:
            err('No variants left for ' + cnf.vcf + ': all rejected and removed.')
            return None, None, None

    # In mutect, running paired analysis on a single sample could lead
    # to a "none" sample column. Removing that column.
    info('get_sample_column_index')
    none_idx = get_sample_column_index(sample.vcf, 'none', suppress_warn=True)
    if none_idx is not None:
        info('Removing the "none" column.')
        def fn(line, i):
            if line and not line.startswith('##'):
                ts = line.split('\t')
                del ts[9 + none_idx]
                return '\t'.join(ts) + '\n'
            return line
        sample.vcf = iterate_file(cnf, sample.vcf, fn, suffix='none_col')

    # Replacing so the main sample goes first (if it is not already)
    main_idx = get_sample_column_index(sample.vcf, sample.name)
    if main_idx:
        info('Moving the main sample column (' + sample.name + ') to the first place.')
        def fn(line, i):
            if line and not line.startswith('##'):
                ts = line.split('\t')
                main_sample_field = ts[9 + main_idx]
                del ts[9 + main_idx]
                ts = ts[:9] + [main_sample_field] + ts[9:]
                return '\t'.join(ts) + '\n'
            return line
        sample.vcf = iterate_file(cnf, sample.vcf, fn, suffix='main_col')

    annotated, anno_vcf_fpath = run_annotators(cnf, sample.vcf, sample.bam)

    return finialize_annotate_file(cnf, anno_vcf_fpath, sample.name, cnf.caller)


def finalize_one(cnf, anno_vcf_fpath):
    msg = ['Annoatation finished for ' + cnf.sample + ':']
    if anno_vcf_fpath:
        msg.append('VCF: ' + anno_vcf_fpath)
        info('Saved final VCF to ' + anno_vcf_fpath)

        if is_gz(anno_vcf_fpath):
            # info(anno_vcf_fpath + ' is in correct gzip format')
            open_gzipsafe(anno_vcf_fpath)
        else:
            warn('Not a gzip:' + anno_vcf_fpath)
            anno_vcf_fpath_ungz = anno_vcf_fpath
            anno_vcf_fpath_gz = anno_vcf_fpath + '.gz'
            if anno_vcf_fpath.endswith('.gz'):
                anno_vcf_fpath_ungz = anno_vcf_fpath.split('.gz')[0]
                anno_vcf_fpath_gz = anno_vcf_fpath
                os.rename(anno_vcf_fpath_gz, anno_vcf_fpath_ungz)
            info('Compressing and indexing with bgzip+tabix again ' + anno_vcf_fpath_ungz)
            anno_vcf_fpath_gz = bgzip_and_tabix(cnf, anno_vcf_fpath_ungz)
            info('Saved VCF again to ' + anno_vcf_fpath_gz)


def is_gz(fpath, mode='rb'):
    try:
        h = gzip.open(fpath)
    except IOError, e:
        return False
    else:
        try:
            h.read(1)
        except IOError, e:
            h.close()
            return False
        else:
            h.close()
            return True


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