#!/usr/bin/env python

import __common

import sys
import shutil
from source.bcbio_structure import BCBioStructure, Sample
from source.file_utils import iterate_file
from source.main import read_opts_and_cnfs, check_system_resources, check_genome_resources
from source.variants.vcf_processing import remove_rejected, extract_sample, fix_chromosome_names, iterate_vcf, \
    read_sample_names_from_vcf, get_sample_column_index
from source.runner import run_one
from source.variants.anno import run_annotators, finialize_annotate_file
from source.utils import info
from source.logger import err, send_email


def main(args):
    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--vcf', '--var'], dict(
                dest='vcf',
                help='variants to annotate')
             ),
            (['--bam'], dict(
                dest='bam',
                help='used to generate some annotations by GATK')
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
        proc_name=BCBioStructure.varannotate_name)

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
    sample = Sample(cnf.name, vcf=cnf.vcf, bam=cnf.bam, genome=cnf.genome)

    sample.vcf = fix_chromosome_names(cnf, sample.vcf)

    if cnf.get('filter_reject'):
        sample.vcf = remove_rejected(cnf, sample.vcf)
        if sample.vcf is None:
            err('No variants left for ' + cnf.vcf + ': all rejected and removed.')
            return None, None, None

    if cnf.get('extract_sample'):
        sample.vcf = extract_sample(cnf, sample.vcf, sample.name)

    # In mutect, running paired analysis on a single sample could lead
    # to a "none" sample column. Removing that column.
    none_idx = get_sample_column_index(sample.vcf, 'none', suppress_warn=True)
    if none_idx is not None:
        info('Removing the "none" column.')
        def fn(line, i):
            if line and not line.startswith('##'):
                ts = line.split('\t')
                del ts[9 + none_idx]
                return '\t'.join(ts) + '\n'
            return line
        sample.vcf = iterate_file(cnf, sample.vcf, fn, 'none_col')

    # Replacing so the main sample goes first (if it is not already)
    main_idx = get_sample_column_index(sample.vcf, sample.name)
    if main_idx != 0:
        info('Moving the main sample column (' + sample.name + ') to the first place.')
        def fn(line, i):
            if line and not line.startswith('##'):
                ts = line.split('\t')
                main_sample_field = ts[9 + main_idx]
                del ts[9 + main_idx]
                ts = ts[:9] + [main_sample_field] + ts[9:]
                return '\t'.join(ts) + '\n'
            return line
        sample.vcf = iterate_file(cnf, sample.vcf, fn, 'main_col')

    annotated, anno_vcf_fpath = run_annotators(cnf, sample.vcf, sample.bam)

    anno_vcf_fpath, anno_tsv_fpath = finialize_annotate_file(cnf, anno_vcf_fpath, sample.name, cnf.caller)

    return anno_vcf_fpath, anno_tsv_fpath


def finalize_one(cnf, anno_vcf_fpath, anno_tsv_fpath):
    msg = ['Annoatation finished for ' + cnf.name + ':']
    if anno_vcf_fpath:
        msg.append('VCF: ' + anno_vcf_fpath)
        info('Saved final VCF to ' + anno_vcf_fpath)
    if anno_tsv_fpath:
        msg.append('TSV: ' + anno_tsv_fpath)
        info('Saved final TSV to ' + anno_tsv_fpath)

    # send_email('\n'.join(msg))


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