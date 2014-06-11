#!/usr/bin/env python

import sys
from source.vcf import filter_rejected

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from source.main import read_opts_and_cnfs, check_system_resources, load_genome_resources
from source.runner import run_one
from source.varannotation.tsv import make_tsv
from source.varannotation.anno import run_annotators
from source.utils import info


def main(args):
    required_keys = ['vcf']
    optional_keys = ['bam']

    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--var', '--vcf'], 'variants.vcf', {
             'dest': 'vcf',
             'help': 'variants to annotate'}),

            (['--bam'], 'align.bam', {
             'dest': 'bam',
             'help': 'used to generate some annotations by GATK'}),

            (['--clinical_reporting'], '', {
             'dest': 'clinical_reporting',
             'help': 'used to generate some annotations by GATK',
             'action': 'store_true',
             'default': False}),
        ],
        required_keys=required_keys,
        optional_keys=optional_keys,
        key_for_sample_name='bam')

    check_system_resources(cnf,
        required=['java', 'perl', 'gatk', 'snpeff', 'snpsift'],
        optional=[])
    load_genome_resources(cnf,
        required=['seq', 'snpeff'],
        optional=['dbsnp', 'cosmic'])

    set_up_snpeff(cnf)

    run_one(cnf, required_keys, optional_keys, process_one, finalize_one)


def set_up_snpeff(cnf):
    if 'clinical_reporting' in cnf:
        if 'snpeff' in cnf:
            cnf['snpeff']['clinical_reporting'] = cnf['clinical_reporting']
        del cnf['clinical_reporting']


def process_one(cnf, vcf_fpath, bam_fpath=None):
    if cnf.get('filter_reject'):
        vcf_fpath = filter_rejected(cnf, vcf_fpath)

    anno_vcf_fpath, anno_maf_fpath = run_annotators(cnf, vcf_fpath, bam_fpath)

    anno_tsv_fpath = None
    if anno_vcf_fpath and 'tsv_fields' in cnf:
        anno_tsv_fpath = make_tsv(cnf, anno_vcf_fpath)

    return anno_vcf_fpath, anno_maf_fpath, anno_tsv_fpath


def finalize_one(cnf, anno_vcf_fpath, anno_maf_fpath, anno_tsv_fpath):
    if anno_vcf_fpath:
        info('Saved final VCF to ' + anno_vcf_fpath)
    if anno_maf_fpath:
        info('Saved final MAF to ' + anno_maf_fpath)
    if anno_tsv_fpath:
        info('Saved final TSV to ' + anno_tsv_fpath)


def finalize_all(cnf, samples, results):
    for (sample_name, cnf), (vcf, maf, tsv) in zip(samples.items(), results):
        if vcf or tsv:
            info(sample_name + ':')
        if vcf:
            info('  ' + vcf)
        if maf:
            info('  ' + maf)
        if tsv:
            info('  ' + tsv)


if __name__ == '__main__':
    main(sys.argv[1:])