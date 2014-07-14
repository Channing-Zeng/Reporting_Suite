#!/usr/bin/env python
import shutil
import sys

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from source.main import read_opts_and_cnfs, check_system_resources, load_genome_resources
from source.variants.vcf_processing import filter_rejected, extract_sample
from source.runner import run_one
from source.variants.tsv import make_tsv
from source.variants.anno import run_annotators
from source.utils import info


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
            (['--clinical_reporting'], dict(
                dest='clinical_reporting',
                help='used to generate some annotations by GATK',
                action='store_true',
                default=None)
             ),
        ],
        required_keys=['vcf'],
        file_keys=['bam', 'vcf'],
        key_for_sample_name='vcf')

    check_system_resources(cnf,
        required=['java', 'perl', 'gatk', 'snpeff', 'snpsift'],
        optional=[])

    load_genome_resources(cnf,
        required=['seq', 'snpeff'],
        optional=['dbsnp', 'cosmic', 'oncomine'])

    set_up_snpeff(cnf)

    # info('Using variants ' + cnf['vcf'])
    # info('Using alignement ' + cnf['bam'])

    run_one(cnf, process_one, finalize_one)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


def set_up_snpeff(cnf):
    if cnf.get('clinical_reporting') is not None:
        if 'snpeff' in cnf:
            cnf['snpeff']['clinical_reporting'] = cnf['clinical_reporting']
        del cnf['clinical_reporting']


def process_one(cnf):
    vcf_fpath = cnf['vcf']
    bam_fpath = cnf['bam']

    if cnf.get('filter_reject'):
        vcf_fpath = filter_rejected(cnf, vcf_fpath)

    if cnf.get('extract_sample'):
        vcf_fpath = extract_sample(cnf, vcf_fpath, cnf['name'])

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