#!/usr/bin/env python
import sys
from source.vcf_read import read_samples_info_and_split

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from source.main import common_main, check_system_resources, load_genome_resources
from source.runner import run_all
from source.varannotation import tsv, anno
from source.utils import info


def main(args):
    required = ['vcf']
    optional = ['bam']

    config, options = common_main(
        'annotation',
        extra_opts=[
            (['--var', '--vcf'], 'variants.vcf', {
             'dest': 'vcf',
             'help': 'variants to annotate'}),

            (['--bam'], 'align.bam', {
             'dest': 'bam',
             'help': 'used to generate some annotations by GATK'}),

            (['--clinical_reporting'], {
             'dest': 'clinical_reporting',
             'help': 'used to generate some annotations by GATK',
             'action': 'set_true',
             'default': False}),
        ],
        required=required)

    if 'clinical_reporting' in options and 'snpeff' in config:
        config['snpeff']['clinical_reporting'] = True

    check_system_resources(config, ['java', 'perl', 'gatk', 'snpeff', 'snpsift'])
    load_genome_resources(config, ['seq', 'dbsnp', 'cosmic', 'snpeff'])

    var_fpath = options.get('vcf')
    if var_fpath:
        print 'Using variants ' + var_fpath
    bam_fpath = options.get('bam')
    if bam_fpath:
        print 'Using bam ' + bam_fpath
    output_dir = options.get('output_dir')
    if output_dir:
        print 'Saving to ' + output_dir

    sample_cnfs_by_name = read_samples_info_and_split(config, options, required + optional)

    try:
        run_all(config, sample_cnfs_by_name, required, optional,
                process_one, finalize_one, finalize_all)
    except KeyboardInterrupt:
        exit()


def process_one(cnf, vcf_fpath, bam_fpath=None):
    anno_vcf_fpath = anno.run_annotators(cnf, vcf_fpath, bam_fpath)
    anno_tsv_fpath = None
    if anno_vcf_fpath and 'tsv_fields' in cnf:
        anno_tsv_fpath = tsv.make_tsv(cnf, anno_vcf_fpath)
    return anno_vcf_fpath, anno_tsv_fpath


def finalize_one(cnf, anno_vcf_fpath, anno_tsv_fpath):
    if anno_vcf_fpath:
        info(cnf['log'], 'Saved final VCF to ' + anno_vcf_fpath)
    if anno_tsv_fpath:
        info(cnf['log'], 'Saved final TSV to ' + anno_tsv_fpath)


def finalize_all(cnf, samples, results):
    for (sample_name, cnf), (vcf, tsv) in zip(samples.items(), results):
        if vcf or tsv:
            info(cnf['log'], sample_name + ':')
        if vcf:
            info(cnf['log'], '  ' + vcf)
        if tsv:
            info(cnf['log'], '  ' + tsv)


if __name__ == '__main__':
    main(sys.argv[1:])