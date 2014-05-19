#!/usr/bin/env python
import sys

from src.main import common_main, read_samples_info_and_split
from src.runner import run_all
from src.tsv import make_tsv
try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper

from src.annotation import run_annotators
from src.my_utils import info


if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    exit('Python 2, versions 2.7 and higher is supported (you are running ' +
     '.'.join(map(str, sys.version_info[:3])) + ')\n')


def main(args):
    config, options = common_main(
        'annotation',
        opts=[
            (['--var', '--vcf'], 'variants.vcf', {
             'dest': 'vcf',
             'help': 'variants to annotate'}),

            (['--bam'], 'align.bam', {
             'dest': 'bam',
             'help': 'used to generate some annotations by GATK'}),
        ])

    var_fpath = options.get('vcf')
    if var_fpath:
        print 'Using variants ' + var_fpath
    bam_fpath = options.get('bam')
    if bam_fpath:
        print 'Using bam ' + bam_fpath
    output_dir = options.get('output_dir')
    if output_dir:
        print 'Saving to ' + output_dir
    
    samples = read_samples_info_and_split(config, options)

    try:
        run_all(config, samples, process_one, finalize_one, finalize_all)
    except KeyboardInterrupt:
        exit()


def process_one(cnf, vcf_fpath):
    anno_vcf_fpath = run_annotators(cnf, vcf_fpath)
    anno_tsv_fpath = None
    if anno_vcf_fpath and 'tsv_fields' in cnf:
        anno_tsv_fpath = make_tsv(cnf, anno_vcf_fpath)
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