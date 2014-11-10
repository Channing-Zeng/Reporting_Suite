#!/usr/bin/env python
import sys
from source.bcbio_structure import BCBioStructure, Sample
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from os.path import abspath, dirname, realpath, pardir, join
from site import addsitedir
source_dir = abspath(dirname(realpath(__file__)))
addsitedir(join(source_dir, 'ext_modules'))

import shutil
from source.main import read_opts_and_cnfs, check_system_resources, load_genome_resources
from source.variants.vcf_processing import remove_rejected, extract_sample, \
     iterate_vcf, tabix_vcf, igvtools_index, get_trasncripts_fpath, fix_chromosome_names
from source.runner import run_one
from source.variants.anno import run_annotators
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

    load_genome_resources(cnf,
        required=['seq'],
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

    cnf.transcripts_fpath = get_trasncripts_fpath(cnf)


def process_one(cnf):
    sample = Sample(cnf.name, vcf=cnf.vcf, bam=cnf.bam)

    sample.vcf = fix_chromosome_names(cnf, sample.vcf)

    if cnf.get('filter_reject'):
        sample.vcf = remove_rejected(cnf, sample.vcf)
        if sample.vcf is None:
            err('No variants left for ' + cnf.vcf + ': all rejected and removed.')
            return None, None, None

    if cnf.get('extract_sample'):
        sample.vcf = extract_sample(cnf, sample.vcf, sample.name)

    anno_vcf_fpath, anno_tsv_fpath, anno_maf_fpath = run_annotators(
        cnf, sample.vcf, sample.bam, sample.name, cnf.transcript_fpath)

    info()
    info('Indexing ' + anno_vcf_fpath)
    igvtools_index(cnf, anno_vcf_fpath)

    return anno_vcf_fpath, anno_tsv_fpath, anno_maf_fpath


def finalize_one(cnf, anno_vcf_fpath, anno_tsv_fpath, anno_maf_fpath):
    msg = ['Annoatation finished for ' + cnf.name + ':']
    if anno_vcf_fpath:
        msg.append('VCF: ' + anno_vcf_fpath)
        info('Saved final VCF to ' + anno_vcf_fpath)
    if anno_tsv_fpath:
        msg.append('TSV: ' + anno_tsv_fpath)
        info('Saved final TSV to ' + anno_tsv_fpath)
    if anno_maf_fpath:
        msg.append('MAF: ' + anno_maf_fpath)
        info('Saved final MAF to ' + anno_maf_fpath)

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