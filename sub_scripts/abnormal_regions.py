#!/usr/bin/env python

import __common

import sys
import shutil
from source.bcbio_structure import BCBioSample
from source.file_utils import adjust_path
from source.logger import send_email
from source.targetcov.abnormal_regions import make_abnormal_regions_reports
from source.main import read_opts_and_cnfs, check_system_resources, check_genome_resources
from source.runner import run_one
from source.utils import info


def main(args):
    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--caller-names'], dict(
                dest='caller_names',
                help='names of variant callers used to create vcfs provided by --vcf')
             ),
            (['--vcfs'], dict(
                dest='vcfs',
                help='filteted variants in VCF, comma-separate, must correspond to caller-names')
             ),
            (['--region-report'], dict(
                dest='vcfs',
                help='filteted variants in VCF, comma-separate, must correspond to caller-names')
             ),
        ],
        required_keys=[],
        file_keys=[],
        key_for_sample_name=None
        )

    check_system_resources(
        cnf,
        required=['samtools', 'bedtools'],
        optional=[])

    check_genome_resources(
        cnf,
        required=['seq', 'exons'],
        optional=['chr_lengths', 'genes'])

    cnf.vcfs = map(adjust_path, cnf.vcfs.split(',') if cnf.vcfs else [])
    cnf.caller_names = cnf.caller_names.split(',') if cnf.caller_names else []
    cnf.vcfs_by_callername = zip(cnf.caller_names, cnf.vcfs)

    run_one(cnf, process_one, finalize_one)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


def process_one(cnf):
    sample = BCBioSample(cnf.name, bam=cnf.bam, bed=cnf.bed)
    return make_abnormal_regions_reports(cnf, sample, cnf.vcfs_by_callername)


def finalize_one(cnf, *abnormal_regions_reports):
    msg = ['Regions with abnormal regions finished for ' + cnf.name + ':']

    if abnormal_regions_reports:
        msg.append('Abnormal region reports: ')
        info('Abnormal region reports:')
        for rep in abnormal_regions_reports:
            msg.append('  ' + rep)
            info('  ' + rep)

    # send_email('\n'.join(msg))


if __name__ == '__main__':
    main(sys.argv)