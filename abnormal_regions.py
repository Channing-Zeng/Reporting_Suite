#!/usr/bin/env python
import sys
from source.bcbio_structure import Sample
from source.file_utils import adjust_path
from source.logger import send_email
from source.targetcov.abnormal_regions import make_abnormal_regions_reports
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

    load_genome_resources(
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
    sample = Sample(cnf.name, bam=cnf.bam, bed=cnf.bed)
    return make_abnormal_regions_reports(cnf, sample, cnf.vcfs_by_callername)


def finalize_one(cnf, *abnormal_regions_reports):
    msg = ['Regions with abnormal regions finished for ' + cnf.name + ':']

    if abnormal_regions_reports:
        msg.append('Abnormal region reports: ')
        info('Abnormal region reports:')
        for rep in abnormal_regions_reports:
            msg.append('  ' + rep)
            info('  ' + rep)

    send_email('\n'.join(msg))


if __name__ == '__main__':
    main(sys.argv)