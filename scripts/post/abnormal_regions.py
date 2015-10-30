#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc

import sys
import shutil
from source import BaseSample
from source.file_utils import adjust_path
from source.prepare_args_and_cnf import check_system_resources
from source.main import read_opts_and_cnfs
from source.targetcov.flag_regions import generate_flagged_regions_report
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
            # (['--region-report'], dict(
            #     dest='vcfs',
            #     help='filteted variants in VCF, comma-separate, must correspond to caller-names')
            #  ),
        ],
        required_keys=[],
        file_keys=[],
        key_for_sample_name=None
    )

    check_system_resources(
        cnf,
        required=['samtools', 'bedtools'],
        optional=[])

    cnf.vcfs = map(adjust_path, cnf.vcfs.split(',') if cnf.vcfs else [])
    cnf.caller_names = cnf.caller_names.split(',') if cnf.caller_names else []
    cnf.vcfs_by_callername = zip(cnf.caller_names, cnf.vcfs)

    process_all(cnf.output_dir)


    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


def process_all(targetcov_dirpath):
    generate_flagged_regions_report()
    pass
    # read all detail reports
    # normalize
    # extrac low-cov
    # report cov and missed vars


class Sample(BaseSample):
    def __init__(self, name, output_dir, **kwargs):
        BaseSample.__init__(self, name, output_dir, **kwargs)


# def process_one(cnf, output_dir):
#     sample = Sample(cnf.name, output_dir, bam=cnf.bam, bed=cnf.bed)
#     return make_flagged_regions_reports(cnf, output_dir, sample)


def finalize_one(cnf, *abnormal_regions_reports):
    msg = ['Regions with abnormal regions finished for ' + cnf.sample + ':']

    if abnormal_regions_reports:
        msg.append('Abnormal region reports: ')
        info('Abnormal region reports:')
        for rep in abnormal_regions_reports:
            msg.append('  ' + rep)
            info('  ' + rep)

    # send_email('\n'.join(msg))


if __name__ == '__main__':
    main(sys.argv)