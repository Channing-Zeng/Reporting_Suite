#!/usr/bin/env python
import bcbio_postproc

import sys
import shutil

import source
from source.clinical_reporting.clinical_reporting import run_sample_clinical_reporting
from source.file_utils import verify_module, verify_file
from source.logger import warn, info
from source.prepare_args_and_cnf import check_system_resources, check_genome_resources
from source.main import read_opts_and_cnfs


def main(args):
    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--targqc-dir'], dict(
                dest='targqc_dirpath',
            )),
            (['--mutations'], dict(
                dest='mutations_fpath',
            )),
            (['--varqc'], dict(
                dest='varqc_json_fpath',
            )),
            (['--varqc-after'], dict(
                dest='varqc_after_json_fpath',
            )),
            (['--target-type'], dict(
                dest='target_type',
            )),
            (['--bed'], dict(
                dest='bed_fpath',
            )),
            (['--seq2c'], dict(
                dest='seq2c_tsv_fpath',
            )),
            (['--project-level-report'], dict(
                dest='project_level_report_fpath',
            )),
            (['--match'], dict(
                dest='match_sample_name',
            )),
        ],
        required_keys=['targqc_dirpath', 'mutations_fpath', 'varqc_json_fpath'],
        file_keys=['varqc_json_fpath', 'varqc_after_json_fpath', 'seq2c_tsv_fpath'],  # do not check mutations_fpath! could be vardict.PASS.txt, vardict-java.PASS.txt, vardict.single.PASS.txt, vardict.paired.PASS.txt,
        dir_keys=['targqc_dirpath'],
        key_for_sample_name=None
    )

    check_genome_resources(cnf)
    check_system_resources(
        cnf,
        required=['samtools', 'bedtools'],
        optional=[])

    html_fpath = run_sample_clinical_reporting(cnf)
    info('Clinical report: ' + html_fpath)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


if __name__ == '__main__':
    main(sys.argv)