#!/usr/bin/env python
import __check_python_version

import sys
import shutil

from source import BaseSample, info
from source.clinical_reporting.clinical_reporting import make_key_genes_reports
from source.prepare_args_and_cnf import check_system_resources
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
        ],
        required_keys=['targqc_dirpath', 'mutations_fpath'],
        file_keys=['mutations_fpath'],
        dir_keys=['targqc_dirpath'],
        key_for_sample_name=None
    )

    check_system_resources(
        cnf,
        required=['samtools', 'bedtools'],
        optional=[])

    sample = BaseSample(cnf.sample, cnf.output_dir, targqc_dirpath=cnf.targqc_dirpath)

    key_genes_report = make_key_genes_reports(cnf, sample)
    if key_genes_report and key_genes_report.tsv_fpath:
        info('Saved to ' + key_genes_report.tsv_fpath)

    '''
    TODO:
    3. Parse mutaions
    4. Select mutations
    5. Make mutaions report
    6. Intersect with Brian's genes, make that report
    7. Gender
    8. Make the total HTML
    '''

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


if __name__ == '__main__':
    main(sys.argv)