#!/usr/bin/env python
import __check_python_version

import sys
import shutil

from source import BaseSample, info, verify_file
from source.clinical_reporting.clinical_reporting import make_key_genes_reports, make_mutations_report
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

    sample = BaseSample(cnf.sample, cnf.output_dir, targqc_dirpath=cnf.targqc_dirpath, clinical_report_dirpath=cnf.output_dir)

    info('Building clinical report for AZ 300 key genes ' + str(cnf.key_genes))
    cnf.key_genes = verify_file(cnf.key_genes, is_critical=True, description='300 AZ key genes')
    with open(cnf.key_genes) as f:
        key_gene_names = set([l.strip() for l in f.readlines()])

    key_genes_report = make_key_genes_reports(cnf, key_gene_names, sample)

    mutations_report = make_mutations_report(cnf, key_gene_names, sample, cnf.mutations_fpath)

    _finalize_all(key_genes_report, mutations_report)

    '''
    TODO:
    6. Intersect with Brian's genes, make that report
    7. Gender
    8. Make the total HTML
    '''

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


def _finalize_all(key_genes_report, mutations_report):
    info()
    if key_genes_report and key_genes_report.tsv_fpath:
        info('Coverage stats: ' + key_genes_report.tsv_fpath)

    if mutations_report and mutations_report.tsv_fpath:
        info('Mutations stats: ' + mutations_report.tsv_fpath)


if __name__ == '__main__':
    main(sys.argv)