#!/usr/bin/env python
import __check_python_version

import sys
import shutil

from source import BaseSample, info, verify_file
from source.clinical_reporting.clinical_reporting import make_key_gene_cov_report, make_mutations_report, \
    make_clinical_html_report, get_target_fraction, get_ave_coverage, get_gender, get_total_variants_number, \
    get_min_coverage
from source.file_utils import verify_module
from source.logger import warn
from source.prepare_args_and_cnf import check_system_resources, check_genome_resources
from source.main import read_opts_and_cnfs
from tools import seq2c_plots


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
            (['--seq2c'], dict(
                dest='seq2c_tsv_fpath',
            )),
        ],
        required_keys=['targqc_dirpath', 'mutations_fpath', 'varqc_json_fpath', 'seq2c_tsv_fpath'],
        file_keys=['mutations_fpath', 'varqc_json_fpath', 'varqc_after_json_fpath', 'seq2c_tsv_fpath'],
        dir_keys=['targqc_dirpath'],
        key_for_sample_name=None
    )

    check_genome_resources(cnf)
    check_system_resources(
        cnf,
        required=['samtools', 'bedtools'],
        optional=[])

    sample = BaseSample(cnf.sample, cnf.output_dir, targqc_dirpath=cnf.targqc_dirpath,
                        clinical_report_dirpath=cnf.output_dir)

    info('Building clinical report for AZ 300 key genes ' + str(cnf.key_genes))
    cnf.key_genes = verify_file(cnf.key_genes, is_critical=True, description='300 AZ key genes')
    with open(cnf.key_genes) as f:
        key_gene_names = set([l.strip() for l in f.readlines() if l.strip() != ''])

    ave_depth = get_ave_coverage(sample, sample.targetcov_json_fpath)
    target_fraction = get_target_fraction(sample, sample.targetcov_json_fpath)
    gender = get_gender(sample, sample.targetcov_json_fpath)
    total_variants = get_total_variants_number(sample, cnf.varqc_json_fpath)

    key_genes_report = make_key_gene_cov_report(cnf, sample, key_gene_names, ave_depth)
    mutations_report = make_mutations_report(cnf, sample, key_gene_names, cnf.mutations_fpath)

    seq2c_plot_fpath = None
    if not cnf.seq2c_tsv_fpath:
        warn('No Seq2C results provided by option --seq2c, skipping plotting Seq2C')
    else:
        if not verify_module('matplotlib'):
            warn('No matplotlib, skipping plotting Seq2C')
        else:
            seq2c_plot_fpath = seq2c_plots.draw_seq2c_plot(cnf, cnf.seq2c_tsv_fpath, sample.name, cnf.output_dir, key_gene_names)

    html_fpath = make_clinical_html_report(cnf, sample, key_genes_report, mutations_report,
        cnf.target_type, ave_depth, target_fraction, gender, total_variants,
        key_gene_names, seq2c_plot_fpath)

    info('Clinical report: ' + html_fpath)

    '''
    TODO:
    6. Intersect with Brian's genes, make that report
    7. Gender
    8. Make the total HTML
    '''

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


if __name__ == '__main__':
    main(sys.argv)