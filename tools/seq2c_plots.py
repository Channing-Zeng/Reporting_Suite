#!/usr/bin/env python
import __check_python_version

from source import info
from source.clinical_reporting.seq2c_plot import draw_seq2c_plot
from source.main import read_opts_and_cnfs
from source.prepare_args_and_cnf import check_system_resources
from source.prepare_args_and_cnf import check_genome_resources
cnv_plot_ending = '.cnv.png'


def main():
    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--seq2c-results'], dict(
                dest='seq2c_tsv_fpath')
             ),
            (['--key-genes'], dict(
                dest='key_genes_fpath')
             ),
        ],
        required_keys=['seq2c_tsv_fpath', 'key_genes'],
        file_keys=['seq2c_tsv_fpath', 'key_genes'],
        key_for_sample_name=None,
    )
    check_system_resources(cnf)
    check_genome_resources(cnf)

    with open(cnf.key_genes_fpath) as f:
        key_gene_names = set([l.strip() for l in f.readlines() if l.strip() != ''])

    plot_fpath = draw_seq2c_plot(cnf, cnf.seq2c_tsv_fpath, cnf.sample, cnf.output_dir, key_gene_names)
    if plot_fpath:
        info('Saved plot to ' + plot_fpath)


if __name__ == '__main__':
    main()