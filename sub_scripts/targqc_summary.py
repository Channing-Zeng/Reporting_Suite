#!/usr/bin/env python

import __check_python_version

import sys
from os.path import abspath
from source.bcbio_structure import BCBioStructure, summary_script_proc_params
from source.file_utils import adjust_path
from source.logger import info
from source.standalone_targqc.summarize import summarize_targqc


def main():
    info(' '.join(sys.argv))
    info()

    cnf, bcbio_structure = summary_script_proc_params(BCBioStructure.targqc_name, BCBioStructure.targqc_summary_dir)

    exons_bed_fpath = adjust_path(cnf.exons) if cnf.exons else adjust_path(cnf.genome.exons)
    info('Exons: ' + exons_bed_fpath)

    bed_fpath = bcbio_structure.bed or cnf.genome.az_exome or exons_bed_fpath
    info('Using amplicons/capture panel ' + abspath(bed_fpath))

    genes_fpath = None
    if cnf.genes:
        genes_fpath = adjust_path(cnf.genes)
        info('Custom genes list: ' + genes_fpath)

    summarize_targqc(cnf, cnf.threads or len(bcbio_structure.samples),
        cnf.output_dir, bcbio_structure.samples, bed_fpath, exons_bed_fpath, genes_fpath)


if __name__ == '__main__':
    main()