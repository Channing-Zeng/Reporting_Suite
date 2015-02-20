#!/usr/bin/env python

import __check_python_version

import sys
from source.bcbio_structure import BCBioStructure, summary_script_proc_params
from source.logger import info
from source.standalone_targqc.summarize import summarize_targqc


def main():
    info(' '.join(sys.argv))
    info()

    cnf, bcbio_structure = summary_script_proc_params(BCBioStructure.targqc_name, BCBioStructure.targqc_summary_dir)

    summarize_targqc(cnf, cnf.output_dir, bcbio_structure.samples, bcbio_structure.bed,
                     cnf.exons if cnf.exons else cnf.genome.exons)


if __name__ == '__main__':
    main()