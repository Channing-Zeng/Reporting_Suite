#!/usr/bin/env python

import __check_python_version

import sys
from source.fastqc.summarize_fastqc import summary_reports
from source.bcbio_structure import BCBioStructure
from source.prepare_args_and_cnf import summary_script_proc_params
from source.logger import info


def main():
    info(' '.join(sys.argv))
    info()

    cnf, bcbio_structure = summary_script_proc_params(BCBioStructure.fastqc_name, BCBioStructure.fastqc_dir)

    summary_reports(cnf, bcbio_structure)


if __name__ == '__main__':
    main()
