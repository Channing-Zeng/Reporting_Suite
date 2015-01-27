#!/usr/bin/env python

from source.bcbio_structure import BCBioStructure
from source.prepare_args_and_cnf import summary_script_proc_params
from source.copy_number import cnv_reports


def main():
    cnf, bcbio_structure = summary_script_proc_params(BCBioStructure.seq2c_name, BCBioStructure.cnv_summary_dir)

    cnv_reports(cnf, bcbio_structure)


if __name__ == '__main__':
    main()

