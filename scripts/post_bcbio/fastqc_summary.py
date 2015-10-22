#!/usr/bin/env python
import bcbio_postproc

import sys
from os.path import join

import source
from source.fastqc.summarize_fastqc import write_fastqc_combo_report
from source.bcbio.bcbio_structure import BCBioStructure, summary_script_proc_params
from source.logger import info, step_greetings


def main():
    info(' '.join(sys.argv))
    info()

    cnf, bcbio_structure = summary_script_proc_params(BCBioStructure.fastqc_name, BCBioStructure.fastqc_dir)

    step_greetings('FastQC summary for all samples')

    final_summary_report_fpath = join(cnf.output_dir, source.fastqc_name + '.html')

    write_fastqc_combo_report(final_summary_report_fpath, bcbio_structure.samples)

    info()
    info('*' * 70)
    info('Fastqc summary:')
    info('  ' + final_summary_report_fpath)


if __name__ == '__main__':
    main()
