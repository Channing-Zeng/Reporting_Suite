#!/usr/bin/env python

import __check_python_version

import sys

from source.bcbio_structure import BCBioStructure, summary_script_proc_params
from source.copy_number import cnv_reports
from source.logger import info


def main():
    info(' '.join(sys.argv))
    info()

    cnf, bcbio_structure = summary_script_proc_params(
        BCBioStructure.seq2c_name,
        BCBioStructure.cnv_summary_dir,
        extra_opts=[
           (['--controls', '-c'], dict(
                dest='controls',
                help='Optional control sample names for Seq2C. For multiple controls, separate them using :')
            ),
           (['--seq2c_opts'], dict(
                dest='seq2c_opts',
                help='Options for the final lr2gene.pl script.')
            ),
           (['--bed', '--capture', '--amplicons'], dict(
                dest='bed',
                help='BED file to run targetSeq and Seq2C analysis on.')
            ),
        ],
    )

    cnv_reports(cnf, bcbio_structure)


if __name__ == '__main__':
    main()

