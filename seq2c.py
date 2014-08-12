#!/usr/bin/env python

from __future__ import print_function
import sys

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from os.path import abspath, dirname, realpath, join
from site import addsitedir
source_dir = abspath(dirname(realpath(__file__)))
addsitedir(join(source_dir, 'ext_modules'))

from source.logger import info
from source.bcbio_structure import BCBioStructure
from source.summary import process_cnf
from source.targetcov.copy_number import cnv_reports


def main():
    info(' '.join(sys.argv))
    info()

    cnf, bcbio_structure = process_cnf(BCBioStructure.targetseq_summary_dir)

    cnv_report_fpath = cnv_reports(cnf, bcbio_structure)

    info()
    info('*' * 70)

    if cnv_report_fpath:
        info('Gene CNV:')
        info('  ' + cnv_report_fpath)


if __name__ == '__main__':
    main()

