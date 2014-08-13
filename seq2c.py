#!/usr/bin/env python

from __future__ import print_function
from genericpath import isdir
import os
import sys
from os import listdir

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from os.path import abspath, dirname, realpath, join, islink
from site import addsitedir
source_dir = abspath(dirname(realpath(__file__)))
addsitedir(join(source_dir, 'ext_modules'))

from source.bcbio_structure import BCBioStructure
from source.summary import summary_script_proc_params
from source.targetcov.copy_number import cnv_reports
from source.file_utils import safe_mkdir


def main():
    cnf, bcbio_structure = summary_script_proc_params(BCBioStructure.seq2c_name)

    cnf.output_dir = join(bcbio_structure.date_dirpath, BCBioStructure.cnv_summary_dir)
    safe_mkdir(cnf.output_dir)

    cnv_reports(cnf, bcbio_structure)


if __name__ == '__main__':
    main()

