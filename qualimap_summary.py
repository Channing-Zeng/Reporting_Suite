#!/usr/bin/env python

import sys
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from os.path import abspath, dirname, realpath, join, pardir
from site import addsitedir
source_dir = abspath(dirname(realpath(__file__)))
addsitedir(join(source_dir, 'ext_modules'))

from source.qualimap.summarize_qualimap import summary_reports
from source.bcbio_structure import BCBioStructure
from source.summary import summary_script_proc_params


def main():
    cnf, bcbio_structure = summary_script_proc_params(BCBioStructure.qualimap_name, BCBioStructure.qualimap_dir)

    summary_reports(cnf, bcbio_structure)


if __name__ == '__main__':
    main()

