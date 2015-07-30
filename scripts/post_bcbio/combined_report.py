#!/usr/bin/env python

import __check_python_version

import sys
from source.bcbio_structure import summary_script_proc_params
from source.project_level_report import make_postproc_project_level_report
from source.logger import info


def main():
    info(' '.join(sys.argv))
    info()

    cnf, bcbio_structure = summary_script_proc_params('project_level_report')

    make_postproc_project_level_report(cnf, bcbio_structure)


if __name__ == '__main__':
    main()