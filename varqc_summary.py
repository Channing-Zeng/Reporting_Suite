#!/usr/bin/env python

from __future__ import print_function
import sys
from source.bcbio_structure import BCBioStructure, load_bcbio_cnf
from source.file_utils import verify_dir

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from os.path import abspath, dirname, realpath, pardir, join, basename
from site import addsitedir
source_dir = abspath(dirname(realpath(__file__)))
addsitedir(join(source_dir, 'ext_modules'))

from os.path import join, pardir
from optparse import OptionParser

from source.reporting import read_sample_names
from source.variants.summarize_qc import make_summary_reports
from source.config import Defaults, Config
from source.main import check_keys, check_inputs, set_up_dirs
from source.logger import info


def main():
    info(' '.join(sys.argv))
    info()

    description = 'This script generates project-level summaries based on per-sample targetcov reports.'

    parser = OptionParser(description=description)
    parser.add_option('-d', dest='bcbio_final_dir', help='Path to bcbio-nextgen final directory (default is pwd)')

    parser.add_option('-v', dest='verbose', action='store_true', help='Verbose')
    parser.add_option('-t', dest='threads', type='int', help='Number of threads for each process')
    parser.add_option('-w', dest='overwrite', action='store_true', help='Overwrite existing results')
    parser.add_option('-o', '--output_dir', dest='output_dir', metavar='DIR',
                      help='output directory (or directory name in case of bcbio final dir)')

    parser.add_option('--runner', dest='qsub_runner', help='Bash script that takes command line as the 1st argument. This script will be submitted to GRID. Default: ' + Defaults.qsub_runner)
    parser.add_option('--work-dir', dest='work_dir', metavar='DIR')
    parser.add_option('--sys-cnf', '--sys-info', '--sys-cfg', dest='sys_cnf', default=Defaults.sys_cnf,
                      help='System configuration yaml with paths to external tools and genome resources (see default one %s)' % Defaults.sys_cnf)
    parser.add_option('--run-cnf', '--run-info', '--run-cfg', dest='run_cnf', default=Defaults.run_cnf,
                      help='Run configuration yaml (see default one %s)' % Defaults.run_cnf)

    (opts, args) = parser.parse_args()
    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)

    cnf.name = cnf['name'] or 'varqc_summary'
    set_up_dirs(cnf)

    if not check_keys(cnf, ['bcbio_final_dir']):
        parser.print_help()
        sys.exit(1)

    cnf.bcbio_final_dir = verify_dir(cnf.bcbio_final_dir)
    if not cnf.bcbio_final_dir:
        sys.exit(1)

    info('BCBio "final" dir: ' + cnf.bcbio_final_dir + ' (set with -d)')

    if 'qsub_runner' in cnf:
        cnf.qsub_runner = join(cnf.sys_cnf, pardir, cnf.qsub_runner)

    if not check_inputs(cnf, file_keys=['qsub_runner'], dir_keys=['bcbio_final_dir']):
        sys.exit(1)

    info()
    info('*' * 70)

    load_bcbio_cnf(cnf)
    bcbio_structure = BCBioStructure(cnf, cnf.bcbio_final_dir, cnf.bcbio_cnf)
    make_summary_reports(cnf, bcbio_structure)


if __name__ == '__main__':
    main()















