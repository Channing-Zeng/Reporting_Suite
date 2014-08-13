#!/usr/bin/env python

import sys
from os.path import join, pardir
from optparse import OptionParser

from source.bcbio_structure import BCBioStructure, load_bcbio_cnf
from source.file_utils import verify_dir
from source.config import Defaults, Config
from source.main import check_keys, check_inputs, set_up_dirs, set_up_work_dir, set_up_log
from source.logger import info, critical


def summary_script_proc_params(name, description=None, extra_opts=list()):
    info(' '.join(sys.argv))
    info()

    description = description or 'This script generates project-level summaries based on per-sample ' + name + ' reports.'

    parser = OptionParser(description=description)
    parser.add_option('-d', dest='bcbio_final_dir', help='Path to bcbio-nextgen final directory (default is pwd)')
    parser.add_option('-v', dest='verbose', action='store_true', help='Verbose')
    parser.add_option('-t', dest='threads', type='int', help='Number of threads for each process')
    parser.add_option('-w', dest='overwrite', action='store_true', help='Overwrite existing results')
    parser.add_option('--reuse', dest='overwrite', help='Reuse intermediate files from previous run', action='store_false')

    parser.add_option('--runner', dest='qsub_runner', help='Bash script that takes command line as the 1st argument. This script will be submitted to GRID. Default: ' + Defaults.qsub_runner)
    parser.add_option('--sys-cnf', '--sys-info', '--sys-cfg', dest='sys_cnf', default=Defaults.sys_cnf, help='System configuration yaml with paths to external tools and genome resources (see default one %s)' % Defaults.sys_cnf)
    parser.add_option('--run-cnf', '--run-info', '--run-cfg', dest='run_cnf', default=Defaults.run_cnf, help='Run configuration yaml (see default one %s)' % Defaults.run_cnf)

    for args, kwargs in extra_opts:
        parser.add_option(*args, **kwargs)

    (opts, args) = parser.parse_args()
    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)
    if not opts.bcbio_final_dir and len(args) > 0:
        cnf.bcbio_final_dir = args[0]
    else:
        critical('Usage: ./' + __file__ + ' <final_dir>')

    _check_args(parser, cnf)

    load_bcbio_cnf(cnf)

    bcbio_structure = BCBioStructure(cnf, cnf.bcbio_final_dir, cnf.bcbio_cnf)
    cnf.work_dir = bcbio_structure.work_dir
    cnf.name = name

    set_up_work_dir(cnf)
    set_up_log(cnf)

    info()
    info('*' * 70)

    return cnf, bcbio_structure


def _check_args(parser, cnf):
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


