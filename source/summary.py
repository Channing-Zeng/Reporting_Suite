#!/usr/bin/env python
from optparse import OptionParser
from os.path import join, pardir, isfile, isdir
from os import getcwd
import sys

from source.bcbio_structure import BCBioStructure, load_bcbio_cnf
from source.file_utils import verify_dir, safe_mkdir, adjust_path, verify_file, remove_quotes
from source.config import Defaults, Config
from source.main import check_keys, check_inputs, set_up_work_dir
from source.logger import info, critical


def add_post_bcbio_args(parser):
    parser.add_option('--sys-cnf', '--sys-info', '--sys-cfg', dest='sys_cnf', help='System configuration yaml with paths to external tools and genome resources (see default one %s)' % Defaults.sys_cnf)
    parser.add_option('--run-cnf', '--run-info', '--run-cfg', dest='run_cnf', help='Run configuration yaml (see default one %s)' % Defaults.run_cnf)
    parser.add_option('-v', dest='verbose', action='store_true', help='Verbose')
    parser.add_option('-t', dest='threads', type='int', help='Number of threads for each process')
    parser.add_option('-w', dest='overwrite', action='store_true', help='Overwrite existing results')
    parser.add_option('--reuse', dest='overwrite', action='store_false', help='Reuse intermediate files from previous run')
    parser.add_option('--runner', dest='qsub_runner', help='Bash script that takes command line as the 1st argument. This script will be submitted to GRID. Default: ' + Defaults.qsub_runner)


def process_post_bcbio_args(parser):
    (opts, args) = parser.parse_args()
    opt_dict = opts.__dict__

    dir_arg = args[0] if len(args) > 0 else getcwd()
    dir_arg = adjust_path(dir_arg)
    if not verify_dir(dir_arg):
        sys.exit(1)
    if isdir(join(dir_arg, 'final')):
        bcbio_final_dir = join(dir_arg, 'final')
    else:
        bcbio_final_dir = dir_arg

    info('BCBio "final" dir: ' + bcbio_final_dir)

    config_dirpath = join(bcbio_final_dir, pardir, 'config')
    for cnf_name in ['run', 'sys']:
        if not opt_dict.get(cnf_name + '_cnf'):
            cnf_fpath = adjust_path(join(config_dirpath, cnf_name + '_info.yaml'))
            if not isfile(cnf_fpath) or not verify_file(cnf_fpath):
                cnf_fpath = Defaults.__dict__[cnf_name + '_cnf']
                # critical('Usage: ' + __file__ + ' BCBIO_FINAL_DIR [--run-cnf YAML_FILE] [--sys-cnf YAML_FILE]')
            opt_dict[cnf_name + '_cnf'] = cnf_fpath

    cnf = Config(opt_dict, opt_dict['sys_cnf'], opt_dict['run_cnf'])
    cnf.bcbio_final_dir = bcbio_final_dir

    if 'qsub_runner' in cnf:
        cnf.qsub_runner = remove_quotes(cnf.qsub_runner)
        cnf.qsub_runner = adjust_path(join(cnf.sys_cnf, pardir, cnf.qsub_runner))
    if not check_inputs(cnf, file_keys=['qsub_runner'], dir_keys=['bcbio_final_dir']):
        sys.exit(1)

    return cnf


def summary_script_proc_params(name, dir, description=None, extra_opts=list()):
    description = description or 'This script generates project-level summaries based on per-sample ' + name + ' reports.'
    parser = OptionParser(description=description)
    add_post_bcbio_args(parser)

    parser.add_option('--log-dir', dest='log_dir')
    parser.add_option('--dir', dest='dir', default=dir)
    parser.add_option('--name', dest='name', default=name)
    for args, kwargs in extra_opts:
        parser.add_option(*args, **kwargs)

    cnf = process_post_bcbio_args(parser)

    load_bcbio_cnf(cnf)
    bcbio_structure = BCBioStructure(cnf, cnf.bcbio_final_dir, cnf.bcbio_cnf, cnf.name)

    cnf.work_dir = bcbio_structure.work_dir
    set_up_work_dir(cnf)
    if cnf.dir:
        cnf.output_dir = join(bcbio_structure.date_dirpath, cnf.dir)
        if not isdir(join(cnf.output_dir, pardir)):
            safe_mkdir(join(cnf.output_dir, pardir))
        if not isdir(cnf.output_dir):
            safe_mkdir(cnf.output_dir)

    info('*' * 70)
    info()

    return cnf, bcbio_structure



