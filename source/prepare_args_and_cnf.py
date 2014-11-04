from optparse import OptionParser
from distutils import file_util
import sys
import os
from os.path import join, pardir, isfile, isdir, expanduser, dirname, abspath
from os import getcwd

from source.bcbio_structure import BCBioStructure, load_bcbio_cnf
from source.file_utils import verify_dir, safe_mkdir, adjust_path, verify_file, adjust_system_path
from source.config import Defaults, Config
from source.main import check_keys, check_inputs, set_up_work_dir, load_genome_resources
from source.logger import info, critical


def add_post_bcbio_args(parser):
    parser.add_option('--sys-cnf', '--sys-info', '--sys-cfg', dest='sys_cnf', help='System configuration yaml with paths to external tools and genome resources (see default one %s)' % Defaults.sys_cnf)
    parser.add_option('--run-cnf', '--run-info', '--run-cfg', dest='run_cnf', help='Run configuration yaml (see default one %s)' % Defaults.run_cnf)
    parser.add_option('-v', dest='verbose', action='store_true', help='Verbose')
    parser.add_option('-t', dest='threads', type='int', help='Number of threads for each process')
    parser.add_option('-w', dest='overwrite', action='store_true', help='Overwrite existing results')
    parser.add_option('--reuse', dest='overwrite', action='store_false', help='Reuse intermediate files from previous run')
    parser.add_option('--runner', dest='qsub_runner', help='Bash script that takes command line as the 1st argument. This script will be submitted to GRID. Default: ' + Defaults.qsub_runner)
    parser.add_option('--project-name', '--project', dest='project_name')
    parser.add_option('--email', dest='email')


def process_post_bcbio_args(parser):
    (opts, args) = parser.parse_args()

    dir_arg = args[0] if len(args) > 0 else getcwd()
    dir_arg = adjust_path(dir_arg)
    if not verify_dir(dir_arg):
        sys.exit(1)

    bcbio_project_dirpath, final_dirpath, config_dirpath = _set_bcbio_dirpath(dir_arg)

    _set_sys_config(opts)

    _set_run_cnf(config_dirpath, opts)

    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)

    if 'qsub_runner' in cnf:
        cnf.qsub_runner = adjust_system_path(cnf.qsub_runner)
    if not check_inputs(cnf, file_keys=['qsub_runner']):
        sys.exit(1)

    load_genome_resources(cnf, required=['seq'])

    bcbio_cnf = load_bcbio_cnf(cnf, config_dirpath)

    return cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath


def summary_script_proc_params(name, dir, description=None, extra_opts=None):
    description = description or 'This script generates project-level summaries based on per-sample ' + name + ' reports.'
    parser = OptionParser(description=description)
    add_post_bcbio_args(parser)

    parser.add_option('--log-dir', dest='log_dir')
    parser.add_option('--dir', dest='dir', default=dir)  # to distinguish VarQC_summary and VarQC_after_summary
    parser.add_option('--name', dest='name', default=name)  # proc name
    for args, kwargs in extra_opts or []:
        parser.add_option(*args, **kwargs)

    cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath = process_post_bcbio_args(parser)

    bcbio_structure = BCBioStructure(cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath, cnf.name)

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


def _set_bcbio_dirpath(dir_arg):
    if isdir(join(dir_arg, 'config')):
        bcbio_project_dirpath = dir_arg
        final_dirpath = None
        config_dirpath = join(bcbio_project_dirpath, 'config')

    elif isdir(abspath(join(dir_arg, pardir, 'config'))):
        bcbio_project_dirpath = abspath(join(dir_arg, pardir))
        final_dirpath = dir_arg
        config_dirpath = join(bcbio_project_dirpath, 'config')

    else:
        critical(
            'No config directory ' + join(dir_arg, 'config') + ' or ' + abspath(join(dir_arg, pardir, 'config')) +
            ', check if you provided a correct path to the bcbio directory.'
            '\nIt can be provided by the first argument for the script, or by changing to it.')

    info('BCBio project directory: ' + bcbio_project_dirpath)
    if final_dirpath: info('final directory: ' + final_dirpath)
    info('Config directory: ' + config_dirpath)
    return bcbio_project_dirpath, final_dirpath, config_dirpath


def _set_sys_config(opts):
    opts.sys_cnf = adjust_path(opts.sys_cnf)
    if opts.sys_cnf:
        if not verify_file(opts.sys_cnf):
            sys.exit(1)
    else:
        import socket
        hostname = socket.gethostname()
        info('hostname: ' + hostname)

        opts.sys_cnf = Defaults.sys_cnfs['us']
        if 'ukap' in hostname:
            opts.sys_cnf = Defaults.sys_cnfs['uk']
        elif 'local' in hostname:
            opts.sys_cnf = Defaults.sys_cnfs['local']
        elif any(name in hostname for name in ['rask', 'blue', 'chara', 'usbod']):
            opts.sys_cnf = Defaults.sys_cnfs['us']
    info('Using ' + opts.sys_cnf)


def _set_run_cnf(config_dirpath, opts):
    provided_run_cnf_fpath = adjust_path(opts.run_cnf)
    project_run_cnf_fpath = adjust_path(join(config_dirpath, 'run_info.yaml'))
    if provided_run_cnf_fpath:
        if not verify_file(provided_run_cnf_fpath):
            sys.exit(1)
        if provided_run_cnf_fpath != project_run_cnf_fpath:
            info('Using ' + provided_run_cnf_fpath + ', copying to ' + project_run_cnf_fpath)
            if isfile(project_run_cnf_fpath):
                try:
                    os.remove(project_run_cnf_fpath)
                except OSError:
                    pass
            if not isfile(project_run_cnf_fpath):
                file_util.copy_file(provided_run_cnf_fpath, project_run_cnf_fpath, preserve_times=False)

    else:  # No configs provided in command line options
        if isfile(project_run_cnf_fpath) and verify_file(project_run_cnf_fpath):
            provided_run_cnf_fpath = project_run_cnf_fpath
        else:
            provided_run_cnf_fpath = Defaults.run_cnf
            if provided_run_cnf_fpath != project_run_cnf_fpath:
                info('Using ' + provided_run_cnf_fpath + ', copying to ' + project_run_cnf_fpath)
                if isfile(project_run_cnf_fpath):
                    try:
                        os.remove(project_run_cnf_fpath)
                    except OSError:
                        pass
                if not isfile(project_run_cnf_fpath):
                    file_util.copy_file(provided_run_cnf_fpath, project_run_cnf_fpath, preserve_times=False)
                    # critical('Usage: ' + __file__ + ' BCBIO_FINAL_DIR [--run-cnf YAML_FILE] [--sys-cnf YAML_FILE]')
    opts.run_cnf = project_run_cnf_fpath
    info()


