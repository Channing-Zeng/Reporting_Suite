from optparse import OptionParser
from distutils import file_util
import sys
import os
from os.path import join, pardir, isfile, isdir, expanduser, dirname, abspath
from os import getcwd

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
    parser.add_option('--project-name', dest='project')
    parser.add_option('--email', dest='email')


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

    config_dirpath = abspath(join(bcbio_final_dir, pardir, 'config'))
    if not verify_dir(config_dirpath):
        critical('No config directory ' + dirname(config_dirpath) + ', check if you provided a correct path ' +
                 'to the bcbio directory. \nIt can be provided by the first argument for the script, nor the default one' +
                 'is the current working directory (' + bcbio_final_dir + ')')

    import socket
    hostname = socket.gethostname()
    info('hostname: ' + hostname)

    Defaults.sys_cnf = Defaults.sys_cnfs['us']
    if 'ukap' in hostname:
        Defaults.sys_cnf = Defaults.sys_cnfs['uk']
    elif 'local' in hostname:
        Defaults.sys_cnf = Defaults.sys_cnfs['local']
    elif any(name in hostname for name in ['rask', 'blue', 'chara', 'usbod']):
        Defaults.sys_cnf = Defaults.sys_cnfs['us']

    for file_basename, cnf_name in zip(['run', 'system'], ['run', 'sys']):
        provided_cnf_fpath = adjust_path(opt_dict.get(cnf_name + '_cnf'))
        project_cnf_fpath = adjust_path(join(config_dirpath, file_basename + '_info.yaml'))

        if provided_cnf_fpath:
            if not verify_file(provided_cnf_fpath):
                sys.exit(1)
            if provided_cnf_fpath != project_cnf_fpath:
                info('Using ' + provided_cnf_fpath + ', coping to ' + project_cnf_fpath)
                if isfile(project_cnf_fpath):
                    os.remove(project_cnf_fpath)
                file_util.copy_file(provided_cnf_fpath, project_cnf_fpath, preserve_times=False)

        else:  # No cnf in opts
            if isfile(project_cnf_fpath) and verify_file(project_cnf_fpath):
                provided_cnf_fpath = project_cnf_fpath
            else:
                provided_cnf_fpath = Defaults.__dict__[cnf_name + '_cnf']
                if provided_cnf_fpath != project_cnf_fpath:
                    info('Using ' + provided_cnf_fpath + ', coping to ' + project_cnf_fpath)
                    if isfile(project_cnf_fpath):
                        os.remove(project_cnf_fpath)
                    file_util.copy_file(provided_cnf_fpath, project_cnf_fpath, preserve_times=False)
                # critical('Usage: ' + __file__ + ' BCBIO_FINAL_DIR [--run-cnf YAML_FILE] [--sys-cnf YAML_FILE]')
            opt_dict[cnf_name + '_cnf'] = provided_cnf_fpath

    info()
    cnf = Config(opt_dict, opt_dict['sys_cnf'], opt_dict['run_cnf'])
    cnf.bcbio_final_dir = bcbio_final_dir

    if 'qsub_runner' in cnf:
        cnf.qsub_runner = expanduser(remove_quotes(cnf.qsub_runner))
        cnf.qsub_runner = adjust_path(join(__file__, pardir, pardir, cnf.qsub_runner))
    if not check_inputs(cnf, file_keys=['qsub_runner'], dir_keys=['bcbio_final_dir']):
        sys.exit(1)

    return cnf


def summary_script_proc_params(name, dir, description=None, extra_opts=None):
    description = description or 'This script generates project-level summaries based on per-sample ' + name + ' reports.'
    parser = OptionParser(description=description)
    add_post_bcbio_args(parser)

    parser.add_option('--log-dir', dest='log_dir')
    parser.add_option('--dir', dest='dir', default=dir)
    parser.add_option('--name', dest='name', default=name)
    for args, kwargs in extra_opts or []:
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



