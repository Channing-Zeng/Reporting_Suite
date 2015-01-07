from optparse import OptionParser
from distutils import file_util
import sys
import os
from os.path import join, pardir, isfile, isdir, expanduser, dirname, abspath, basename
from os import getcwd

from source.bcbio_structure import BCBioStructure, load_bcbio_cnf
from source.file_utils import verify_dir, safe_mkdir, adjust_path, verify_file, adjust_system_path
from source.config import defaults, Config
from source.main import check_keys, check_inputs, set_up_work_dir, check_genome_resources
from source.logger import info, critical, warn


def add_post_bcbio_args(parser):
    parser.add_option('--sys-cnf', '--sys-info', '--sys-cfg', dest='sys_cnf', help='System configuration yaml with paths to external tools and genome resources (see default one %s)' % defaults['sys_cnf'])
    parser.add_option('--run-cnf', '--run-info', '--run-cfg', dest='run_cnf', help='Run configuration yaml (see default one %s)' % defaults['run_cnf'])
    parser.add_option('-v', dest='verbose', action='store_true', help='Verbose')
    parser.add_option('-t', dest='threads', type='int', help='Number of threads for each process')
    parser.add_option('--reuse', dest='reuse_intermediate', action='store_true', help='Reuse intermediate non-empty files in the work dir from previous run')
    parser.add_option('--runner', dest='qsub_runner', help='Bash script that takes command line as the 1st argument. This script will be submitted to GRID. Default: ' + defaults['qsub_runner'])
    parser.add_option('--project-name', '--project', dest='project_name')
    parser.add_option('--email', dest='email')
    parser.add_option('--bed', dest='bed', help='BED file to run targetSeq and Seq2C analysis on.')
    parser.add_option('--exons', '--exome', dest='exons', help='Exons BED file to make targetSeq exon/amplicon regions reports.')


def process_post_bcbio_args(parser):
    (opts, args) = parser.parse_args()

    dir_arg = args[0] if len(args) > 0 else getcwd()
    dir_arg = adjust_path(dir_arg)
    if not verify_dir(dir_arg):
        sys.exit(1)

    bcbio_project_dirpath, final_dirpath, config_dirpath = _set_bcbio_dirpath(dir_arg)

    _set_sys_config(config_dirpath, opts)

    _set_run_config(config_dirpath, opts)

    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)

    if 'qsub_runner' in cnf:
        cnf.qsub_runner = adjust_system_path(cnf.qsub_runner)
    if not check_inputs(cnf, file_keys=['qsub_runner']):
        sys.exit(1)

    bcbio_cnf = load_bcbio_cnf(cnf, config_dirpath)

    check_genome_resources(cnf)

    return cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath


def summary_script_proc_params(name, dir_name=None, description=None, extra_opts=None):
    description = description or 'This script generates project-level summaries based on per-sample ' + name + ' reports.'
    parser = OptionParser(description=description)
    add_post_bcbio_args(parser)

    parser.add_option('--log-dir', dest='log_dir')
    parser.add_option('--dir', dest='dir_name', default=dir_name, help='Optional - to distinguish VarQC_summary and VarQC_after_summary')
    parser.add_option('--name', dest='name', default=name, help='Procedure name')
    for args, kwargs in extra_opts or []:
        parser.add_option(*args, **kwargs)

    cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath = process_post_bcbio_args(parser)

    bcbio_structure = BCBioStructure(cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath, cnf.name)

    cnf.output_dir = join(bcbio_structure.date_dirpath, cnf.dir_name) if cnf.dir_name else None
    cnf.work_dir = bcbio_structure.work_dir
    set_up_work_dir(cnf)

    info('*' * 70)
    info()

    return cnf, bcbio_structure


def _set_bcbio_dirpath(dir_arg):
    final_dirpath, bcbio_project_dirpath, config_dirpath = None, None, None

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


def _set_sys_config(config_dirpath, opts):
    provided_cnf_fpath = adjust_path(opts.sys_cnf)
    # provided in commandline?
    if provided_cnf_fpath:
        if not verify_file(provided_cnf_fpath):
            sys.exit(1)
        # alright, in commandline.
        opts.sys_cnf = provided_cnf_fpath

    else:
        # or probably in config dir?
        fpaths_in_config = [
            abspath(join(config_dirpath, fname))
            for fname in os.listdir(config_dirpath)
            if fname.startswith('system_info') and fname.endswith('.yaml')]

        if len(fpaths_in_config) > 1:
            critical('More than one YAML file containing run_info in name found in the config '
                     'directory ' + config_dirpath + ': ' + ' '.join(fpaths_in_config))

        elif len(fpaths_in_config) == 1:
            opts.sys_cnf = fpaths_in_config[0]
            if not verify_file(opts.sys_cnf):
                sys.exit(1)
            # alright, in config dir.

        else:
            # detect system and use default system config
            import socket
            hostname = socket.gethostname()
            info('hostname: ' + hostname)

            opts.sys_cnf = defaults['sys_cnfs']['us']
            if 'ukap' in hostname:
                opts.sys_cnf = defaults['sys_cnfs']['uk']
            elif 'cniclhpc' in hostname:
                opts.sys_cnf = defaults['sys_cnfs']['china']
            elif 'local' in hostname or 'Home' in hostname:
                opts.sys_cnf = defaults['sys_cnfs']['local']
            elif any(name in hostname for name in ['rask', 'blue', 'chara', 'usbod']):
                opts.sys_cnf = defaults['sys_cnfs']['us']

    info('Using ' + opts.sys_cnf)


def _set_run_config(config_dirpath, opts):
    provided_cnf_fpath = adjust_path(opts.run_cnf)
    # provided in commandline?
    if provided_cnf_fpath:
        if not verify_file(provided_cnf_fpath):
            sys.exit(1)

        # alright, in commandline. copying over to config dir.
        opts.run_cnf = provided_cnf_fpath
        project_run_cnf_fpath = adjust_path(join(config_dirpath, basename(provided_cnf_fpath)))
        if provided_cnf_fpath != project_run_cnf_fpath:
            info('Using ' + provided_cnf_fpath + ', copying to ' + project_run_cnf_fpath)
            if isfile(project_run_cnf_fpath):
                try:
                    os.remove(project_run_cnf_fpath)
                except OSError:
                    pass

            if not isfile(project_run_cnf_fpath):
                run_info_fpaths_in_config = [
                    abspath(join(config_dirpath, fname))
                    for fname in os.listdir(config_dirpath)
                    if fname.startswith('run_info') and fname.endswith('.yaml')]

                if len(run_info_fpaths_in_config) > 0:
                    warn('Warning: there are run_info files in config directory ' + config_dirpath + '. '
                         'Provided config will be copied there and can cause ambigity in future.')

                file_util.copy_file(provided_cnf_fpath, project_run_cnf_fpath, preserve_times=False)

    else:  # no configs provided in command line options
        run_info_fpaths_in_config = [
            abspath(join(config_dirpath, fname))
            for fname in os.listdir(config_dirpath)
            if fname.startswith('run_info') and fname.endswith('.yaml')]

        if len(run_info_fpaths_in_config) > 1:
            critical('More than one YAML file containing run_info in name found in the config '
                     'directory ' + config_dirpath + ': ' + ' '.join(run_info_fpaths_in_config))

        elif len(run_info_fpaths_in_config) == 1:
            opts.run_cnf = run_info_fpaths_in_config[0]
            if not verify_file(opts.run_cnf):
                sys.exit(1)
            # alright, in config dir.

        elif len(run_info_fpaths_in_config) == 0:
            info('No YAMLs containing run_info in name found in the config directory ' +
                 config_dirpath + ', using the default one.')

            # using default one.
            opts.run_cnf = defaults['run_cnf']
            project_run_cnf_fpath = adjust_path(join(config_dirpath, basename(opts.run_cnf)))
            info('Using ' + opts.run_cnf + ', copying to ' + project_run_cnf_fpath)
            if isfile(project_run_cnf_fpath):
                try:
                    os.remove(project_run_cnf_fpath)
                except OSError:
                    pass
            if not isfile(project_run_cnf_fpath):
                file_util.copy_file(opts.run_cnf, project_run_cnf_fpath, preserve_times=False)

    info('Using ' + opts.run_cnf)


