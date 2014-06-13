#!/usr/bin/env python
import sys
import os
from os.path import isdir, dirname, join, realpath, expanduser, basename, abspath
from optparse import OptionParser
from shutil import rmtree

from source import logger
from source.logger import info, err, critical
from source.bcbio_utils import which, file_exists
from source.utils import verify_file, safe_mkdir, verify_dir, verify_module
from source.ngscat.bed_file import verify_bam, verify_bed


if verify_module('yaml'):
    from yaml import load
    try:
        from yaml import CDumper as Dumper, CLoader as Loader
    except ImportError:
        from yaml import Dumper, Loader
else:
    critical('Cannot import module yaml.')


default_sys_config_path = 'system_info_Waltham.yaml'
default_run_config_path = 'run_info.yaml'


def read_opts_and_cnfs(extra_opts, key_for_sample_name):
    options = [
        (['-o', '--output_dir'], 'DIR', {
             'dest': 'output_dir',
             'help': 'output directory (or directory name in case of bcbio final dir)'}),

        (['-s', '--sample'], 'NAME', {
             'dest': 'name',
             'help': 'sample name (default is part of name of the first parameter prior to the first - or .'}),

        # (['--bcbio-final-dir'], 'DIR', {
        #     'dest': 'bcbio_final_dir',
        #     'help': 'bcbio-nextgen "final" directory (assumes bcbio project structure)'}),
        #
        # (['-ss', '--samples'], 'SAMPLES.txt', {
        #     'dest': 'samples',
        #     'help': 'list of samples (assumes bcbio project structure)'}),

        (['-t', '--nt', '--threads'], 'N', {
             'dest': 'threads',
             'type': 'int',
             'help': 'number of threads'}),

        (['-w', '--overwrite'], 'False', {
             'dest': 'overwrite',
             'help': 'do not reuse intermediate files from previous run',
             'action': 'store_true',
             'default': False}),

        (['--sys-cnf'], 'SYS_CNF.yaml', {
             'dest': 'sys_cnf',
             'help': 'system configuration yaml with paths to external tools and genome resources '
                     '(see default one %s)' % default_sys_config_path,
             'default': default_sys_config_path}),

        (['--run-cnf'], 'RUN_CNF.yaml', {
             'dest': 'run_cnf',
             'help': 'run configuration yaml (see default one %s)' % default_run_config_path,
             'default': default_run_config_path}),
    ] + extra_opts

    # opts_line = ' '.join(args[0] + ' ' + example for args, example, _ in options)
    # format_params = {
    #     'script': __file__,
    #     'opts_line': opts_line,
    #     'sys_cnf': default_sys_config_path,
    #     'run_cnf': default_run_config_path}
    parser = OptionParser(
        # usage=(
        #     'python {script} [{sys_cnf}] [{run_cnf}] {opts_line}\n'
        #     'or\n'
        #     'python {script} [{sys_cnf}] {run_cnf}\n'
        #     'or\n'
        #     'python {script} {sys_cnf}'.format(**format_params))
    )
    for args, _, kwargs in options:
        parser.add_option(*args, **kwargs)

    (opt_obj, args) = parser.parse_args()
    opts = dict((k, v) for k, v in opt_obj.__dict__.items() if v is not None)

    if opts['overwrite']:
        opts['reuse_intermediate'] = False
        del opts['overwrite']

    cnf = load_configs(opt_obj.sys_cnf, opt_obj.run_cnf)

    cnf.update(opts)

    # if 'bcbio_final_dir' in opts:
    #     if 'samples' not in opts:
    #         critical('Error: you specified bcbio_final_dir, but not samples. '
    #                  'You have to specify a file with a list of sample names with '
    #                  'the --samples option as well.')
    #     opts['bcbio_final_dir'] = abspath(expanduser(opts['bcbio_final_dir']))
    #     if not verify_dir(opts['bcbio_final_dir']):
    #         sys.exit(1)
    #
    # if 'samples' in opts:
    #     if 'samples' not in opts:
    #         critical('Error: you specified bcbio_final_dir, but not samples. '
    #                  'You have to specify a file with a list of sample names with '
    #                  'the --samples option as well.')
    #     opts['samples'] = abspath(expanduser(opts['samples']))
    #     if not verify_file(opts['samples']):
    #         sys.exit(1)

    assert key_for_sample_name
    cnf['name'] = \
        cnf.get('name') or \
        basename(dirname(cnf[key_for_sample_name])) or \
        basename(cnf[key_for_sample_name])

    set_up_dirs(cnf)

    return cnf


def check_inputs(cnf, required_keys, file_keys):
    to_exit = False

    def _verify_input(_key):
        if not verify_file(cnf[_key], _key):
            return False
        if 'bam' in _key and not verify_bam(cnf[_key]):
            return False
        if 'bed' in _key and not verify_bed(cnf[_key]):
            return False
        return True

    for key in required_keys:
        if not key or key not in cnf:
            to_exit = True
            err('Error: ' + key + ' must be provided in options or '
                'in ' + cnf['run_config_fpath'] + '.')

    for key in file_keys:
        if key and key in cnf:
            if not _verify_input(key):
                to_exit = True

    if to_exit:
        sys.exit(1)


def _fill_config_from_defaults(cnf, defaults):
    for key in defaults:
        if key in cnf:
            if isinstance(cnf[key], dict) and isinstance(defaults[key], dict):
                _fill_config_from_defaults(cnf[key], defaults[key])
        else:
            cnf[key] = defaults[key]
    return cnf


def load_configs(sys_cnf, run_cnf):
    defaults_config_path = join(dirname(dirname(realpath(__file__))), join('source', 'defaults.yaml'))

    sys_config_path = abspath(sys_cnf)
    run_config_path = abspath(run_cnf)
    info('Using ' + sys_config_path + ' as a system configuration file.')
    info('Using ' + run_config_path + ' as a run configuration file.')
    info()

    to_exit = False
    for f in [sys_config_path, run_config_path]:
        if not verify_file(f):
            to_exit = True
    if to_exit:
        sys.exit(1)

    for f in [sys_config_path, run_config_path]:
        if not f.endswith('.yaml'):
            err(f + ' does not end with .yaml, maybe incorrect parameter?\n')
            to_exit = True
    if to_exit:
        sys.exit(1)

    dft_cnf = load(open(defaults_config_path), Loader=Loader)
    sys_cnf = load(open(sys_config_path), Loader=Loader)
    run_cnf = load(open(run_config_path), Loader=Loader)

    cnf = dict(run_cnf.items() + sys_cnf.items())
    _fill_config_from_defaults(cnf, dft_cnf)

    cnf['system_config_path'] = sys_config_path
    cnf['run_config_path'] = run_config_path

    info('Loaded system config ' + sys_config_path)
    info('Loaded run config ' + run_config_path)
    return cnf


def input_fpaths_from_cnf(cnf, required_inputs, optional_inputs):
    input_fpaths = []
    for key in required_inputs:
        input_fpaths.append(cnf[key])
    for key in optional_inputs:
        if key in cnf:
            input_fpaths.append(cnf[key])
    return input_fpaths


def check_system_resources(cnf, required=list(), optional=list()):
    to_exit = False

    for program in required:
        if not which(program):
            resources = cnf.get('resources', None)
            if not resources:
                critical('No "resources" section in system config.')

            data = resources.get(program)
            if data is None:
                err(program + ' is required. Specify path in system config or in your environment.')
                to_exit = True
            else:
                data['path'] = expanduser(data['path'])
                if not isdir(data['path']) and not file_exists(data['path']):
                    err(data['path'] + ' does not exist.')
                    to_exit = True

    for program in optional:
        resources = cnf.get('resources', None)
        if not resources:
            break

        data = resources.get(program)
        if data is None:
            continue
        else:
            data['path'] = expanduser(data['path'])
            if not isdir(data['path']) and not file_exists(data['path']):
                err(data['path'] + ' does not exist.')
                to_exit = True

    if to_exit:
        exit()


def load_genome_resources(cnf, required=list(), optional=list()):
    if 'genome' not in cnf:
        critical('"genome" is not specified in run config.')
    if 'genomes' not in cnf:
        critical('"genomes" section is not specified in system config.')
    genome_name = cnf['genome']
    if genome_name not in cnf['genomes']:
        critical(genome_name + ' is not in "genomes section" of system config.')
    genome_cnf = cnf['genomes'][genome_name].copy()

    to_exit = False

    if 'seq' in required:
        required.remove('seq')

        if 'seq' not in genome_cnf:
            err('Please, provide path to the reference file (seq).')
            to_exit = True

        genome_cnf['seq'] = abspath(expanduser(genome_cnf['seq']))
        if not verify_file(genome_cnf['seq'], 'Reference seq'):
            to_exit = True

    if 'snpeff' in required:
        required.remove('snpeff')
        if 'snpeff' in genome_cnf:
            genome_cnf['snpeff'] = abspath(expanduser(genome_cnf['snpeff']))
            if not verify_dir(genome_cnf['snpeff'], 'snpeff'):
                to_exit = True

    for key in required:  # 'dbsnp', 'cosmic', 'dbsnfp', '1000genomes':
        if key not in genome_cnf:
            err('Please, provide path to ' + key + ' in system config genome section.')
            to_exit = True
        else:
            genome_cnf[key] = abspath(expanduser(genome_cnf[key]))
            if not verify_file(genome_cnf[key], key):
                to_exit = True

    for key in optional:
        if key not in genome_cnf:
            continue
        else:
            genome_cnf[key] = abspath(expanduser(genome_cnf[key]))
            if not verify_file(genome_cnf[key], key):
                to_exit = True

    if to_exit:
        sys.exit(1)

    cnf['genome'] = genome_cnf
    del cnf['genomes']
    genome_cnf['name'] = genome_name

    info('Loaded resources for ' + genome_cnf['name'])


def set_up_dirs(cnf):
    cnf['output_dir'] = cnf.get('output_dir') or os.getcwd()
    cnf['output_dir'] = realpath(expanduser(cnf['output_dir']))

    safe_mkdir(cnf['output_dir'], 'output_dir')
    info('Saving into ' + cnf['output_dir'])

    work_dir_name = 'work_' + cnf['name']
    cnf['work_dir'] = join(cnf['output_dir'], work_dir_name)
    if not cnf.get('reuse_intermediate') and isdir(cnf['work_dir']):
        rmtree(cnf['work_dir'])

    safe_mkdir(cnf['work_dir'], 'working directory')

    cnf['log'] = join(cnf['work_dir'], cnf['name'] + '_log.txt')
    logger.log_fname = cnf['log']