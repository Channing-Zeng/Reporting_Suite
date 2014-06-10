#!/usr/bin/env python
import sys
import os
from os.path import isdir, dirname, join, realpath, expanduser, basename, abspath
from optparse import OptionParser
from source import logger

from source.utils import err, critical, verify_file, \
    info, safe_mkdir, verify_dir, verify_module

if verify_module('yaml'):
    from yaml import load
    try:
        from yaml import CDumper as Dumper, CLoader as Loader
    except ImportError:
        from yaml import Dumper, Loader
else:
    critical('Cannot import module yaml.')

from source.bcbio_utils import which, file_exists


def read_opts_and_cnfs(extra_opts, required_keys, optional_keys, fpaths_keys=None):
    run_config_name = 'run_info.yaml'
    system_config_name = 'system_info_Waltham.yaml'

    basic_opts = [
        (['-o', '--output_dir'], 'DIR', {
         'dest': 'output_dir',
         'help': 'output directory'}),

        (['--sample'], 'NAME', {
         'dest': 'name',
         'help': 'sample name (default is first input file directory name'}),

        (['-t', '--nt', '--threads'], 'N', {
         'dest': 'threads',
         'help': 'number of threads'}),

        (['--overwrite'], 'False', {
         'dest': 'overwrite',
         'help': 'do not reuse intermediate files from previous run',
         'action': 'store_true',
         'default': False}),
    ]

    options = basic_opts + extra_opts
    opts_line = ' '.join(args[0] + ' ' + example for args, example, _ in options)
    filename = __file__
    parser = OptionParser(
        usage=(
            'python {filename} [{system_config_name}] [{run_config_name}] {opts_line}\n'
            'or\n'
            'python {filename} [{system_config_name}] {run_config_name}\n'
            'or\n'
            'python {filename} {system_config_name}'.format(**locals())
        ))
    for args, _, kwargs in options:
        parser.add_option(*args, **kwargs)

    (opt_obj, args) = parser.parse_args()
    opts = dict((k, v) for k, v in opt_obj.__dict__.items() if v is not None)

    if opts['overwrite']:
        opts['reuse_intermediate'] = False
        del opts['overwrite']

    cnf = _load_configs(run_config_name, system_config_name, args)

    cnf.update(opts)

    assert required_keys
    cnf['name'] = cnf.get('name') or basename(dirname(cnf[required_keys[0]]))

    set_up_dirs(cnf)

    check_inputs(cnf, required_keys, optional_keys, fpaths_keys)

    return cnf


def check_inputs(cnf, required_keys, optional_keys, fpaths_keys=None):
    if fpaths_keys is None:
        fpaths_keys = required_keys + optional_keys

    to_exit = False

    for key in required_keys:
        if not key or key not in cnf:
            to_exit = True
            err('Error: ' + key + ' must be provided in options or ' + cnf['run_config_fpath'] + '.')
        if key in fpaths_keys:
            if not verify_file(cnf[key], key):
                to_exit = True

    for key in optional_keys:
        if key and key in fpaths_keys:
            if not verify_file(cnf[key], key):
                to_exit = True

    if to_exit:
        sys.exit(1)

    info('Input:')
    for key in required_keys + optional_keys:
        info('  ' + key + ': ' + cnf[key])


def _load_configs(run_config_name, system_config_name, args):
    system_config_path = join(dirname(dirname(realpath(__file__))), system_config_name)
    run_config_path = join(dirname(dirname(realpath(__file__))), run_config_name)
    if len(args) < 1:
        err('Notice: using ' + run_config_name + ' as a default run configutation file.\n\n')
    else:
        run_config_path = args[0]

    if len(args) < 2:
        err('Notice: using system_info_rask.yaml as a default tools configutation file.\n\n')
    else:
        system_config_path = args[0]
        run_config_path = args[1]

    if not os.path.isfile(system_config_path):
        critical(system_config_path + ' does not exist or is a directory.\n')
    if not os.path.isfile(run_config_path):
        critical(run_config_path + ' does not exist or is a directory.\n')

    to_exit = False
    if not system_config_path.endswith('.yaml'):
        err(system_config_path + ' does not end with .yaml, maybe incorrect parameter?\n')
        to_exit = True
    if not run_config_path.endswith('.yaml'):
        err(run_config_path + ' does not end with .yaml, maybe incorrect parameter?\n')
        to_exit = True

    if to_exit:
        sys.exit(1)

    sys_cnf = load(open(system_config_path), Loader=Loader)
    run_cnf = load(open(run_config_path), Loader=Loader)
    info('Loaded system config ' + system_config_path)
    info('Loaded run config ' + run_config_path)

    cnf = dict(run_cnf.items() + sys_cnf.items())

    cnf['system_config_path'] = system_config_path
    cnf['run_config_path'] = run_config_path

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
    cnf['work_dir'] = join(cnf['output_dir'], 'work')
    safe_mkdir(cnf['work_dir'], 'working directory')

    cnf['log'] = join(cnf['work_dir'], cnf['name'] + '_log.txt')
    logger.log_fname = cnf['log']