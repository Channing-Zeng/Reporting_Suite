#!/usr/bin/env python
import sys
import os
from os.path import isdir, dirname, join, realpath, expanduser
from optparse import OptionParser

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


def read_opts_and_cnfs(extra_opts, required_inputs):
    run_config_name = 'run_info.yaml'
    system_config_name = 'system_info_Waltham.yaml'

    opts = ' '.join(args[0] + ' ' + example for args, example, _ in extra_opts)

    filename = __file__
    parser = OptionParser(
        usage=(
            'python {filename} [{system_config_name}] [{run_config_name}] {opts} [--output_dir dir]\n'
            'or\n'
            'python {filename} [{system_config_name}] {run_config_name}'.format(**locals())
        )
    )

    for args, _, kwargs in extra_opts:
        parser.add_option(*args, **kwargs)
    parser.add_option('-o', '--output_dir', dest='output_dir', metavar='DIR')
    parser.add_option('-t', '--nt', dest='threads', help='number of threads')
    parser.add_option('-w', dest='rewrite', action='store_true', default=False,
                      help='do not reuse intermediate files from previous run')
    parser.add_option('--np', '--no-parallel', dest='no_parallel', action='store_true', default=False)
    (options, args) = parser.parse_args()
    options = options.__dict__

    cnf = _load_configs(run_config_name, system_config_name, args)

    if options.get('output_dir'):
        cnf['output_dir'] = options.get('output_dir')
    if cnf.get('output_dir'):
        cnf['output_dir'] = expanduser(cnf['output_dir'])
        _set_up_dirs(cnf)

    if (len([k for k in required_inputs if k in options]) < len(required_inputs) and
        not cnf.get('details')):
        critical(
            'Error: provide input files with command line options '
            'or by specifying in the "details" section in run config.')

    cnf['filter_reject'] = cnf.get('filter_reject', False)
    cnf['split_samples'] = cnf.get('split_samples', False)

    if options.get('threads'):
        if 'gatk' in cnf:
            gatk_opts = cnf['gatk'].get('options', [])
            new_opts = []
            for opt in gatk_opts:
                if opt.startswith('-nt '):
                    new_opts.append('-nt ' + options.get('threads'))
                else:
                    new_opts.append(opt)
            cnf['gatk']['options'] = new_opts

    if options.get('rewrite'):
        cnf['reuse_intermediate'] = False

    if options.get('no_parallel'):
        cnf['parallel'] = False

    return cnf, options


def _load_configs(run_config_name, system_config_name, args):
    system_config_path = join(dirname(dirname(realpath(__file__))), system_config_name)
    run_config_path = join(dirname(dirname(realpath(__file__))), run_config_name)
    if len(args) < 1:
        sys.stderr.write('Notice: using ' + run_config_name + ' as a default'
                         ' run configutation file.\n\n')
    else:
        run_config_path = args[0]

    if len(args) < 2:
        sys.stderr.write('Notice: using system_info_rask.yaml as a default'
                         ' tools configutation file.\n\n')
    else:
        system_config_path = args[0]
        run_config_path = args[1]

    if not os.path.isfile(system_config_path):
        exit(system_config_path + ' does not exist or is a directory.\n')
    if not os.path.isfile(run_config_path):
        exit(run_config_path + ' does not exist or is a directory.\n')

    to_exit = False
    if not system_config_path.endswith('.yaml'):
        sys.stderr.write(system_config_path + ' does not end with .yaml,'
                                              ' maybe incorrect parameter?\n')
        to_exit = True
    if not run_config_path.endswith('.yaml'):
        sys.stderr.write(run_config_path + ' does not end with .yaml,'
                                           ' maybe incorrect parameter?\n')
        to_exit = True

    if to_exit:
        exit(1)

    _check_system_tools()

    sys_cnf = load(open(system_config_path), Loader=Loader)
    run_cnf = load(open(run_config_path), Loader=Loader)
    info('Loaded system config ' + system_config_path)
    info('Loaded run config ' + run_config_path)

    return dict(run_cnf.items() + sys_cnf.items())


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
                critical(cnf.get('log'), 'No "resources" section in system config.')

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

        genome_cnf['seq'] = expanduser(genome_cnf['seq'])
        if not verify_file(genome_cnf['seq'], 'Reference seq'):
            to_exit = True

    if 'snpeff' in required:
        required.remove('snpeff')
        if 'snpeff' in genome_cnf:
            genome_cnf['snpeff'] = expanduser(genome_cnf['snpeff'])
            if not verify_dir(genome_cnf['snpeff'], 'snpeff'):
                to_exit = True

    for f in required:  # 'dbsnp', 'cosmic', 'dbsnfp', '1000genomes':
        if f not in genome_cnf:
            err('Please, provide path to ' + f  + ' in system config genome section.')
            to_exit = True
        else:
            genome_cnf[f] = expanduser(genome_cnf[f])
            if not verify_file(genome_cnf[f], f):
                to_exit = True

    for f in optional:
        if f not in genome_cnf:
            continue
        else:
            genome_cnf[f] = expanduser(genome_cnf[f])
            if not verify_file(genome_cnf[f], f):
                to_exit = True

    if to_exit:
        exit(1)

    cnf['genome'] = genome_cnf
    del cnf['genomes']
    genome_cnf['name'] = genome_name

    info('Loaded resources for ' + genome_cnf['name'])


def _check_system_tools():
    to_exit = False
    if not which('java'):
        err('\n* Warning: Java not found. You may want to run "module load java", '
            'or better ". /group/ngs/bin/bcbio-prod.sh"\n')
        to_exit = True

    if not which('perl'):
        err('\n* Warning: Perl not found. You may want to run "module load perl", '
            'or better ". /group/ngs/bin/bcbio-prod.sh"\n')

    # print ''
    # print 'In Waltham, run this as well:'
    # print '   export PATH=$PATH:/group/ngs/source/snpEff/snpEff3.5/scripts'
    # print '   export PERL5LIB=$PERL5LIB:/opt/az/local/bcbio-nextgen/'
    #       'stable/0.7.6/tooldir/lib/perl5/site_perl'
    if to_exit:
        exit()


def _set_up_dirs(cnf):
    cnf['output_dir'] = expanduser(cnf['output_dir'])
    output_dirpath = realpath(cnf['output_dir'])
    safe_mkdir(output_dirpath, 'output_dir')
    work_dirpath = join(output_dirpath, 'work')
    safe_mkdir(work_dirpath, 'working directory')
    cnf['work_dir'] = work_dirpath
    cnf['log'] = join(cnf['work_dir'], 'log.txt')