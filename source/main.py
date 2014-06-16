#!/usr/bin/env python
import sys
from os.path import isdir, join, realpath, expanduser, basename, abspath
from optparse import OptionParser
from shutil import rmtree

from source import logger
from source.config import Config, Defaults
from source.logger import info, err, critical
from source.utils_from_bcbio import which, file_exists
from source.utils import verify_file, safe_mkdir, verify_dir, make_tmpdir
from source.ngscat.bed_file import verify_bam, verify_bed


def read_opts_and_cnfs(extra_opts,
                       key_for_sample_name,
                       required_keys,
                       file_keys=list(),
                       dir_keys=list()):
    options = extra_opts + [
        (['-o', '--output_dir'], dict(
             dest='output_dir',
             metavar='DIR',
             help='output directory (or directory name in case of bcbio final dir)')
         ),
        (['-s', '--sample'], dict(
             dest='name',
             metavar='NAME',
             help='sample name (default is part of name of the first parameter prior to the first - or .')
         ),
        (['-t', '--nt', '--threads'], dict(
             dest='threads',
             type='int',
             help='number of threads')
         ),
        (['-w', '--overwrite'], dict(
             dest='overwrite',
             help='do not reuse intermediate files from previous run',
             action='store_true')
         ),
        (['--sys-cnf'], dict(
             dest='sys_cnf',
             metavar='SYS_CNF.yaml',
             default=Defaults.sys_cnf,
             help='system configuration yaml with paths '
                  'to external tools and genome resources '
                  '(see default one %s)' % Defaults.sys_cnf)
         ),
        (['--run-cnf'], dict(
             dest='run_cnf',
             metavar='RUN_CNF.yaml',
             default=Defaults.run_cnf,
             help='run configuration yaml (see default one %s)'
                  % Defaults.run_cnf)
         ),
    ]

    parser = OptionParser()
    for args, kwargs in options:
        parser.add_option(*args, **kwargs)

    (opts, args) = parser.parse_args()
    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)

    if not check_keys(cnf, required_keys):
        parser.print_help()
        sys.exit(1)

    if not check_inputs(cnf, file_keys, dir_keys):
        sys.exit(1)

    assert key_for_sample_name and cnf[key_for_sample_name]
    if key_for_sample_name not in cnf:
        critical('Error: ' + key_for_sample_name + ' must be provided '
                 'in options or in ' + cnf.run_cnf + '.')

    key_fname = basename(cnf.get(key_for_sample_name))
    cnf.name = cnf.get('name') or key_fname.split('.')[0].split('-')[0]

    set_up_dirs(cnf)

    return cnf


def check_keys(cnf, required_keys):
    to_exit = False

    for key in required_keys:
        if key not in cnf or not cnf[key]:
            to_exit = True
            err('Error: "' + key + '" must be provided in options or '
                'in ' + cnf.run_cnf + '.')
    return not to_exit


def check_inputs(cnf, file_keys=list(), dir_keys=list()):
    to_exit = False

    def _verify_input_file(_key):
        if not verify_file(cnf[_key], _key):
            return False
        if 'bam' in _key and not verify_bam(cnf[_key]):
            return False
        if 'bed' in _key and not verify_bed(cnf[_key]):
            return False
        return True

    for key in file_keys:
        if key and key in cnf and cnf.key:
            if not _verify_input_file(key):
                to_exit = True
            else:
                cnf[key] = abspath(expanduser(cnf[key]))

    for key in dir_keys:
        if key and key in cnf and cnf.key:
            if not verify_dir(cnf[key], key):
                to_exit = True
            else:
                cnf[key] = abspath(expanduser(cnf[key]))

    return not to_exit


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
            if cnf.resources is None:
                critical('No "resources" section in system config.')

            data = cnf.resources.get(program)
            if data is None:
                err(program + ' is required. Specify path in system config or in your environment.')
                to_exit = True
            else:
                data['path'] = expanduser(data['path'])
                if not isdir(data['path']) and not file_exists(data['path']):
                    err(data['path'] + ' does not exist.')
                    to_exit = True

    for program in optional:
        resources = cnf.get('resources')
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
    genome_name = cnf.genome
    if genome_name not in cnf.genomes:
        critical(genome_name + ' is not in "genomes section" of system config.')
    genome_cnf = cnf.genomes[genome_name].copy()

    to_exit = False

    for key in genome_cnf:
        genome_cnf[key] = expanduser(genome_cnf[key])

    for key in required:  # 'dbsnp', 'cosmic', 'dbsnfp', '1000genomes':
        if key not in genome_cnf:
            if key == 'seq':
                err('Please, provide path to the reference file (seq).')
            else:
                err('Please, provide path to ' + key + ' in system config genome section.')
            to_exit = True
        else:
            genome_cnf[key] = abspath(expanduser(genome_cnf[key]))
            if key == 'snpeff':
                if not verify_dir(genome_cnf['snpeff'], 'snpeff'):
                    to_exit = True
            else:
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

    cnf.genome = genome_cnf
    del cnf.__dict__['genomes']
    genome_cnf['name'] = genome_name

    info('Loaded resources for ' + genome_cnf['name'])


def set_up_dirs(cnf):
    """ Creates output_dir, work_dir; sets up log
    """
    cnf.output_dir = realpath(expanduser(cnf.output_dir))

    safe_mkdir(cnf.output_dir, 'output_dir')
    info('Saving into ' + cnf.output_dir)

    work_dir_name = 'work_' + cnf.name
    cnf.work_dir = join(cnf.output_dir, work_dir_name)
    if not cnf.reuse_intermediate and isdir(cnf.work_dir):
        rmtree(cnf.work_dir)

    safe_mkdir(cnf.work_dir, 'working directory')

    cnf.log = join(cnf.work_dir, cnf.name + '_log.txt')
    logger.log_fname = cnf.log