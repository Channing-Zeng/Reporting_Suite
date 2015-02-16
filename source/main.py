#!/usr/bin/env python
import sys
from os.path import isdir, join, realpath, expanduser, basename, abspath, dirname, pardir
from optparse import OptionParser
from shutil import rmtree
from source.bcbio_structure import ungzip_if_needed

from source.file_utils import verify_file, verify_dir, adjust_path, remove_quotes, adjust_system_path, \
    verify_obj_by_path
from source import logger
from source.config import Config, defaults
from source.logger import info, err, critical
from source.file_utils import which, file_exists, safe_mkdir
from source.ngscat.bed_file import verify_bam, verify_bed


code_base_path = abspath(join(dirname(abspath(__file__)), pardir))


def read_opts_and_cnfs(extra_opts,
                       key_for_sample_name,
                       required_keys,
                       file_keys=list(),
                       dir_keys=list(),
                       description=None,
                       extra_msg=None,
                       proc_name=None):
    options = extra_opts + [
        (['-o', '--output_dir'], dict(
             dest='output_dir',
             metavar='DIR',
             help='output directory (or directory name in case of bcbio final dir)')
         ),
        (['-s', '--sample', '--name'], dict(
             dest='name',
             metavar='NAME',
             help='sample name (default is part of name of the first parameter prior to the first - or .')
         ),
        (['-c', '--caller'], dict(
             dest='caller',
             metavar='CELLR',
             help='variant caller name (default is part of name of the first parameter between the first - and following .')
         ),
        (['-t', '--nt', '--threads'], dict(
             dest='threads',
             type='int',
             help='number of threads')
         ),
        (['--clean'], dict(
             dest='keep_intermediate',
             help='do not store work directory',
             action='store_false')
         ),
        (['--reuse'], dict(
             dest='reuse_intermediate',
             help='reuse intermediate non-empty files in the work dir from previous run',
             action='store_true')
         ),
        (['--sys-cnf'], dict(
             dest='sys_cnf',
             metavar='SYS_CNF.yaml',
             default=defaults['sys_cnfs']['uk'],
             help='System configuration file with paths to external tools and genome resources. The default is  '
                  '(see default one %s)' % defaults['sys_cnf'])
         ),
        (['--run-cnf'], dict(
             dest='run_cnf',
             metavar='RUN_CNF.yaml',
             default=defaults['run_cnf'],
             help='Customised run details: list of annotations/QC metrics/databases/filtering criteria. '
                  'The default is %s' % defaults['run_cnf'])
         ),
        (['--work-dir'], dict(dest='work_dir', metavar='DIR')),
        (['--log-dir'], dict(dest='log_dir', metavar='DIR')),
        (['--proc-name'], dict(dest='proc_name')),
        (['--project-name'], dict(dest='project_name')),
        (['--genome'], dict(dest='genome', default=defaults['genome'])),
        (['--email'], dict(dest='email')),
    ]

    parser = OptionParser(description=description)
    for args, kwargs in options:
        parser.add_option(*args, **kwargs)

    (opts, args) = parser.parse_args()
    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)

    if not check_keys(cnf, required_keys):
        parser.print_help()
        sys.exit(1)

    if not check_inputs(cnf, file_keys, dir_keys):
        sys.exit(1)

    if cnf.name:
        cnf.name = remove_quotes(cnf.name)
    else:
        if not key_for_sample_name or not cnf[key_for_sample_name]:
            if cnf.name:
                critical('Error: ' + (key_for_sample_name or 'key_for_sample_name') + ' must be provided '
                         'in options or in ' + cnf.run_cnf + '.')
        key_fname = basename(cnf[key_for_sample_name])
        cnf.name = key_fname.split('.')[0]

    if cnf.caller:
        cnf.caller = remove_quotes(cnf.caller)
    elif key_for_sample_name and cnf[key_for_sample_name]:
        key_fname = basename(cnf[key_for_sample_name])
        try:
            cnf.caller = cnf.caller or key_fname.split('.')[0].split('-')[1]
        except:
            cnf.caller = ''
    else:
        cnf.caller = ''

    cnf.proc_name = cnf.proc_name or proc_name
    set_up_dirs(cnf)
    info(' '.join(sys.argv))
    info()

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
        cnf[_key] = adjust_path(cnf[_key])
        if not verify_file(cnf[_key], _key):
            return False
        if 'bam' in _key and not verify_bam(cnf[_key]):
            return False
        if 'bed' in _key and not verify_bed(cnf[_key]):
            return False
        return True

    for key in file_keys:
        if key and key in cnf and cnf[key]:
            if not _verify_input_file(key):
                to_exit = True
            else:
                cnf[key] = adjust_path(cnf[key])

    for key in dir_keys:
        if key and key in cnf and cnf[key]:
            cnf[key] = adjust_system_path(cnf[key])
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
                data['path'] = adjust_system_path(data['path'])
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
            data['path'] = adjust_system_path(data['path'])
            if not isdir(data['path']) and not file_exists(data['path']):
                err(data['path'] + ' does not exist.')
                to_exit = True

    if to_exit:
        exit()


def check_genome_resources(cnf):
    if not cnf.genomes:
        critical('"genomes" section is not specified in system config.')

    info('Checking paths in the genomes sections in ' + cnf.sys_cnf)
    info()

    for build_name, genome_cnf in cnf.genomes.items():
        info(build_name)
        for key in genome_cnf.keys():
            if isinstance(genome_cnf[key], basestring):
                genome_cnf[key] = adjust_system_path(genome_cnf[key])

            if not verify_obj_by_path(genome_cnf[key], key):
                if not genome_cnf[key].endswith('.gz') and verify_file(genome_cnf[key] + '.gz'):
                    gz_fpath = genome_cnf[key] + '.gz'
                    if verify_file(gz_fpath):
                        info(key + ': ' + gz_fpath)
                        genome_cnf[key] = gz_fpath
                else:
                    err('   err: no ' + genome_cnf[key] + (' and .gz' if not genome_cnf[key].endswith('gz') else ''))
            else:
                info(key + ': ' + genome_cnf[key])
        info()
        genome_cnf['name'] = build_name

    cnf.genome = cnf.genomes[cnf.genome]

    info('Checked genome resources.')
    info('*' * 70)
    info()


def set_up_dirs(cnf):
    """ Creates output_dir, work_dir; sets up log
    """
    cnf.output_dir = adjust_path(cnf.output_dir)
    safe_mkdir(cnf.output_dir, 'output_dir')
    info('Saving into ' + cnf.output_dir)

    set_up_work_dir(cnf)
    set_up_log(cnf)


def set_up_work_dir(cnf):
    if not cnf.work_dir:
        work_dir_name = 'work_' + cnf.name
        cnf.work_dir = join(cnf.output_dir, work_dir_name)
        # if not cnf.reuse_intermediate and isdir(cnf.work_dir):
        #     rmtree(cnf.work_dir)
    else:
        cnf.work_dir = adjust_path(cnf.work_dir)

    safe_mkdir(cnf.work_dir, 'working directory')


def set_up_log(cnf):
    log_fname = cnf.proc_name + '_' if cnf.proc_name else ''
    log_fname += cnf.name + '_log.txt'

    cnf.log = join(cnf.work_dir, log_fname)
    logger.log_fpath = cnf.log
    logger.smtp_host = cnf.smtp_host
    logger.proc_name = cnf.proc_name
    logger.address = remove_quotes(cnf.email) if cnf.email else ''
    logger.project_name = cnf.project_name
    logger.project_fpath = cnf.output_dir