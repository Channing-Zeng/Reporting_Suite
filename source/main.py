#!/usr/bin/env python
import os
import sys
from os.path import join, basename, abspath, dirname, pardir
from optparse import OptionParser, SUPPRESS_HELP
from traceback import format_exc

from source import logger
from source.file_utils import remove_quotes
from source.config import Config, defaults
from source.logger import info, err, critical
from source.prepare_args_and_cnf import set_up_dirs, check_dirs_and_files, check_keys_presence, determine_run_cnf, \
    determine_sys_cnf


code_base_path = abspath(join(dirname(abspath(__file__)), pardir))


def read_opts_and_cnfs(extra_opts,
                       key_for_sample_name=None,
                       required_keys=list(),
                       file_keys=list(),
                       dir_keys=list(),
                       description='',
                       extra_msg=None,
                       proc_name=None,
                       fpath_for_sample_name=None,
                       with_output_file=False,
                       with_output_dir=True):
    options = extra_opts
    if with_output_file:
        options += [
            (['-o', '--output-file'], dict(
                 dest='output_file',
                 metavar='FILE',
                 help='Output file')
             )]
    elif with_output_dir:
        options += [
            (['-o', '--output-dir'], dict(
                 dest='output_dir',
                 metavar='DIR',
                 help='Output directory (or directory name in case of bcbio final dir)',
                 default=os.getcwd())
             )]

    options += [
        (['-s', '--sample', '--name'], dict(
             dest='sample',
             metavar='NAME',
             help='Sample name (default is part of name of the first parameter prior to the first - or .')
         ),
        (['-c', '--caller'], dict(
             dest='caller',
             metavar='CALLER_NAME',
             help='Variant caller name (default is part of name of the first parameter between the first - and following .')
         ),
        (['-t', '--nt', '--threads'], dict(
             dest='threads',
             type='int',
             help='Number of threads')
         ),
        (['--clean'], dict(  # do not keep work directory
             dest='keep_intermediate',
             action='store_false',
             help=SUPPRESS_HELP)
         ),
        (['--debug'], dict(
             dest='debug',
             action='store_true',
             default=False,
             help=SUPPRESS_HELP)
         ),
        (['--reuse'], dict(
             dest='reuse_intermediate',
             help='reuse intermediate non-empty files in the work dir from previous run',
             action='store_true')
         ),
        (['--sys-cnf'], dict(
             dest='sys_cnf',
             metavar='SYS_CNF.yaml',
             help='System configuration file with paths to external tools and genome resources. The default is  '
                  '(see default one %s)' % defaults['sys_cnf'])
         ),
        (['--run-cnf'], dict(
             dest='run_cnf',
             metavar='RUN_CNF.yaml',
             default=defaults['run_cnf_exome_seq'],
             help='Customised run details: list of annotations/QC metrics/databases/filtering criteria. '
                  'The default is %s' % defaults['run_cnf_exome_seq'])
         ),
        (['--work-dir'], dict(dest='work_dir', metavar='DIR', help=SUPPRESS_HELP)),
        (['--log-dir'], dict(dest='log_dir', metavar='DIR', help=SUPPRESS_HELP)),
        (['--proc-name'], dict(dest='proc_name', help=SUPPRESS_HELP)),
        (['--project-name'], dict(dest='project_name')),
        (['--no-check'], dict(dest='no_check', action='store_true', help=SUPPRESS_HELP)),
        (['--genome'], dict(dest='genome')),
        (['--email'], dict(dest='email', help=SUPPRESS_HELP)),
        (['--done-marker'], dict(dest='done_marker', help=SUPPRESS_HELP)),
    ]

    parser = OptionParser(description=description)
    for args, kwargs in options:
        parser.add_option(*args, **kwargs)

    req_keys_usage = ''
    if required_keys:
        req_keys_usage = '\nRequired options:'
    for args, kwargs in options:
        try:
            if kwargs['dest'] in required_keys:
                req_keys_usage += '\n  ' + '/'.join(args)
        except:
            err(format_exc())
            pass
    parser.set_usage(parser.get_usage() + req_keys_usage)

    (opts, args) = parser.parse_args()
    logger.is_debug = opts.debug

    run_cnf = determine_run_cnf(opts, is_wgs=not opts.__dict__.get('bed'))
    cnf = Config(opts.__dict__, determine_sys_cnf(opts), run_cnf)

    errors = check_keys_presence(cnf, required_keys)
    if errors:
        parser.print_help()
        critical(errors)
    file_keys = [k for k in file_keys if k in required_keys]
    dir_keys = [k for k in dir_keys if k in required_keys]
    errors = check_dirs_and_files(cnf, file_keys, dir_keys)
    if errors:
        critical(errors)

    if cnf.sample:
        cnf.sample = remove_quotes(cnf.sample)
    else:
        if not fpath_for_sample_name:
            if not key_for_sample_name:
                critical('Error: --sample must be provided in options.')

            fpath_for_sample_name = cnf[key_for_sample_name]
            if not fpath_for_sample_name:
                critical('Error: --sample or ' + (str(key_for_sample_name)) + ' must be provided in options.')

            key_fname = basename(cnf[key_for_sample_name])
            cnf.sample = key_fname.split('.')[0]

    if cnf.caller:
        cnf.caller = remove_quotes(cnf.caller)
    # elif key_for_sample_name and cnf[key_for_sample_name]:
    #     key_fname = basename(cnf[key_for_sample_name])
    #     try:
    #         cnf.caller = cnf.caller or key_fname.split('.')[0].split('-')[1]
    #     except:
    #         cnf.caller = ''
    else:
        cnf.caller = None

    cnf.proc_name = cnf.proc_name or proc_name
    set_up_dirs(cnf)
    info(' '.join(sys.argv))
    info()

    return cnf
