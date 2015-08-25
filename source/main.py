#!/usr/bin/env python
import os
import sys
from os.path import isdir, join, realpath, expanduser, basename, abspath, dirname, pardir, isfile
from optparse import OptionParser
from shutil import rmtree
from source.bcbio_structure import ungzip_if_needed

from source.file_utils import verify_file, verify_dir, adjust_path, remove_quotes, adjust_system_path, \
    verify_obj_by_path
from source import logger
from source.config import Config, defaults
from source.logger import info, err, critical
from source.prepare_args_and_cnf import set_up_dirs, check_inputs, check_keys, determine_run_cnf, \
    determine_sys_cnf


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
             dest='sample',
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
        (['--work-dir'], dict(dest='work_dir', metavar='DIR')),
        (['--log-dir'], dict(dest='log_dir', metavar='DIR')),
        (['--proc-name'], dict(dest='proc_name')),
        (['--project-name'], dict(dest='project_name')),
        (['--genome'], dict(dest='genome', default=defaults['genome'])),
        (['--email'], dict(dest='email')),
        (['--done-marker'], dict(dest='done_marker')),
    ]

    parser = OptionParser(description=description)
    for args, kwargs in options:
        parser.add_option(*args, **kwargs)

    (opts, args) = parser.parse_args()

    run_cnf = determine_run_cnf(opts, is_wgs=not opts.__dict__.get('bed'))
    cnf = Config(opts.__dict__, determine_sys_cnf(opts), run_cnf)

    errors = check_keys(cnf, required_keys)
    if errors:
        parser.print_help()
        critical(errors)
    errors = check_inputs(cnf, file_keys, dir_keys)
    if errors:
        critical(errors)

    if cnf.sample:
        cnf.sample = remove_quotes(cnf.sample)
    else:
        if not key_for_sample_name or not cnf[key_for_sample_name]:
            if cnf.sample:
                critical('Error: ' + (key_for_sample_name or 'key_for_sample_name') +
                    ' must be provided in options or in ' + cnf.run_cnf + '.')
        key_fname = basename(cnf[key_for_sample_name])
        cnf.sample = key_fname.split('.')[0]

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
