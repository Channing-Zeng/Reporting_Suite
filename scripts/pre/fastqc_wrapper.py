#!/usr/bin/env python
import __check_python_version

from os.path import join, basename, dirname, splitext, isfile
from itertools import izip
import os
import sys
from source import info, verify_file
from optparse import OptionParser
from source.calling_process import call
from source.config import Config
from source.fastqc.fastq_utils import downsample
from source.file_utils import adjust_path, verify_dir, workdir, open_gzipsafe
from source.logger import critical, warn, err
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_reuse_marker_genome, determine_sys_cnf, \
    determine_run_cnf
from source.tools_from_cnf import get_system_path, get_script_cmdline
from source.webserver.ssh_utils import connect_to_server


class Sample:
    def __init__(self, name, l_fpath=None, r_fpath=None, fastqc_html_fpath=None):
        self.name = name
        self.l_fpath = l_fpath
        self.r_fpath = r_fpath
        self.fastqc_html_fpath = None


REPORTING_SUITE_PATH_CLARITY = '/opt/gls/clarity/customextensions/reporting/Reporting_Suite-master'
REPORTING_SUITE_PATH_WALTHAM = '/group/ngs/src/az.reporting-dev'

PYTHON_PATH_CLARITY = '/usr/local/bin/python'
PYTHON_PATH_WALTHAM = '/group/ngs/src/bcbio-nextgen/0.8.7/rhel6-x64/anaconda/bin/python'


def main():
    info(' '.join(sys.argv))
    info()

    description = 'This script runs preprocessing.'

    parser = OptionParser(description=description)
    parser.add_option('-1', dest='left_reads_fpath', help='Left reads fpath')
    parser.add_option('-2', dest='right_reads_fpath', help='Right reads fpath')
    parser.add_option('--sample', dest='sample_name', help='Sample name')
    parser.add_option('-o', dest='output_dir', help='Output directory path')
    parser.add_option('--downsample-to', dest='downsample_to', default=None, type='int',
        help='Downsample reads to avoid excessive processing times with large files. '
            'Default is 1 million. Set to 0 to turn off downsampling.')
    add_cnf_t_reuse_prjname_reuse_marker_genome(parser)
    (opts, args) = parser.parse_args()

    if not opts.left_reads_fpath or not opts.right_reads_fpath or not opts.output_dir:
        parser.print_usage()

    verify_file(opts.left_reads_fpath, is_critical=False)
    left_reads_fpath = adjust_path(opts.left_reads_fpath)
    verify_file(opts.right_reads_fpath, is_critical=False)
    right_reads_fpath = adjust_path(opts.right_reads_fpath)
    output_dirpath = adjust_path(opts.output_dir) if opts.output_dir else critical('Please, specify output directory with -o')
    verify_dir(dirname(output_dirpath), description='output_dir', is_critical=True)

    left_reads_fpath, right_reads_fpath, output_dirpath =\
        map(_proc_path, [left_reads_fpath, right_reads_fpath, output_dirpath])

    ssh = connect_to_server(server_url='blue.usbod.astrazeneca.net', username='klpf990', password='123qweasd')
    fastqc_py = get_script_cmdline(None, 'python', 'scripts/pre/fastqc.py')
    fastqc_py = fastqc_py.replace(REPORTING_SUITE_PATH_CLARITY, REPORTING_SUITE_PATH_WALTHAM)
    fastqc_py = fastqc_py.replace(PYTHON_PATH_CLARITY, PYTHON_PATH_WALTHAM)

    cmdl = '{fastqc_py} -1 {left_reads_fpath} -2 {right_reads_fpath} -o {output_dirpath}'
    if opts.sample_name:
        cmdl += ' --sample {opts.sample_name}'
    if opts.downsample_to:
        cmdl += ' --downsample-to ' + str(int(opts.downsample_to))
    cmdl = cmdl.format(**locals())
    cmdl += ' 2>&1'
    info(cmdl)
    stdin, stdout, stderr = ssh.exec_command(cmdl)
    for l in stdout:
        err(l, ending='')
    info()
    ssh.close()


def _proc_path(path):
    starts = {'/mnt/Datasets': '/ngs/oncology/datasets',
              '/mnt/HiSeq': '/ngs/oncology/datasets/HiSeq/',
              '/mnt/MiSeq': '/ngs/oncology/datasets/MiSeq/'}
    if not any(path.startswith(s) for s in starts.keys()):
        critical('Error: path ' + path + ' has to start with something from ' + str(starts.keys()))
    for k, v in starts.iteritems():
        path = path.replace(k, v)
    return path


if __name__ == '__main__':
    main()
