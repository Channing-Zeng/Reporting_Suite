#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc


from os.path import join, basename, dirname, splitext, isfile
from itertools import izip
import os
import sys
from source import logger
from source import info, verify_file
from optparse import OptionParser
from source.calling_process import call
from source.config import Config
from source.fastqc.fastq_utils import downsample
from source.file_utils import adjust_path, verify_dir, workdir, open_gzipsafe, safe_mkdir
from source.logger import critical, warn, err
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug, determine_sys_cnf, \
    determine_run_cnf
from source.tools_from_cnf import get_system_path


class Sample:
    def __init__(self, name, l_fpath=None, r_fpath=None, fastqc_html_fpath=None):
        self.name = name
        self.l_fpath = l_fpath
        self.r_fpath = r_fpath
        self.fastqc_html_fpath = None


def main():
    info(' '.join(sys.argv))
    info()

    description = 'This script runs preprocessing.'

    parser = OptionParser(description=description)
    parser.add_option('-1', dest='left_reads_fpath', help='Left reads fpath')
    parser.add_option('-2', dest='right_reads_fpath', help='Right reads fpath')
    parser.add_option('--sample', dest='sample_name', help='Sample name')
    parser.add_option('--suffix', dest='suffix', default='subset', help='Output files suffix')
    parser.add_option('-o', dest='output_dir', help='Output directory path')
    parser.add_option('--downsample-to', dest='downsample_to', default=5e5, type='int',
        help='Downsample reads to avoid excessive processing times with large files. '
             'Default is 1 million. Set to 0 to turn off downsampling.')
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser)
    (opts, args) = parser.parse_args()
    logger.is_debug = opts.debug

    cnf = Config(opts.__dict__, determine_sys_cnf(opts), determine_run_cnf(opts))
    left_reads_fpath = verify_file(opts.left_reads_fpath, is_critical=True)
    right_reads_fpath = verify_file(opts.right_reads_fpath, is_critical=True) if opts.right_reads_fpath else None
    output_dirpath = adjust_path(opts.output_dir) if opts.output_dir else critical('Please, specify output directory with -o')
    safe_mkdir(output_dirpath)
    verify_dir(dirname(output_dirpath), description='output_dir', is_critical=True)

    with workdir(cnf):
        info('Downsampling to ' + str(cnf.downsample_to))
        downsample(
            cnf, left_reads_fpath, right_reads_fpath, cnf.downsample_to,
            output_dir=cnf.output_dir, suffix=cnf.suffix)


if __name__ == '__main__':
    main()
