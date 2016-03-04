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
    parser.add_option('-o', dest='output_dir', help='Output directory path')
    parser.add_option('--downsample-to', dest='downsample_to', default=None, type='int',
        help='Downsample reads to avoid excessive processing times with large files. '
             'Default is 1 million. Set to 0 to turn off downsampling.')
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser, threads=1)
    (opts, args) = parser.parse_args()
    logger.is_debug = opts.debug

    cnf = Config(opts.__dict__, determine_sys_cnf(opts), determine_run_cnf(opts))
    left_reads_fpath = verify_file(opts.left_reads_fpath, is_critical=True)
    right_reads_fpath = verify_file(opts.right_reads_fpath, is_critical=True)
    output_dirpath = adjust_path(opts.output_dir) if opts.output_dir else critical('Please, specify output directory with -o')
    verify_dir(dirname(output_dirpath), description='output_dir', is_critical=True)

    with workdir(cnf):
        sample_name = cnf.sample_name
        if not sample_name:
            sample_name = _get_sample_name(left_reads_fpath, right_reads_fpath)
        results_dirpath = run_fastq(cnf, sample_name, left_reads_fpath, right_reads_fpath, output_dirpath, downsample_to=cnf.downsample_to)

    verify_dir(results_dirpath, is_critical=True)
    info()
    info('*' * 70)
    info('Fastqc results:')
    info('  ' + results_dirpath)


def run_fastq(cnf, sample_name, l_r_fpath, r_r_fpath, output_dirpath, downsample_to=1e7):
    fastqc = get_system_path(cnf, 'fastqc', is_critical=True)
    java = get_system_path(cnf, 'java', is_critical=True)

    if downsample_to:
        info('Downsampling to ' + str(downsample_to))
        l_fpath, r_fpath = downsample(cnf, l_r_fpath, r_r_fpath, downsample_to, output_dir=cnf.work_dir)

    # Joining fastq files to run on a combination
    fastqc_fpath = join(cnf.work_dir, sample_name + '.fq')
    info('Combining fastqs, writing to ' + fastqc_fpath)
    with open(fastqc_fpath, 'w') as out:
        out.write(open_gzipsafe(l_r_fpath).read())
        out.write(open_gzipsafe(r_r_fpath).read())

    # Running FastQC
    info('Running FastQC')
    tmp_dirpath = join(cnf.work_dir, 'FastQC_' + sample_name + '_tmp')
    safe_mkdir(tmp_dirpath)
    cmdline = '{fastqc} --dir {tmp_dirpath} --extract -o {output_dirpath} -f fastq -j {java} {fastqc_fpath}'.format(**locals())
    call(cnf, cmdline)

    # Cleaning and getting report
    sample_fastqc_dirpath = join(output_dirpath, sample_name + '.fq_fastqc')
    if isfile(sample_fastqc_dirpath + '.zip'):
        os.remove(sample_fastqc_dirpath + '.zip')
    fastqc_html_fpath = join(sample_fastqc_dirpath, 'fastqc_report.html')
    verify_file(fastqc_html_fpath, is_critical=True)

    return sample_fastqc_dirpath


def _get_sample_name(lp, rp):
    return splitext(''.join(lc if lc == rc else '' for lc, rc in izip(lp, rp)))[0]


if __name__ == '__main__':
    main()
