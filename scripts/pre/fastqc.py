#!/usr/bin/env python
from itertools import izip
import __check_python_version

from os.path import join, basename, dirname, splitext
import sys
from source import info, verify_file
from optparse import OptionParser
from source.fastqc.summarize_fastqc import write_fastqc_combo_report
from source.file_utils import adjust_path, verify_dir
from source.logger import critical
from source.tools_from_cnf import get_system_path


class Sample:
    def __init__(self, name, l_fpath=None, r_fpath=None, fastqc_html_fpath=None):
        self.name = name
        self.l_fpath = l_fpath
        self.r_fpath = r_fpath
        self.fastqc_html_fpath = None


def _get_sample_name(lp, rp):
    return splitext(''.join(lc if lc == rc else '' for lc, rc in izip(lp, rp)))[0]


def run_fastq(sample_name, left_reads_fpath, right_reads_fpath):
    sample = Sample(sample_name, left_reads_fpath, right_reads_fpath)

    fastqc = get_system_path(cnf, 'fastqc')
    if not fastqc:
        err('FastQC is not found, cannot make reports')
    else:
        safe_mkdir(fastqc_dirpath)
        fastqc_jobs = []
        for sample in samples:
            j = run_fastqc(cnf, sample, fastqc_dirpath)
            fastqc_jobs.append(j)
        wait_for_jobs(fastqc_jobs)

    for s in samples:
        sample_fastqc_dirpath = join(fastqc_dirpath, s.name + '.fq_fastqc')
        s.fastqc_html_fpath = join(sample_fastqc_dirpath, 'fastqc_report.html')
        verify_file(s.fastqc_html_fpath, is_critical=True, silent=True)
        if os.path.exists(sample_fastqc_dirpath + '.zip'):
            os.remove(sample_fastqc_dirpath + '.zip')

    return results_dirpath


def main():
    info(' '.join(sys.argv))
    info()

    description = 'This script runs preprocessing.'

    parser = OptionParser(description=description)
    parser.add_option('-r1', dest='left_reads_fpath', help='Left reads fpath')
    parser.add_option('-r2', dest='right_reads_fpath', help='Right reads fpath')
    parser.add_option('-o', dest='output_dir', help='Output directory path')
    (opts, args) = parser.parse_args()

    left_reads_fpath = verify_file(opts.left_reads_fpath, is_critical=True)
    right_reads_fpath = verify_file(opts.right_reads_fpath, is_critical=True)
    output_dir = opts.output_dir if opts.output_dir else critical('Please, specify output directory with -o')
    verify_dir(dirname(output_dir), is_critical=True)

    sample_name = _get_sample_name(left_reads_fpath, right_reads_fpath)
    results_dirpath = run_fastq(sample_name, left_reads_fpath, right_reads_fpath)

    verify_dir(results_dirpath, is_critical=True)

    info()
    info('*' * 70)
    info('Fastqc results:')
    info('  ' + results_dirpath)


if __name__ == '__main__':
    main()
