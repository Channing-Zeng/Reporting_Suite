#!/usr/bin/env python

from __future__ import print_function
from optparse import OptionParser
import sys
from os.path import join
from source.config import Defaults, Config
from source.main import check_keys, check_inputs
from source.summarize import summarize_qc
from source.utils import verify_file
from source.utils_from_bcbio import file_exists
from source.logger import info, critical

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))


def main():
    description = 'This script generates project-level summaries based on per-sample targetcov reports.'

    parser = OptionParser(description=description)
    parser.add_option('-d', dest='bcbio_final_dir', help='Path to bcbio-nextgen final directory (default is pwd)')
    parser.add_option('-s', dest='samples', help='List of samples (default is samples.txt in bcbio final directory)')
    parser.add_option('-n', dest='base_name', default='VarQC', help='Name of targetcov directory inside sample folder. (default is TargetCov)')

    parser.add_option('-v', dest='verbose', action='store_true', help='Verbose')
    parser.add_option('-t', dest='threads', type='int', help='Number of threads for each process')
    parser.add_option('-w', dest='overwrite', action='store_true', help='Overwrite existing results')

    parser.add_option('--runner', dest='qsub_runner', help='sh script for qsub that accepts command line as the 1st argument ' + Defaults.qsub_runner)
    parser.add_option('--sys-cnf', dest='sys_cnf', default=Defaults.sys_cnf, help='system configuration yaml with paths to external tools and genome resources (see default one %s)' % Defaults.sys_cnf)
    parser.add_option('--run-cnf', dest='run_cnf', default=Defaults.run_cnf, help='run configuration yaml (see default one %s)' % Defaults.run_cnf)

    (opts, args) = parser.parse_args()
    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)

    if not cnf.samples:
        cnf.samples = join(cnf.bcbio_final_dir, 'samples.txt')

    info('BCBio "final" dir: ' + cnf.bcbio_final_dir + ' (set with -d)')
    info('Samples: ' + cnf.samples + ' (set with -s)')

    if not check_keys(cnf, ['bcbio_final_dir', 'samples']):
        parser.print_help()
        sys.exit(1)

    if not check_inputs(cnf, file_keys=['samples', 'qsub_runner'], dir_keys=['bcbio_final_dir']):
        sys.exit(1)

    info()
    info('*' * 70)

    summarize_varqc_report(cnf.bcbio_final_dir, cnf.samples, cnf.base_name)


def summarize_varqc_report(out_dirpath, samples_fpath, report_basedir):
    summary_report_fpath = join(out_dirpath, 'varqc.summary.txt')
    report_suffix = '.varqc.txt'

    report_fpaths = []

    with open(samples_fpath) as f:
        for line in f:
            sample_name = line.strip()

            report_fpath = join(out_dirpath, sample_name, report_basedir, sample_name + report_suffix)
            info(sample_name + ': ' + report_fpath)

            if verify_file(report_fpath):
                report_fpaths.append(report_fpath)
            else:
                critical(report_fpath + ' does not exist.')

    summarize_qc(report_fpaths, summary_report_fpath, report_suffix)

    info()
    info('*' * 70)
    info('Result: ' + summary_report_fpath)


if __name__ == '__main__':
    main()
