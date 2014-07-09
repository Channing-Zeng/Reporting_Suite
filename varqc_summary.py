#!/usr/bin/env python

from __future__ import print_function
import sys

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from os.path import join
from optparse import OptionParser

from source.reporting import read_sample_names
from source.variants.summarize_qc import make_summary_reports
from source.config import Defaults, Config
from source.main import check_keys, check_inputs, set_up_dirs
from source.logger import info


def main():
    description = 'This script generates project-level summaries based on per-sample targetcov reports.'

    parser = OptionParser(description=description)
    parser.add_option('-d', dest='bcbio_final_dir', help='Path to bcbio-nextgen final directory (default is pwd)')
    parser.add_option('-s', dest='samples', help='List of samples (default is samples.txt in bcbio final directory)')
    parser.add_option('-n', dest='base_name', default='varqc', help='Name of variant QC directory inside sample folder. (default is VarQC)')
    parser.add_option('--vcf-suf', dest='vcf_suf', help='Suffix to choose VCF files (mutect, ensembl, freebayes, etc). Multiple comma-separated values allowed.')
    parser.add_option('-o', '--output_dir', dest='output_dir', metavar='DIR', help='output directory (or directory name in case of bcbio final dir)')

    parser.add_option('-v', dest='verbose', action='store_true', help='Verbose')
    parser.add_option('-t', dest='threads', type='int', help='Number of threads for each process')
    parser.add_option('-w', dest='overwrite', action='store_true', help='Overwrite existing results')

    parser.add_option('--runner', dest='qsub_runner', help='Bash script that takes command line as the 1st argument. This script will be submitted to GRID. Default: ' + Defaults.qsub_runner)
    parser.add_option('--sys-cnf', dest='sys_cnf', default=Defaults.sys_cnf, help='system configuration yaml with paths to external tools and genome resources (see default one %s)' % Defaults.sys_cnf)
    parser.add_option('--run-cnf', dest='run_cnf', default=Defaults.run_cnf, help='run configuration yaml (see default one %s)' % Defaults.run_cnf)

    parser.add_option('--work-dir', dest='work_dir', metavar='DIR')

    (opts, args) = parser.parse_args()
    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)

    cnf.name = cnf['name'] or 'varqc_summary'
    set_up_dirs(cnf)

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

    sample_names = read_sample_names(cnf['samples'])

    make_summary_reports(cnf, sample_names)


if __name__ == '__main__':
    main()















