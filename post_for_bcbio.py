#!/usr/bin/env python
import sys

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from optparse import OptionParser
from os.path import join, abspath, dirname, pardir, isdir

from source.config import Defaults, Config
from source.logger import info
from source.main import check_system_resources, check_inputs, check_keys
from source.bcbio_runner import run_on_bcbio_final_dir
from source.utils_from_bcbio import safe_mkdir


def main():
    info(' '.join(sys.argv))
    info()

    description = 'This script runs reporting suite on the bcbio final directory.'

    parser = OptionParser(description=description)
    parser.add_option('-d', '-o', dest='bcbio_final_dir', help='Path to bcbio-nextgen final directory (default is pwd)')
    parser.add_option('-s', '--samples', dest='samples', help='List of samples (default is samples.txt in bcbio final directory)')
    parser.add_option('-b', '--bed', dest='bed', help='BED file')
    parser.add_option('--vcf-suf', '--vcf-suffix', dest='vcf_suf', help='Suffix to choose VCF files (mutect, ensembl, freebayes, etc). Multiple comma-separated values allowed.')
    parser.add_option('--qualimap', dest='qualimap', action='store_true', default=Defaults.qualimap, help='Run QualiMap in the end')

    parser.add_option('-v', dest='verbose', action='store_true', help='Verbose')
    parser.add_option('-t', dest='threads', type='int', help='Number of threads for each process')
    parser.add_option('-w', dest='overwrite', action='store_true', help='Overwrite existing results')
    parser.add_option('--reuse', dest='overwrite', action='store_false', help='Reuse intermediate results in work directory for subroutines')

    parser.add_option('--runner', dest='qsub_runner', help='Bash script that takes command line as the 1st argument. This script will be submitted to GRID. Default: ' + Defaults.qsub_runner)
    parser.add_option('--sys-cnf', dest='sys_cnf', default=Defaults.sys_cnf, help='system configuration yaml with paths to external tools and genome resources (see default one %s)' % Defaults.sys_cnf)
    parser.add_option('--run-cnf', dest='run_cnf', default=Defaults.run_cnf, help='run configuration yaml (see default one %s)' % Defaults.run_cnf)

    (opts, args) = parser.parse_args()
    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)

    if not check_keys(cnf, ['bcbio_final_dir', 'bed', 'vcf_suf']):
        parser.print_help()
        sys.exit(1)

    if not check_inputs(cnf, dir_keys=['bcbio_final_dir']):
        sys.exit(1)

    if isdir(join(cnf.bcbio_final_dir, 'final')):
        cnf.bcbio_final_dir = join(cnf.bcbio_final_dir, 'final')

    if not cnf.samples:
        cnf.samples = join(cnf.bcbio_final_dir, 'samples.txt')

    cur_dirpath = dirname(abspath(__file__))
    cnf.qsub_runner = abspath(join(cur_dirpath, cnf.qsub_runner))

    if not check_inputs(cnf, file_keys=['samples', 'bed', 'qsub_runner'], dir_keys=['bcbio_final_dir']):
        sys.exit(1)

    info('BCBio "final" dir: ' + cnf.bcbio_final_dir + ' (set with -d)')
    info('Samples: ' + cnf.samples + ' (set with -s)')
    info('Capture/amplicons BED file: ' + cnf.bed + ' (set with -b)')
    info('Suffix(es) to choose VCF files: ' + cnf.vcf_suf + ' (set with --vcf-suf)')
    info()
    info('*' * 70)

    if opts.qualimap and 'QualiMap' not in cnf.steps:
        cnf.steps.append('QualiMap')

    check_system_resources(cnf, required=['qsub'])

    vcf_sufs = cnf['vcf_suf'].split(',')

    cnf['work_dir'] = join(cnf.bcbio_final_dir, pardir, 'work', 'post_processing')
    if not isdir(cnf.work_dir):
        safe_mkdir(cnf.work_dir)
    # cnf['bcbio_cnf'] = join(cnf.bcbio_final_dir, pardir, 'config')

    run_on_bcbio_final_dir(cnf, cnf.bcbio_final_dir, cnf.samples, cnf.bed, vcf_sufs)


if __name__ == '__main__':
    main()









