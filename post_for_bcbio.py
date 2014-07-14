#!/usr/bin/env python
import sys
from source.utils_from_bcbio import safe_mkdir

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from optparse import OptionParser
from os.path import join, pardir, isdir
from os import listdir

from source.config import Defaults, Config, load_yaml_config
from source.logger import info, critical
from source.main import check_system_resources, check_inputs, check_keys
from source.bcbio_runner import run_on_bcbio_final_dir


def main():
    description = 'This script runs reporting suite on the bcbio final directory.'

    parser = OptionParser(description=description)
    parser.add_option('-d', '-o', dest='bcbio_final_dir', help='Path to bcbio-nextgen final directory (default is pwd)')
    parser.add_option('-b', '--bed', dest='bed', help='BED file')
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
    if not opts.bcbio_final_dir and len(args) > 0:
        cnf.bcbio_final_dir = args[0]

    if not check_keys(cnf, ['bcbio_final_dir']):
        parser.print_help()
        sys.exit(1)

    if not check_inputs(cnf, dir_keys=['bcbio_final_dir']):
        sys.exit(1)

    if isdir(join(cnf.bcbio_final_dir, 'final')):
        cnf.bcbio_final_dir = join(cnf.bcbio_final_dir, 'final')

    if 'qsub_runner' in cnf:
        cnf.qsub_runner = join(cnf.sys_cnf, pardir, cnf.qsub_runner)

    if not check_inputs(cnf, file_keys=['qsub_runner'], dir_keys=['bcbio_final_dir']):
        sys.exit(1)

    info('BCBio "final" dir: ' + cnf.bcbio_final_dir)
    info()
    info('*' * 70)

    if opts.qualimap and 'QualiMap' not in cnf.steps:
        cnf.steps.append('QualiMap')

    check_system_resources(cnf, required=['qsub'])

    cnf.work_dir = join(cnf.bcbio_final_dir, pardir, 'work', 'post_processing')
    if not isdir(cnf.work_dir):
        safe_mkdir(cnf.work_dir)

    load_bcbio_cnf(cnf)
    # if cnf.vcf_suf:
    #     vcf_sufs = cnf['vcf_suf'].split(',')
    # else:
    #     vcf_sufs = 'mutect'

    run_on_bcbio_final_dir(cnf, cnf.bcbio_final_dir, cnf.bcbio_cnf)


def load_bcbio_cnf(cnf):
    bcbio_config_dirpath = join(cnf.bcbio_final_dir, pardir, 'config')
    yaml_files = [join(bcbio_config_dirpath, fname) for fname in listdir(bcbio_config_dirpath) if fname.endswith('.yaml')]
    if len(yaml_files) > 1:
        critical('More than one YAML file in config directory: ' + ' '.join(yaml_files))
    if len(yaml_files) == 0:
        critical('No YAML file in config directory.')
    cnf.bcbio_cnf = load_yaml_config(yaml_files[0])


if __name__ == '__main__':
    main()









