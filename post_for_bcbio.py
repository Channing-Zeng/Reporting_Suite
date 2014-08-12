#!/usr/bin/env python
import sys
from source.bcbio_runner import BCBioRunner
from source.bcbio_structure import BCBioStructure, load_bcbio_cnf

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from os.path import abspath, dirname, realpath, pardir, join, basename
from site import addsitedir
source_dir = abspath(dirname(realpath(__file__)))
addsitedir(join(source_dir, 'ext_modules'))

from optparse import OptionParser
from os.path import join, pardir, isdir, basename, splitext, abspath
from os import listdir

from source.file_utils import safe_mkdir, adjust_path, remove_quotes
from source.config import Defaults, Config, load_yaml_config
from source.logger import info, critical
from source.main import check_system_resources, check_inputs, check_keys, load_genome_resources


def main():
    description = 'This script runs reporting suite on the bcbio final directory.'

    parser = OptionParser(description=description)
    parser.add_option('-d', dest='bcbio_final_dir', help='Path to bcbio-nextgen final directory (default is pwd)')
    parser.add_option('--qualimap', dest='qualimap', action='store_true', default=Defaults.qualimap, help='Run QualiMap in the end')
    parser.add_option('-v', dest='verbose', action='store_true', help='Verbose output')
    # parser.add_option('-t', dest='threads', type='int', help='Number of threads for each process')
    parser.add_option('-w', dest='overwrite', action='store_true', help='Overwrite existing results')
    parser.add_option('--reuse', dest='overwrite', action='store_false', help='Reuse intermediate results in work directory for subroutines')

    parser.add_option('--sys-cnf', '--sys-info', '--sys-cfg', dest='sys_cnf', default=Defaults.sys_cnf, help='system configuration yaml with paths to external tools and genome resources (see default one %s)' % Defaults.sys_cnf)
    parser.add_option('--run-cnf', '--run-info', '--run-cfg', dest='run_cnf', default=Defaults.run_cnf, help='run configuration yaml (see default one %s)' % Defaults.run_cnf)

    (opts, args) = parser.parse_args()
    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)
    if not opts.bcbio_final_dir and len(args) > 0:
        cnf.bcbio_final_dir = args[0]
    else:
        critical('Usage: ./post_for_bcbio.py <final_dir>')

    if not check_keys(cnf, ['bcbio_final_dir']):
        parser.print_help()
        sys.exit(1)

    if not check_inputs(cnf, dir_keys=['bcbio_final_dir']):
        sys.exit(1)

    if isdir(join(cnf.bcbio_final_dir, 'final')):
        cnf.bcbio_final_dir = join(cnf.bcbio_final_dir, 'final')

    if 'qsub_runner' in cnf:
        cnf.qsub_runner = remove_quotes(cnf.qsub_runner)
        cnf.qsub_runner = adjust_path(join(cnf.sys_cnf, pardir, cnf.qsub_runner))
    if not check_inputs(cnf, file_keys=['qsub_runner'], dir_keys=['bcbio_final_dir']):
        sys.exit(1)

    info('BCBio "final" dir: ' + cnf.bcbio_final_dir)

    if opts.qualimap and 'QualiMap' not in cnf.steps:
        cnf.steps.append('QualiMap')

    check_system_resources(cnf, required=['qsub'])

    load_genome_resources(cnf, required=['seq'])

    load_bcbio_cnf(cnf)

    info()
    info('*' * 70)

    bcbio_structure = BCBioStructure(cnf, cnf.bcbio_final_dir, cnf.bcbio_cnf)
    bcbio_runner = BCBioRunner(cnf, bcbio_structure)
    bcbio_runner.post_jobs()


if __name__ == '__main__':
    main()









