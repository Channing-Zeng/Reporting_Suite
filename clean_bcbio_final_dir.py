#!/usr/bin/env python
import sys

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from os.path import dirname, realpath, join, isdir, abspath
from site import addsitedir
source_dir = abspath(dirname(realpath(__file__)))
addsitedir(join(source_dir, 'ext_modules'))

from optparse import OptionParser

from source.config import Defaults, Config
from source.logger import info, critical
from source.main import check_inputs, check_keys
from source.bcbio_structure import BCBioStructure, load_bcbio_cnf


def main():
    # description = 'Clean reporting-postprocessing results.'
    #
    # parser = OptionParser(description=description)
    # parser.add_option('-d', dest='bcbio_final_dir', help='Path to bcbio-nextgen final directory (default is pwd)')
    # parser.add_option('-v', dest='verbose', action='store_true', help='Verbose output')
    # parser.add_option('--sys-cnf', '--sys-info', '--sys-cfg', dest='sys_cnf', default=Defaults.sys_cnf, help='system configuration yaml with paths to external tools and genome resources (see default one %s)' % Defaults.sys_cnf)
    # parser.add_option('--run-cnf', '--run-info', '--run-cfg', dest='run_cnf', default=Defaults.run_cnf, help='run configuration yaml (see default one %s)' % Defaults.run_cnf)
    #
    # (opts, args) = parser.parse_args()
    # cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)
    # if not opts.bcbio_final_dir and len(args) > 0:
    #     cnf.bcbio_final_dir = args[0]
    # else:
    #     critical('Usage: ' + __file__ + ' <final_dir>')
    #
    # if not check_keys(cnf, ['bcbio_final_dir']):
    #     parser.print_help()
    #     sys.exit(1)
    #
    # if not check_inputs(cnf, dir_keys=['bcbio_final_dir']):
    #     sys.exit(1)
    #
    # if isdir(join(cnf.bcbio_final_dir, 'final')):
    #     cnf.bcbio_final_dir = join(cnf.bcbio_final_dir, 'final')
    #
    # info('BCBio "final" dir: ' + cnf.bcbio_final_dir)
    #
    # load_bcbio_cnf(cnf)
    #
    # info()
    # info('*' * 70)

    bcbio_structure = BCBioStructure(cnf, cnf.bcbio_final_dir, cnf.bcbio_cnf)
    bcbio_structure.clean()


if __name__ == '__main__':
    main()









