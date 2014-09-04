#!/usr/bin/env python
from genericpath import isfile
import sys
from os import getcwd
from source.summary import process_post_bcbio_args, add_post_bcbio_args

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

from source.bcbio_runner import BCBioRunner
from source.file_utils import safe_mkdir, adjust_path, remove_quotes, verify_dir, verify_file
from source.config import Defaults, Config, load_yaml_config
from source.logger import info, critical
from source.main import check_system_resources, check_inputs, check_keys, load_genome_resources
from source.bcbio_structure import BCBioStructure, load_bcbio_cnf


def main():
    info(' '.join(sys.argv))
    info()
    description = 'This script runs reporting suite on the bcbio final directory.'

    parser = OptionParser(description=description)
    add_post_bcbio_args(parser)

    parser.add_option('--qualimap', dest='qualimap', action='store_true', default=Defaults.qualimap, help='Run QualiMap in the end')
    parser.add_option('--load-mongo', '--mongo-loader', dest='load_mongo', action='store_true', default=Defaults.load_mongo, help='Load to Mongo DB')
    cnf = process_post_bcbio_args(parser)

    if cnf.qualimap and 'QualiMap' not in cnf.steps: cnf.steps.append('QualiMap')
    if cnf.load_mongo and 'MongoLoader' not in cnf.steps: cnf.steps.append('MongoLoader')

    check_system_resources(cnf, required=['qsub'])

    load_genome_resources(cnf, required=['seq'])

    load_bcbio_cnf(cnf)
    info()
    info('*' * 70)

    bcbio_structure = BCBioStructure(cnf, cnf.bcbio_final_dir, cnf.bcbio_cnf)
    bcbio_runner = BCBioRunner(cnf, bcbio_structure, cnf.bcbio_cnf)
    bcbio_runner.post_jobs()


if __name__ == '__main__':
    main()









