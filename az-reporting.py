#!/usr/bin/env python
import sys
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from os.path import dirname, realpath, join
from site import addsitedir
source_dir = dirname(realpath(__file__))
addsitedir(join(source_dir, 'ext_modules'))

from optparse import OptionParser
from source.bcbio_runner import BCBioRunner
from source.config import Defaults
from source.logger import info
from source.main import check_system_resources, load_genome_resources
from source.bcbio_structure import BCBioStructure, load_bcbio_cnf
from source.prepare_args_and_cnf import add_post_bcbio_args, process_post_bcbio_args
from source.variants.vcf_processing import get_trasncripts_fpath


def main():
    info(' '.join(sys.argv))
    info()
    description = 'This script runs reporting suite on the bcbio final directory.'

    parser = OptionParser(description=description)
    add_post_bcbio_args(parser)

    parser.add_option('--load-mongo', '--mongo-loader', dest='load_mongo', action='store_true', default=Defaults.load_mongo, help='Load to Mongo DB')
    parser.add_option('--datahub-path', dest='datahub_path', help='DataHub directory path to upload final MAFs and CNV (can be remote).')
    # parser.add_option('--email', dest='email', help='Email to send notifications on errors and finished jobs.')
    
    cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath = process_post_bcbio_args(parser)

    if cnf.steps and cnf.load_mongo and 'MongoLoader' not in cnf.steps:
        cnf.steps.append('MongoLoader')

    check_system_resources(cnf, required=['qsub'], optional='transcripts_fpath')

    info()
    info('*' * 70)

    bcbio_structure = BCBioStructure(cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath)

    get_trasncripts_fpath(cnf)

    bcbio_runner = BCBioRunner(cnf, bcbio_structure, cnf.bcbio_cnf)
    bcbio_runner.post_jobs()


if __name__ == '__main__':
    main()









