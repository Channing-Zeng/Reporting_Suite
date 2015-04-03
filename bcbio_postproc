#!/usr/bin/env python

import sub_scripts.__check_python_version  # checking for python version and adding site dirs inside
import sys
from optparse import OptionParser
from source.bcbio_runner import BCBioRunner
from source.config import defaults
from source.logger import info
from source.bcbio_structure import BCBioStructure, process_post_bcbio_args
from source.prepare_args_and_cnf import add_post_bcbio_args, check_system_resources
from source.variants.vcf_processing import get_trasncripts_fpath


def main():
    info(' '.join(sys.argv))
    info()
    description = 'This script runs reporting suite on the bcbio final directory.'

    parser = OptionParser(description=description)
    add_post_bcbio_args(parser)

    parser.add_option('--load-mongo', '--mongo-loader', dest='load_mongo', action='store_true', default=defaults['load_mongo'], help='Load to Mongo DB')
    parser.add_option('--datahub-path', dest='datahub_path', help='DataHub directory path to upload final MAFs and CNV (can be remote).')
    parser.add_option('--email', dest='email', help='E-mail address to send notifications on errors and finished jobs.')
    parser.add_option('--reannotate', dest='reannotate', action='store_true', default=False, help='re-annotate BED file with gene names')
    parser.add_option('--count-dups', dest='count_dups', action='store_true', default=False, help='count duplicates in coverage metrics')

    cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath = process_post_bcbio_args(parser)

    if cnf.min_freq is not None:
        info('Min freq for filtering is %f' % cnf.min_freq)

    if cnf.steps and cnf.load_mongo and 'MongoLoader' not in cnf.steps:
        cnf.steps.append('MongoLoader')

    check_system_resources(cnf, required=['qsub'], optional='transcripts_fpath')

    info()
    info('*' * 70)

    bcbio_structure = BCBioStructure(cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath)

    bcbio_runner = BCBioRunner(cnf, bcbio_structure, cnf.bcbio_cnf)
    bcbio_runner.post_jobs()


if __name__ == '__main__':
    main()









