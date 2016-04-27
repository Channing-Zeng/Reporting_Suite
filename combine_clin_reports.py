#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc

from collections import defaultdict
import os
import sys
from os.path import join
from optparse import OptionParser

from source.bcbio.bcbio_structure import BCBioStructure, process_post_bcbio_args
from source.clinical_reporting.combine_reports import run_combine_clinical_reports
from source.config import defaults
from source.logger import info, critical, err
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug, set_up_log
from source.file_utils import safe_mkdir, adjust_path


def main():
    info(' '.join(sys.argv))
    info()
    description = 'This script makes clinical reports based on multiple bcbio projects.'

    parser = OptionParser(description=description)
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser, threads=1)

    parser.add_option('--email', dest='email', help='E-mail address to send notifications on errors and finished jobs.')
    parser.add_option('--jira', dest='jira', help='JIRA case path')
    parser.add_option('--bed', '--capture', '--amplicons', dest='bed', help='BED file to run targetSeq and Seq2C analysis on.')
    parser.add_option('--exons', '--exome', dest='exons', help='Exons BED file to make targetSeq exon/amplicon regions reports.')
    parser.add_option('-o', dest='output_dir', help='Output directory for report combining.')

    cnf, bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths, tags, is_wgs_in_bcbio, is_rnaseq \
        = process_post_bcbio_args(parser)
    is_wgs = cnf.is_wgs = cnf.is_wgs or is_wgs_in_bcbio

    if len(bcbio_project_dirpaths) < 2:
        critical('Usage: ' + __file__ + ' project_bcbio_path [project_bcbio_path] [-o output_dir]')

    info()
    info('*' * 70)
    bcbio_structures = []
    for bcbio_project_dirpath, bcbio_cnf, final_dirpath in zip(
            bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths):
        bs = BCBioStructure(cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath,
                            is_wgs=is_wgs, is_rnaseq=is_rnaseq)
        bcbio_structures.append(bs)

    if not cnf.project_name:
        cnf.project_name = bcbio_structures[0].project_name.replace('_WGS', '').replace('WGS', '')

    if cnf.output_dir is None:
        cnf.output_dir = join(os.getcwd(), cnf.project_name)
    safe_mkdir(cnf.output_dir)

    cnf.log_dir = join(cnf.output_dir, 'log')
    info('log_dirpath: ' + cnf.log_dir)
    safe_mkdir(cnf.log_dir)
    set_up_log(cnf, 'combine_clin_reports', cnf.project_name, cnf.output_dir)

    cnf.work_dir = cnf.work_dir or adjust_path(join(cnf.output_dir, 'work'))
    safe_mkdir(cnf.work_dir)

    # shared_sample_names = set(s.name for bs in bcbio_structures for s in bs.samples)
    # if not shared_sample_names:
    #    sample_names = [bs.project_name + ': ' + ', '.join(s.name for s in bs.samples) for bs in bcbio_structures]
    #    critical('Not shared samples in projects.\n' + '\n'.join(sample_names))

    # info('Shared samples: ' + ', '.join(shared_sample_names))

    info('')
    info('*' * 70)
    run_combine_clinical_reports(cnf, bcbio_structures)


if __name__ == '__main__':
    main()
