#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc

from collections import defaultdict
import os
import sys
from os.path import join
from optparse import OptionParser

from source.bcbio.bcbio_structure import BCBioStructure, process_post_bcbio_args
from source.bcbio.bcbio_runner import BCBioRunner
from source.config import defaults
from source.logger import info, critical, err
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug, check_system_resources, set_up_log
from source.file_utils import safe_mkdir, adjust_path

from ngs_reporting.combine_reports import run_clinical_target2wgs


def main():
    info(' '.join(sys.argv))
    info()
    description = 'This script makes paired WGS-target clincial reports based on 2 bcbio projects.'

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

    if len(bcbio_project_dirpaths) < 2 or len(bcbio_project_dirpaths) > 2:
        critical('Usage: ' + __file__ + ' wgs_project_project_bcbio_path '
                 'targetseq_project_bcbio_path [-o output_dir]')

    info()
    info('*' * 70)
    bcbio_structures = []
    for bcbio_project_dirpath, bcbio_cnf, final_dirpath in zip(
            bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths):
        bs = BCBioStructure(cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath,
                            is_wgs=is_wgs, is_rnaseq=is_rnaseq)
        bcbio_structures.append(bs)

    trg_bs = next((bs for bs in bcbio_structures if bs.bed), None)
    wgs_bs = next((bs for bs in bcbio_structures if not bs.bed), None)
    if not trg_bs and not wgs_bs:
        critical('One of the projects must be targeted, and one must be WGS')
    if not trg_bs:
        critical('One of the projects must be targeted.')
    if not wgs_bs:
        critical('One of the projects must be WGS.')

    if not cnf.project_name:
        cnf.project_name = wgs_bs.project_name.replace('_WGS', '').replace('WGS', '')

    if cnf.output_dir is None:
        cnf.output_dir = join(os.getcwd(), cnf.project_name)
    safe_mkdir(cnf.output_dir)

    cnf.log_dir = join(cnf.output_dir, 'log')
    info('log_dirpath: ' + cnf.log_dir)
    safe_mkdir(cnf.log_dir)
    set_up_log(cnf, 'clinical_target2wgs', cnf.project_name, cnf.output_dir)

    cnf.work_dir = cnf.work_dir or adjust_path(join(cnf.output_dir, 'work'))
    safe_mkdir(cnf.work_dir)

    shared_sample_names = set(s.name for s in wgs_bs.samples) & set(s.name for s in trg_bs.samples)
    if not shared_sample_names:
        critical('Not shared samples in target and WGS projects.\n'
                 'Target: ' + ', '.join(s.name for s in trg_bs.samples) +
                 'WGS: ' + ', '.join(s.name for s in wgs_bs.samples))
    info('Shared samples: ' + ', '.join(shared_sample_names))

    info('')
    info('*' * 70)
    run_clinical_target2wgs(cnf.genome, wgs_bs, trg_bs, shared_sample_names, cnf.output_dir)


if __name__ == '__main__':
    main()
