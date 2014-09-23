#!/usr/bin/env python
import sys
from source.bcbio_structure import Sample
from source.file_utils import adjust_path
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from os.path import abspath, dirname, realpath, pardir, join
from site import addsitedir
source_dir = abspath(dirname(realpath(__file__)))
addsitedir(join(source_dir, 'ext_modules'))

import shutil
from source.main import read_opts_and_cnfs, check_system_resources, load_genome_resources
from source.config import Defaults
from source.targetcov.cov import run_targetcov_reports
from source.runner import run_one
from source.utils import info


def main(args):
    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--bam'], dict(
                dest='bam',
                help='path to the BAM file')
             ),
            (['--bed', '--capture', '--amplicons'], dict(
                dest='bed',
                help='capture panel/amplicons')
             ),
            (['--caller-names'], dict(
                dest='caller_names',
                help='names of variant callers used to create vcfs provided by --vcf')
             ),
            (['--vcfs'], dict(
                dest='vcfs',
                help='filteted variants in VCF, comma-separate, must correspond to caller-names')
             ),
            (['--padding'], dict(
                dest='padding',
                help='integer indicating the number of bases to extend each target region up and down-stream. '
                     'Default is ' + str(Defaults.coverage_reports['padding']),
                type='int')
             ),
            (['--reports'], dict(
                dest='report_types',
                metavar=Defaults.coverage_reports['report_types'],
                help='Comma-separated report names.')
             ),
            (['--depth-thresholds'], dict(
                dest='depth_thresholds',
                metavar='A,B,C',
                help='Default: ' + ','.join(map(str,
                      Defaults.coverage_reports['depth_thresholds']))),
             ),
        ],
        required_keys=['bam', 'bed'],
        file_keys=['bam', 'bed'],
        key_for_sample_name='bam')

    check_system_resources(
        cnf,
        required=['samtools', 'bedtools'],
        optional=[])

    load_genome_resources(
        cnf,
        required=['seq', 'exons'],
        optional=['chr_lengths', 'genes'])

    if cnf['report_types']:
        cnf['coverage_reports']['report_types'] = cnf['report_types']
        del cnf['report_types']

    info('Using alignement ' + cnf['bam'])
    info('Using amplicons/capture panel ' + cnf['bed'])

    cnf.vcfs = map(adjust_path, cnf.vcfs.split(',') if cnf.vcfs else [])
    cnf.caller_names = cnf.caller_names.split(',') if cnf.caller_names else []
    cnf.vcfs_by_callername = zip(cnf.caller_names, cnf.vcfs)

    run_one(cnf, process_one, finalize_one)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


def process_one(cnf):
    sample = Sample(cnf.name, cnf.bam, cnf.bed)

    return run_targetcov_reports(cnf, sample, cnf.vcfs_by_callername)


def finalize_one(cnf, summary_report_txt_path, summary_report_json_path, summary_report_html_path,
                 gene_report_fpath, abnormal_regions_reports):
    if summary_report_txt_path:
        info('Summary report: ' + summary_report_txt_path)
    if gene_report_fpath:
        info('Exons coverage report:')
        info('   ' + gene_report_fpath)
    if abnormal_regions_reports:
        info('Abnormal region reports:')
        for rep in abnormal_regions_reports:
            info('  ' + rep)


if __name__ == '__main__':
    main(sys.argv)