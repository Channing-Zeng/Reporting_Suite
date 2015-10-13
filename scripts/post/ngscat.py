#!/usr/bin/python
import __check_python_version

import sys
from source.bcbio_structure import BCBioStructure
from source.prepare_args_and_cnf import check_genome_resources, check_system_resources
from source.targetcov.bam_and_bed_utils import index_bam
from os.path import isfile
from source.logger import err

try:
    if 'matplotlib' not in sys.modules:
        import matplotlib
        matplotlib.use('Agg')  # non-GUI backend
except ImportError:
    err('Cannot import matplotlib')
    pass

import shutil
from os.path import dirname

from source.logger import info, critical, err
from source.main import read_opts_and_cnfs
from source.config import defaults
from source.runner import run_one
from source.ngscat import config, ngscat_main
from source.ngscat.bam_file import filter_unmapped_reads, process_pysam_bug


def main():
    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--bam'], dict(
                dest='bam',
                help='path to the BAM file')
             ),
            (['--extra_bam'], dict(
                dest='extra_bam',
                help='additional BAM file to compare')
             ),
            (['--bed', '--capture', '--amplicons'], dict(
                dest='bed',
                help='capture panel/amplicons')
             ),
            (['--padding'], dict(
                dest='padding',
                help='integer indicating the number of bases to extend each target region up and down-stream',
                type='int',
                default=defaults['coverage_reports']['padding'])
             ),
            (['--saturation'], dict(
                dest='saturation',
                help='Y/n to indicate whether saturation curve should be calculated',
                metavar='{y,n}',
                default=defaults['coverage_reports']['saturation'])
             ),
            (['--depth-list'], dict(
                dest='depthlist',
                help='will only be used in case --saturation is "y". Comma separated list of real numbers '
                     '(do not leave spaces between) indicating the number of millions of reads to simulate '
                     'for the saturation curve. E.g.: 1,5,10 would indicate 1*10^6, 5*10^6 and 10*10^6.',
                metavar='<float,float,..>',
                default=defaults['coverage_reports']['depthlist'])
             ),
            (['--one_feature'], dict(
                dest='feature',
                metavar='<str>',
                help='Use this option if just one of the graphs/statistics should be calculated. '
                     'String indicating one of the following features: '
                     '{%s}' % (', '.join(defaults['coverage_reports']['availablefeatures'])))),
        ],
        required_keys=['bam', 'bed'],
        file_keys=['bam', 'bed', 'extra_bam'],
        key_for_sample_name='bam',
        proc_name=BCBioStructure.ngscat_name)

    check_system_resources(cnf,
        required=['samtools', 'bedtools'],
        optional=[])

    check_genome_resources(cnf)

    if 'coverage_reports' not in cnf:
        critical('No coverage_reports section in the report, cannot run NGScat.')

    process_parameters(cnf)

    # TODO: cnf['parallel'] = False?
    cnf['parallel'] = False
    cnf['name'] = dirname(cnf['bam'])  # TODO: remove

    run_one(cnf, process_one, finalize_one)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


def process_parameters(cnf):
    config.cnf = cnf

    if cnf['saturation'] not in ['y', 'n']:
        critical('Incorrect value for --saturation parameter. Please indicate "y" or "n".')

    try:
        if cnf['depthlist'] != 'auto':
            cnf['depthlist'] = map(float, cnf['depthlist'].split(','))
    except ValueError:
        critical('Invalid values for --depth_list option. Please, provide a comma separated '
                 'list of values without leaving spaces, e.g.: 1.0,2.5,10.0')

    if cnf.feature:
        cnf['feature'] = cnf['feature'].lower()

    if cnf.feature and cnf.feature not in defaults['ngscat']['availablefeatures']:
        critical('%s is not available. Please, check that the selected '
                 'feature is one of the following: %s'
                 % (cnf['feature'], ', '.join(defaults['ngscat']['availablefeatures'])))

    if 'extend' in cnf:
        cnf['extend'] = int(cnf['extend'])
    else:
        cnf['extend'] = None


def process_one(cnf):
    bams = [cnf['bam']]
    if 'extra_bam' in cnf:
        bams.append(cnf['extra_bam'])

    for bam in bams:
        if not isfile(bam + '.bai'):
            info('Indexing bam ' + bam)
            index_bam(cnf, bam)

    filtered_bams = [process_pysam_bug(cnf, filter_unmapped_reads(cnf, bam)) for bam in bams]

    report_fpath = ngscat_main.ngscat(
        cnf, filtered_bams, cnf['bed'], cnf['output_dir'], cnf['genome'].get('seq'),
        cnf['saturation'] == 'y', int(cnf.get('threads') or '1'), cnf['extend'],
        cnf['depthlist'], cnf['coverage_reports']['depth_thresholds'],
        cnf['feature'])

    return [report_fpath]


def finalize_one(cnf, report_fpath):
    if report_fpath:
        info('ngsCAT report: ' + report_fpath)


if __name__ == '__main__':
    main()
