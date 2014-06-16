#!/usr/bin/python

import sys
from source.config import Defaults

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

if not 'matplotlib' in sys.modules:
    import matplotlib
    matplotlib.use('Agg')  # non-GUI backend

from os.path import join, expanduser, isdir, abspath, dirname

from source.main import read_opts_and_cnfs, check_system_resources, load_genome_resources, check_inputs
from source.utils import verify_file, critical, info
from source.runner import run_one
from source.ngscat import config, ngscat_main
from source.ngscat.bam_file import filter_unmapped_reads


def main():
    required_keys = ['bam', 'bed']
    optional_keys = ['extra_bam', 'extend', 'saturation', 'depthlist', 'feature']

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
                default=Defaults.coverage_reports['padding'])
             ),
            (['--saturation'], dict(
                dest='saturation',
                help='Y/n to indicate whether saturation curve should be calculated',
                metavar='{y,n}',
                default=Defaults.ngscat['saturation'])
             ),
            (['--depth-list'], dict(
                dest='depthlist',
                help='will only be used in case --saturation is "y". Comma separated list of real numbers '
                     '(do not leave spaces between) indicating the number of millions of reads to simulate '
                     'for the saturation curve. E.g.: 1,5,10 would indicate 1*10^6, 5*10^6 and 10*10^6.',
                metavar='<float,float,..>',
                default=Defaults.ngscat['depthlist'])
             ),
            (['--one_feature'], dict(
                dest='feature',
                metavar='<str>',
                help='Use this option if just one of the graphs/statistics should be calculated. '
                     'String indicating one of the following features: '
                     '{%s}' % (', '.join(Defaults.ngscat['availablefeatures'])))),
        ],
        required_keys=['bam', 'bed'],
        file_keys=['bam', 'bed', 'extra_bam'],
        key_for_sample_name='bam')

    check_system_resources(cnf,
        required=['samtools', 'bedtools'],
        optional=[])

    load_genome_resources(cnf,
        required=['chr_lengths'],
        optional=[])

    #load_genome_resources(cnf, ['chr_lengths']) #TODO: check whether it needed at all!
    #chr_len_fpath = cnf.get('chr_lengths') or cnf['genome'].get('chr_lengths')
    #if not chr_len_fpath:
    #    critical('Specify chromosome lengths for the genome'
    #             ' in system info or in run info.')
    #config.CHR_LENGTHS = chr_len_fpath

    if 'coverage_reports' not in cnf:
        critical('No coverage_reports section in the report, cannot run NGScat.')

    process_parameters(cnf)

    # TODO: cnf['parallel'] = False?
    cnf['parallel'] = False
    cnf['name'] = dirname(cnf['bam'])  # TODO: remove

    run_one(cnf, process_one, finalize_one)


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

    if cnf.feature and cnf.feature not in Defaults.ngscat['availablefeatures']:
        critical('%s is not available. Please, check that the selected '
                 'feature is one of the following: %s'
                 % (cnf['feature'], ', '.join(Defaults.ngscat['availablefeatures'])))

    if 'extend' in cnf:
        cnf['extend'] = int(cnf['extend'])
    else:
        cnf['extend'] = None


def process_one(cnf):
    bams = [cnf['bam']]
    if 'extra_bam' in cnf:
        bams.append(cnf['extra_bam'])

    filtered_bams = [filter_unmapped_reads(cnf, bam) for bam in bams]

    report_fpath = ngscat_main.ngscat(
        cnf, filtered_bams, cnf['bed'], cnf['output_dir'], cnf['genome'].get('seq'),
        cnf['saturation'] == 'y', int(cnf.get('threads', '1')), cnf['extend'],
        cnf['depthlist'], cnf['coverage_reports']['depth_thresholds'],
        cnf['feature'])

    return [report_fpath]


def finalize_one(cnf, report_fpath):
    if report_fpath:
        info('ngsCAT report: ' + report_fpath)


if __name__ == '__main__':
    main()
