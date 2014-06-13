#!/usr/bin/python

import sys

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


def main():
    required_keys = ['bam', 'bed']
    optional_keys = ['extra_bam', 'extend', 'saturation', 'depthlist', 'feature']

    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--bam'], 'align.bam', {
                'dest': 'bam',
                'help': 'Path to the bam file. Required!'}),

            (['--extra_bam'], 'align2.bam', {
                'dest': 'extra_bam',
                'help': 'Additional bam file to compare.'}),

            (['--bed'], 'target_regions.bed', {
                'dest': 'bed',
                'help': 'Path to the bed file containing the target regions. Required!'}),

            (['--extend_target'], '<int>', {
                'dest': 'extend',
                'help': 'Integer indicating the number of bases to extend each target region up and down-stream'}),

            (['--saturation'], '{y,n}', {
                'dest': 'saturation',
                'help': 'Y/n to indicate whether saturation curve should be calculated',
                'default': 'n'}),

            (['--depth_list'], '<float,float,..>', {
                'dest': 'depthlist',
                'help': 'Will only be used in case --saturation is "y". Comma separated list of real numbers '
                        '(do not leave spaces between) indicating the number of millions of reads to simulate '
                        'for the saturation curve. E.g.: 1,5,10 would indicate 1*10^6, 5*10^6 and 10*10^6.',
                'default': 'auto'}),

            (['--one_feature'], '<str>', {
                'dest': 'feature',
                'help': "Use this option if just one of the graphs/statistics should be calculated. "
                        "String indicating one of the following features: "
                        "{%s}" % (', '.join(config.availablefeatures))}),
        ],
        key_for_sample_name='bam')

    check_system_resources(cnf,
        required=['samtools', 'bedtools'],
        optional=[])

    load_genome_resources(cnf,
        required=['chr_lengths'],
        optional=[])

    check_inputs(cnf,
        required_keys=['bam', 'bed'],
        file_keys=['bam', 'bed', 'extra_bam'])

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

    run_one(cnf, required_keys, optional_keys, process_one, finalize_one)


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

    if 'feature' in cnf:
        cnf['feature'] = cnf['feature'].lower()
    else:
        cnf['feature'] = None

    if cnf['feature'] and cnf['feature'] not in config.availablefeatures:
        critical('%s is not available. Please, check that the selected '
                 'feature is one of the following: %s'
                 % (cnf['feature'], ', '.join(config.availablefeatures)))

    if 'extend' in cnf:
        cnf['extend'] = int(cnf['extend'])
    else:
        cnf['extend'] = None


def process_one(cnf):
    bam_file = cnf['bam']

    bams = [cnf['bam']]
    if 'extra_bam' in cnf:
        bams.append(cnf['extra_bam'])

    filtered_bams = [bam_file.filter_unmapped_reads(cnf, bam) for bam in bams]

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
