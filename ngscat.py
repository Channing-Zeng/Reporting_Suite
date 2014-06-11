#!/usr/bin/python

import sys

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

if not 'matplotlib' in sys.modules:
    import matplotlib

    matplotlib.use('Agg')  # non-GUI backend

from source.ngscat import config, bed_file, bam_file, ngscat_main

from source.main import read_opts_and_cnfs, check_system_resources, load_genome_resources
from source.utils import verify_file, critical, step_greetings, info
from source.runner import run_one
from os.path import join, expanduser, isdir, abspath, dirname


def verify_bam(fpath):
    if not verify_file(fpath):
        return False
    if not fpath.endswith('.bam'):
        sys.stderr.write(fpath + ' must have .bam extension. '
                                 'Please, make sure that the bam file is appropriately formatted.\n')
        return False
    textchars = ''.join(map(chr, [7, 8, 9, 10, 12, 13, 27] + range(0x20, 0x100)))
    is_binary_string = lambda bytes: bool(bytes.translate(None, textchars))
    if not is_binary_string(open(fpath).read(3)):
        sys.stderr.write(fpath + ' must be a binary file. '
                                 'Please, make sure that the bam file is appropriately formatted.\n')
        return False
    return True


def verify_bed(fpath):
    if not verify_file(fpath):
        return False
    err = bed_file.bed_file(fpath).checkformat()
    if err:
        sys.stderr.write('ERROR: incorrect bed file format (' + fpath + '): ' + err + '\n')
        return False
    return True


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
                'help': 'Additional bam file'}),

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
        required_keys=required_keys,
        optional_keys=optional_keys,
        key_for_sample_name='bam')

    check_system_resources(cnf, ['samtools', 'bedtools'])
    load_genome_resources(cnf)
    #load_genome_resources(cnf, ['chr_lengths']) #TODO: check whether it needed at all!
    #chr_len_fpath = cnf.get('chr_lengths') or cnf['genome'].get('chr_lengths')
    #if not chr_len_fpath:
    #    critical('Specify chromosome lengths for the genome'
    #             ' in system info or in run info.')
    #config.CHR_LENGTHS = chr_len_fpath

    if 'coverage_reports' not in cnf:
        critical('No coverage_reports section in the report, cannot run NGScat.')

    read_samples_info(cnf)
    config.cnf = cnf

    # TODO: cnf['parallel'] = False?
    cnf['parallel'] = False
    cnf['name'] = dirname(cnf['bam'])  # TODO: remove

    run_one(cnf, required_keys, optional_keys, process_one, finalize_one)


def read_samples_info(cnf):
    info('')
    info('Processing input details...')

    if not verify_bam(cnf['bam']) or not verify_bed(cnf['bed']):
        exit(1)

    if cnf['extra_bam'] and not verify_bam(cnf['extra_bam']):
        exit(1)

    if cnf['saturation'] not in ['y', 'n']:
        critical('Incorrect value for --saturation parameter. Please indicate "y" or "n".')

    try:
        cnf['depthlist'] = cnf['depthlist'] if cnf['depthlist'] == 'auto' \
            else map(float, cnf['depthlist'].split(','))
    except ValueError:
        critical('Invalid values for --depth_list option. Please, provide a comma separated '
                 'list of values without leaving spaces, e.g.: 1.0,2.5,10.0')

    cnf['feature'] = None if 'feature' not in cnf else cnf['feature'].lower()
    if cnf['feature'] is not None and cnf['feature'] not in config.availablefeatures:
        critical('%s is not available. Please, check that the selected '
                 'feature is one of the following: %s'
                 % (cnf['feature'], ', '.join(config.availablefeatures)))

    cnf['extend'] = None if 'extend' not in cnf is None else int(cnf['extend'])
    info('writing to output dir ' + cnf['output_dir'])


def process_one(cnf, *inputs):
    filtered_bams = []
    bams = [cnf['bam']]
    if 'extra_bam' in cnf:
        bams.append(cnf['extra_bam'])
    for bam in bams:
        filtered_bams.append(bam_file.filter_unmapped_reads(cnf, bam))

    report_fpath = ngscat_main.ngscat(cnf, filtered_bams, cnf['bed'], cnf['output_dir'], cnf['genome'].get('seq'),
                                      cnf['saturation'] == 'y', int(cnf.get('threads', '1')), cnf['extend'],
                                      cnf['depthlist'], cnf['coverage_reports']['depth_thresholds'],
                                      cnf['feature'])
    return [report_fpath]


def finalize_one(cnf, report_fpath):
    if report_fpath:
        info('ngsCAT report: ' + report_fpath)


if __name__ == '__main__':
    main()
