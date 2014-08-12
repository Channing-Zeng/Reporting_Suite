#!/usr/bin/env python

from __future__ import print_function
import sys
from source.file_utils import verify_dir

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from os.path import abspath, dirname, realpath, join, pardir
from site import addsitedir
source_dir = abspath(dirname(realpath(__file__)))
addsitedir(join(source_dir, 'ext_modules'))

from optparse import OptionParser
from source.config import Defaults, Config
from source.logger import step_greetings, info
from source.main import check_keys, check_inputs, set_up_dirs
from source.reporting import read_sample_names, get_per_sample_fpaths_for_bcbio_final_dir, \
    summarize, write_summary_reports, Metric, Record


def main():
    info(' '.join(sys.argv))
    info()

    description = 'This script generates project-level summaries based on per-sample ngsCAT reports.'

    parser = OptionParser(description=description)
    parser.add_option('-d', dest='bcbio_final_dir', help='Path to bcbio-nextgen final directory (default is pwd)')
    parser.add_option('-s', dest='samples', help='List of samples (default is samples.txt in bcbio final directory)')
    parser.add_option('-n', dest='base_name', default='qualimap',
                      help='Name of Qualimap directory inside sample folder. (default is qualimap)')
    parser.add_option('-o', '--output_dir', dest='output_dir', metavar='DIR',
                      help='output directory (or directory name in case of bcbio final dir)')

    parser.add_option('-v', dest='verbose', action='store_true', help='Verbose')
    parser.add_option('-t', dest='threads', type='int', help='Number of threads for each process')
    parser.add_option('-w', dest='overwrite', action='store_true', help='Overwrite existing results')

    parser.add_option('--runner', dest='qsub_runner',
                      help='Bash script that takes command line as the 1st argument. This script will be submitted '
                           'to GRID. Default: ' + Defaults.qsub_runner)
    parser.add_option('--work-dir', dest='work_dir', metavar='DIR')

    parser.add_option('--sys-cnf', '--sys-info', '--sys-cfg', dest='sys_cnf', default=Defaults.sys_cnf,
                      help='System configuration yaml with paths to external tools and genome resources (see default one %s)' % Defaults.sys_cnf)
    parser.add_option('--run-cnf', '--run-info', '--run-cfg', dest='run_cnf', default=Defaults.run_cnf,
                      help='Run configuration yaml (see default one %s)' % Defaults.run_cnf)

    (opts, args) = parser.parse_args()
    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)

    cnf.name = cnf['name'] or 'QualiMap'
    set_up_dirs(cnf)

    if not check_keys(cnf, ['bcbio_final_dir']):
        parser.print_help()
        sys.exit(1)

    cnf.bcbio_final_dir = verify_dir(cnf.bcbio_final_dir)
    if not cnf.bcbio_final_dir:
        sys.exit(1)

    if not cnf.samples:
        cnf.samples = join(cnf.bcbio_final_dir, 'samples.txt')

    info('BCBio "final" dir: ' + cnf.bcbio_final_dir + ' (set with -d)')
    info('Samples: ' + cnf.samples + ' (set with -s)')

    if not check_keys(cnf, ['samples']):
        parser.print_help()
        sys.exit(1)

    if 'qsub_runner' in cnf:
        cnf.qsub_runner = join(cnf.sys_cnf, pardir, cnf.qsub_runner)

    if not check_inputs(cnf, file_keys=['samples', 'qsub_runner'], dir_keys=['bcbio_final_dir']):
        sys.exit(1)

    sample_names = read_sample_names(cnf['samples'])

    sum_report_fpaths = summary_reports(cnf, sample_names)

    info()
    info('*' * 70)
    info('Summary:')
    for fpath in sum_report_fpaths:
        if fpath:
            info('  ' + fpath)


def summary_reports(cnf, sample_names):
    step_greetings('QualiMap statistics for all samples')

    sample_html_reports, sample_names = get_per_sample_fpaths_for_bcbio_final_dir(
        cnf['bcbio_final_dir'], sample_names, cnf['base_name'], 'qualimapReport.html', raw_ending=True)

    sum_report = summarize(sample_names, sample_html_reports, _parse_qualimap_sample_report)

    final_summary_report_fpaths = write_summary_reports(
        cnf['output_dir'], cnf['work_dir'], sum_report, 'qualimap', 'QualiMap statistics')

    return final_summary_report_fpaths


def __to_dict(metrics):
    return {m.name: m for m in metrics}

METRICS = __to_dict([
    Metric('Number of reads',                               'Reads',                            'Total number of reads'),
    Metric('Mapped reads',                                  'Mapped',                           'Number of mapped reads'),
    Metric('Unmapped reads',                                'Unmapped',                         'Number of unmapped reads',               quality='Less is better'),
    Metric('Paired reads',                                  'Paired',                           'Total number of paired reads'),
    Metric('Mapped reads, only first in pair',              'Mapped, only first',               'Number of mapped reads, only first in pair'),
    Metric('Mapped reads, only second in pair',             'Mapped, only second',              'Number of mapped reads, only second in pair'),
    Metric('Mapped reads, both in pair',                    'Mapped, both',                     'Number of mapped reads, both in pair'),
    Metric('Mapped reads, singletons',                      'Mapped, singletons',               'Number of mapped reads, singletons'),

    # TODO: split into 3 metrics
    Metric('Read min/max/mean length',                      'Read min/max/mean length',         'PLACEHOLDER for three separate metrics'),
    Metric('Read min length',                               'Read min length',                  'Read min length'),
    Metric('Read max length',                               'Read max length',                  'Read max length'),
    Metric('Read mean length',                              'Read mean length',                 'Read mean length'),

    Metric('Clipped reads',                                 'Clipped reads',                    'Number of clipped reads',                quality='Less is better'),
    Metric('Duplication rate',                              'Duplication',                      'Duplication rate'),
    Metric('Mapped reads (on target)',                      'Mapped (on target)',               'Number of mapped reads inside of regions'),
    Metric('Mapped reads, only first in pair (on target)',  'Mapped, only first (on target)',   'Number of mapped reads inside of regions, only first in pair'),
    Metric('Mapped reads, only second in pair (on target)', 'Mapped, only second (on target)',  'Number of mapped reads inside of regions, only second in pair'),
    Metric('Mapped reads, both in pair (on target)',        'Mapped, both (on target)',         'Number of mapped reads inside of regions, both in pair'),
    Metric('Mapped reads, singletons (on target)',          'Mapped, singletons (on target)',   'Number of mapped reads inside of regions, singletons'),

    Metric('Coverage Mean',                                 'Cov. mean',                        'Coverage mean, inside of regions'),
    Metric('Coverage Standard Deviation',                   'Cov. std. dev.',                   'Coverage std. dev., inside of regions',  quality='Less is better'),
    Metric('Mean Mapping Quality',                          'Mean mapping quality',             'Mean mapping quality, inside of regions'),
    Metric('Total reads with indels',                       'Indels',                           'Total reads with indels, inside of regions'),
    Metric('Insertions',                                    'Insertions',                       'Insertions, inside of regions'),
    Metric('Deletions',                                     'Deletions',                        'Deletions, inside of regions'),
    Metric('Homopolymer indels',                            'Homopolymer indels',               'Percentage of homopolymer indels, inside of regions')
])

ALLOWED_UNITS = ['%']


def _parse_qualimap_sample_report(report_fpath):
    records = []

    def __get_value(line):
        ## examples:
        # <td class=column1>Paired reads</td>
        # <td class=column2>80,244 / 99.89%</td>
        crop_left = line.split('>')
        if len(crop_left) < 2:
            return None
        crop_right = crop_left[1].split('<')
        return crop_right[0].strip()

    def __fill_record(metric_name, line):
        val = __get_value(line)
        val = val.replace(' ', '').replace(',', '')

        if metric_name == 'Read min/max/mean length':  # special case
            for metric_infix, value in zip(['min', 'max', 'mean'], val.split('/')):
                record = Record(METRICS['Read ' + metric_infix + ' length'])
                record.value = value
                records.append(record)
            return

        record = Record(METRICS[metric_name])
        num_chars = []
        unit_chars = []
        i = 0
        while i < len(val) and (val[i].isdigit() or val[i] in ['.']):
            num_chars += val[i]
            i += 1
        while i < len(val):
            unit_chars += val[i]
            i += 1
        val_num = ''.join(num_chars)
        val_unit = ''.join(unit_chars)

        if val_unit and val_unit in ALLOWED_UNITS:
            record.metric.unit = val_unit
        try:
            val = int(val_num)
        except ValueError:
            try:
                val = float(val_num)
            except ValueError:  # it is a string
                val = val_num + val_unit
        record.value = val
        if val_unit.startswith('/'):  # for values like "80,220 / 99.86%"
            record.meta = val_unit[1:]
        records.append(record)

    sections = {'start':            'Summary',
                'on target':        'Globals (inside of regions)',
                'coverage':         'Coverage (inside of regions)',
                'mapping quality':  'Mapping Quality',
                'finish':           'Input data and parameters'}
    on_target_stats_suffix = ' (on target)'
    coverage_stats_prefix = 'Coverage '
    with open(report_fpath) as f:
        cur_section = None
        cur_metric_name = None
        for line in f:
            for name, pattern in sections.items():
                if line.find(pattern) != -1:
                    cur_section = name
                    break
            if cur_section is None:
                continue
            if cur_section == 'finish':
                break

            if line.find('class=column1') != -1:
                value = __get_value(line)
                if cur_section == 'on target':
                    value += on_target_stats_suffix
                elif cur_section == 'coverage':
                    value = coverage_stats_prefix + value
                if value in METRICS.keys():
                    cur_metric_name = value
            if cur_metric_name and line.find('class=column2') != -1:
                __fill_record(cur_metric_name, line)
                cur_metric_name = None
    return records


if __name__ == '__main__':
    main()

