from os.path import join
from source.reporting import summarize, write_summary_reports, Metric, Record
from source.logger import step_greetings, info
from source.bcbio_structure import BCBioStructure


def summary_reports(cnf, bcbio_structure):
    step_greetings('QualiMap statistics for all samples')

    html_by_sample = bcbio_structure.get_qualimap_report_fpaths_by_sample()
    sum_report = summarize(cnf, html_by_sample, _parse_qualimap_sample_report, '')

    final_summary_report_fpaths = write_summary_reports(
        cnf.output_dir,
        cnf.work_dir,
        sum_report,
        join(cnf.output_dir, BCBioStructure.qualimap_name),
        'QualiMap statistics')

    info()
    info('*' * 70)
    info('Summary:')
    for fpath in final_summary_report_fpaths:
        if fpath:
            info('  ' + fpath)

    return final_summary_report_fpaths


METRICS = Metric.to_dict([
    Metric('Number of reads',                               'Reads',                            'Total number of reads'),
    Metric('Mapped reads',                                  'Mapped',                           'Number of mapped reads'),
    Metric('Unmapped reads',                                'Unmapped',                         'Number of unmapped reads',               quality='Less is better'),
    Metric('Paired reads',                                  'Paired',                           'Total number of paired reads'),
    Metric('Mapped reads, only first in pair',              'Mapped, only first',               'Number of mapped reads, only first in pair'),
    Metric('Mapped reads, only second in pair',             'Mapped, only second',              'Number of mapped reads, only second in pair'),
    Metric('Mapped reads, both in pair',                    'Mapped, both',                     'Number of mapped reads, both in pair'),
    Metric('Mapped reads, singletons',                      'Mapped, singletons',               'Number of mapped reads, singletons'),

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


def _parse_qualimap_sample_report(cnf, report_fpath):
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