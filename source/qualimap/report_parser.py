from source.reporting import Metric, Record, MetricStorage, ReportSection


metric_storage = MetricStorage(
    sections=[
        # ReportSection('basic_metrics', 'General', [
            # Metric('Number of reads',                               'Reads ngscat',                       'Total number of reads'),
            # Metric('Mapped reads',                                  'Mapped ngscat',                      'Number of mapped reads'),
            # Metric('Unmapped reads',                                'Unmapped ngscat',                    'Number of unmapped reads',               quality='Less is better'),
        # ]),

        # ReportSection('on_off_metrics', 'ON/OFF target', [
            # Metric('Mapped reads, only first in pair',              'Mapped, 1st',                 'Number of mapped reads, only first in pair'),
            # Metric('Mapped reads, only second in pair',             'Mapped, 2nd',                 'Number of mapped reads, only second in pair'),
            # Metric('Mapped reads, both in pair',                    'Mapped, both',                'Number of mapped reads, both in pair'),
            # Metric('Mapped reads, singletons',                      'Mapped, single',              'Number of mapped reads, singletons'),
            #
            # Metric('Mapped reads (on target)',                      'Mapped (on trg)',             'Number of mapped reads inside of regions'),
            # Metric('Mapped reads, only first in pair (on target)',  'Mapped, 1st (on trg)',        'Number of mapped reads inside of regions, only first in pair'),
            # Metric('Mapped reads, only second in pair (on target)', 'Mapped, 2nd (on trg)',        'Number of mapped reads inside of regions, only second in pair'),
            # Metric('Mapped reads, both in pair (on target)',        'Mapped, both (on trg)',       'Number of mapped reads inside of regions, both in pair'),
            # Metric('Mapped reads, singletons (on target)',          'Mapped, single (on trg)',     'Number of mapped reads inside of regions, singletons')
        # ]),

        # ReportSection('depth_metrics', 'Target coverage depth', [
            # Metric('Coverage Mean',                                 'Cov. mean',                   'Coverage mean, inside of regions'),
            # Metric('Coverage Standard Deviation',                   'Cov. std. dev.',              'Coverage std. dev., inside of regions',  quality='Less is better')
        # ]),

        ReportSection('reads', 'Reads', [
            # Metric('Paired reads',                                  'Paired',                      'Total number of paired reads'),
            Metric('Read min length',                               'Min len',                     'Read min length'),
            Metric('Read max length',                               'Max len',                     'Read max length'),
            Metric('Read mean length',                              'Ave len',                     'Read mean length'),
        ]),

        ReportSection('qualimap', 'Qualimap metrics, inside the regions (unless it is a WGS study)', [
            Metric('Mean Mapping Quality',                          'Mean MQ',                     'Mean mapping quality, inside of regions'),
            Metric('Mismatches',                                    'Mismatches',                  'Mismatches, inside of regions', quality='Less is better'),  # added in Qualimap v.2.0
            Metric('Insertions',                                    'Insertions',                  'Insertions, inside of regions', quality='Less is better'),
            Metric('Deletions',                                     'Deletions',                   'Deletions, inside of regions', quality='Less is better'),
            Metric('Homopolymer indels',                            'Homopolymer indels',          'Percentage of homopolymer indels, inside of regions', quality='Less is better'),
            # Metric('Duplication rate',                              'Duplication rate',            'Duplication rate (inside of regions)')
        ])
    ]
)


ALLOWED_UNITS = ['%']


def parse_qualimap_sample_report(report_fpath):
    records = []

    def __get_td_tag_contents(line):
        ## examples:
        # <td class=column1>Paired reads</td>
        # <td class=column2>80,244 / 99.89%</td>
        crop_left = line.split('>')
        if len(crop_left) < 2:
            return None
        crop_right = crop_left[1].split('<')
        return crop_right[0].strip()

    def __fill_record(metric_name, line):
        val = __get_td_tag_contents(line)
        val = val.replace(' ', '').replace(',', '')

        if metric_name == 'Read min/max/mean length':  # special case
            for metric_infix, value in zip(['min', 'max', 'mean'], val.split('/')):
                metric = metric_storage.get_metric('Read ' + metric_infix + ' length')
                rec = Record(metric, value)
                records.append(rec)
        else:
            metric = metric_storage.get_metric(metric_name)
            if not metric:
                return

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
                metric.unit = val_unit
            try:
                val = int(val_num)
                if val_unit == '%':
                    val = float(val) / 100
            except ValueError:
                try:
                    val = float(val_num)
                    if val_unit == '%':
                        val /= 100
                except ValueError:  # it is a string
                    val = val_num + val_unit

            if not metric_storage.get_metric(metric_name):
                return None

            rec = Record(metric_storage.get_metric(metric_name), val)
            if val_unit.startswith('/'):  # for values like "80,220 / 99.86%"
                rec.meta = val_unit[1:]
            records.append(rec)

    sections = {'start':            'Summary',
                'on target':        'Globals (inside of regions)',
                'coverage':         'Coverage (inside of regions)',
                'mapping quality':  'Mapping Quality',
                'finish':           'Coverage across reference'}  # plots are starting from this line
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
                cur_metric_name = __get_td_tag_contents(line)
                if cur_section == 'on target':
                    cur_metric_name += on_target_stats_suffix
                elif cur_section == 'coverage':
                    cur_metric_name = coverage_stats_prefix + cur_metric_name
                elif not metric_storage.get_metric(cur_metric_name):  # special case for Duplication rate and Clipped reads (Qualimap v.1 and v.2 difference)
                    if metric_storage.get_metric(cur_metric_name + on_target_stats_suffix):  # extra 'on target' metrics
                        cur_metric_name += on_target_stats_suffix

            if cur_metric_name and line.find('class=column2') != -1:
                __fill_record(cur_metric_name, line)
                cur_metric_name = None
    return records