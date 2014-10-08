from collections import OrderedDict, defaultdict
from itertools import repeat, izip, chain
from os.path import join, relpath
from ext_modules.simplejson import load, dump, JSONEncoder, dumps
import datetime

from source.bcbio_structure import VariantCaller, Sample
from source.logger import critical, info, err
from source.html_reporting.html_saver import write_html_report


def to_dict_by_name(objects):
    return OrderedDict((o.name, o) for o in objects)


class Record:
    def __init__(self,
                 metric=None,
                 value=None,
                 meta=None,
                 html_fpath=None):
        self.metric = metric
        self.value = value
        self.meta = meta or dict()
        self.html_fpath = html_fpath

    def format(self):
        return self.metric.format(self.value)

    @staticmethod
    def load(data):
        data['metric'] = Metric.load(data['metric'])
        return Record(**data)


class Metric:
    def __init__(self,
                 name=None,
                 short_name=None,
                 description=None,
                 quality='More is better',  # "More is better", "Less is better", "Equal"
                 unit='',
                 common=False):
        self.name = name
        self.short_name = short_name
        self.description = description
        self.quality = quality
        self.common = common
        self.unit = unit

    def format(self, value):
        if value is None:
            return '-'

        name = self.name
        unit = self.unit

        if isinstance(value, basestring):
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    # assert False, 'Strange value ' + str(value)
                    return '{value}{unit}'.format(**locals())

        if isinstance(value, int):
            return '{value:,}{unit}'.format(**locals())

        if isinstance(value, float):
            presision = 2
            for i in range(10, 2, -1):
                if value < 1./(10**i):
                    presision = i + 1
            return '{value:.{presision}f}{unit}'.format(**locals())

        if isinstance(value, list):
            return ','.join(list)

        return '-'

    def __repr__(self):
        return self.name

    @staticmethod
    def load(data):
        return Metric(**data)


# noinspection PyClassHasNoInit
class Report:
    def flatten(self, sections=None):
        raise NotImplementedError()

    def save_txt(self, output_dirpath, base_fname, sections=None):
        return write_txt_rows(self.flatten(sections), output_dirpath, base_fname)

    def save_tsv(self, output_dirpath, base_fname, sections=None):
        return write_tsv_rows(self.flatten(sections), output_dirpath, base_fname)

    def save_html(self, output_dirpath, base_fname, caption='', type_=None):
        class Encoder(JSONEncoder):
            def default(self, o):
                if isinstance(o, (VariantCaller, Sample)):
                    return o.for_json()
                return o.__dict__

        json = dumps(dict(
            date=datetime.datetime.now().strftime('%d %B %Y, %A, %H:%M:%S'),
            data_outside_reports={},
            report=self,
            type_=None,
        ), separators=(',', ':'), cls=Encoder)
        return write_html_report(json, output_dirpath, base_fname, caption)

    @staticmethod
    def _append_value_to_row(sample_report, row, metric):
        try:
            row.append(next(
                r.metric.format(r.value)
                for r in sample_report.records
                if r.metric.name == metric.name))
        except StopIteration:
            row.append('-')  # if no record for the metric
        return row


class SampleReport(Report):
    def __init__(self, sample=None, html_fpath=None,
                 records=None, metric_storage=None,
                 report_name='', plots=None, json_fpath=None,
                 **kwargs):
        self.sample = sample
        self.html_fpath = html_fpath
        self.records = records or []
        self.metric_storage = metric_storage
        self.report_name = report_name
        self.plots = plots or []  # TODO: make real JS plots, not just included PNG
        self.json_fpath = json_fpath
        self.display_name = sample.name

    def set_display_name(self, name):
        self.display_name = name
        return self

    def add_record(self, metric_name, value, meta=None):
        metric = self.metric_storage.get_metric(metric_name.strip())
        assert metric, metric_name
        rec = Record(metric, value, meta)
        self.records.append(rec)
        info(metric_name + ': ' + rec.format())
        return rec

    def flatten(self, sections=None):
        rows = ['Sample', self.display_name]
        for metric in self.metric_storage.get_metrics(sections):
            row = [metric.name]
            Report._append_value_to_row(self, row, metric)
            rows.append(row)
        return rows

    def save_html(self, output_dirpath, base_fname, caption='', type_=None):
        return Report.save_html(self, output_dirpath, base_fname, caption=caption, type_='SampleReport')

    def __repr__(self):
        return self.display_name + (', ' + self.report_name if self.report_name else '')

    def dump(self, fpath):
        with open(fpath, 'w') as f:
            dump(self, f, default=lambda o: o.__dict__, indent=4)

    @staticmethod
    def load(data, sample=None, bcbio_structure=None):
        data['sample'] = sample or Sample.load(data['sample'], bcbio_structure)
        data['records'] = [Record.load(d) for d in data['records']]
        data['metric_storage'] = MetricStorage.load(data['metric_storage'])

        return SampleReport(**data)


class SquareSampleReport(SampleReport):
    def __init__(self, *args, **kwargs):
        SampleReport.__init__(self, *args, **kwargs)

    def flatten(self, sections=None):
        rows = []

        row = []
        for m in self.metric_storage.get_metrics(sections):
            if m.name in [r.metric.name for r in self.records]:
                row.append(m.name)
        rows.append(row)

        for i in range(len(next(r for r in self.records if isinstance(r.value, list)).value)):
            row = []
            for m in self.metric_storage.get_metrics(sections):
                try:
                    r = next((r for r in self.records if r.metric.name == m.name), None)
                    if r:
                        if m.name in self.metric_storage.general_section.metrics_by_name:
                            val = r.value
                        elif not r.value:
                            val = None
                        else:
                            val = r.value[i]
                        row.append(r.metric.format(val))
                except StopIteration:
                    row.append('-')  # if no record for the metric
            rows.append(row)
        return rows

    def save_html(self, output_dirpath, base_fname, caption='', type_=None):
        return None
        # sample_reports = []
        # fr = FullReport(self.report_name, sample_reports, self.metric_storage)

        # for i in range(len(next(r for r in self.records if isinstance(r.value, list)).value)):row = []
        #     row = []
        #     for m in self.metric_storage.get_metrics(self.metric_storage.sections):
        #         try:
        #             r = next((r for r in self.records if r.metric.name == m.name), None)
        #             if r:
        #                 if m.name in self.metric_storage.general_section.metrics_by_name:
        #                     val = r.value
        #                 elif not r.value:
        #                     val = None
        #                 else:
        #                     val = r.value[i]
        #                 row.append(r.metric.format(val))
        #         except StopIteration:
        #             row.append('-')  # if no record for the metric

            # records = []
            # sr = SampleReport(self.sample, self.html_fpath, records=None, metric_storage=None,
            #      report_name='', plots=None, json_fpath=None,
            #      )
            # return Report.save_html(self, output_dirpath, base_fname, caption=caption, type_='SquareSampleReport')

    # def add_record(self, metric_name, value, meta=None):
    #     raise NotImplementedError
    #
    # def add_row(self, row):
    #     # self.metric_storage.get_metric(metric_name.strip())
    #     assert metric, metric_name
    #     rec = Record(metric, value, meta)
    #     self.records.append(rec)
    #     info(metric_name + ': ' + rec.format())
    #     return rec


class FullReport(Report):
    def __init__(self, name='', sample_reports=None, metric_storage=None):
        self.name = name
        self.sample_reports = sample_reports or []
        self.metric_storage = metric_storage
        if metric_storage:
            for sample_report in sample_reports:
                sample_report.metric_storage = metric_storage

        elif sample_reports and sample_reports[0].metric_storage:
            self.metric_storage = sample_reports[0].metric_storage
            for sample_report in sample_reports:
                sample_report.metric_storage = metric_storage

    def get_common_records(self):
        common_records = list()
        if self.sample_reports:
            sample_report = self.sample_reports[0]
            for record in sample_report.records:
                if record.metric.common:
                    common_records.append(record)
        return common_records

    def flatten(self, sections=None):
        rows = [['Sample'] + [rep.display_name for rep in self.sample_reports]]

        if len(self.sample_reports) == 0:
            err('No sample reports found: summary will not be produced.')
            return []

        for metric in self.metric_storage.get_metrics(sections):
            row = [metric.name]
            for sr in self.sample_reports:
                Report._append_value_to_row(sr, row, metric)
            rows.append(row)
        return rows

    @staticmethod
    def construct_from_sample_report_jsons(samples, bcbio_structure, jsons_by_sample, htmls_by_sample, output_dirpath):
        full_report = FullReport()
        metric_storage = None
        for sample in samples:
            if sample.name in jsons_by_sample:
                with open(jsons_by_sample[sample.name]) as f:
                    data = load(f, object_pairs_hook=OrderedDict)
                    sample_report = SampleReport.load(data, sample, bcbio_structure)
                    sample_report.html_fpath = relpath(htmls_by_sample[sample.name], output_dirpath) \
                        if sample.name in htmls_by_sample else None
                    full_report.sample_reports.append(sample_report)
                    metric_storage = metric_storage or sample_report.metric_storage

        for sample_report in full_report.sample_reports:
            sample_report.metric_storage = metric_storage

        full_report.metric_storage = metric_storage

        return full_report

    def save_into_files(self, output_dirpath, base_fname, caption, sections=None):
        return \
            self.save_txt(output_dirpath, base_fname, sections), \
            self.save_tsv(output_dirpath, base_fname, sections), \
            self.save_html(output_dirpath, base_fname, caption)

    def save_html(self, output_dirpath, base_fname, caption='', type_=None):
        return Report.save_html(self, output_dirpath, base_fname, caption=caption, type_='FullReport')

    def __repr__(self):
        return self.name


class ReportSection:
    def __init__(self, name='', title='', metrics=None, **kwargs):
        self.name = name
        self.title = title
        self.metrics = metrics or []
        self.metrics_by_name = dict((m.name, m) for m in metrics)

    def add_metric(self, metric):
        self.metrics.append(metric)
        self.metrics_by_name[metric.name] = metric

    def __repr__(self):
        return self.name + ', ' + self.title

    @staticmethod
    def load(data):
        data['metrics'] = [Metric(**m) for m in data['metrics']]
        return ReportSection(**data)


class MetricStorage:
    def __init__(self,
                 metrics_list=None,
                 general_section=None,
                 sections=None,
                 sections_by_name=None,
                 **kwargs):
        self.sections_by_name = OrderedDict()
        self.sections = []
        self.general_section = general_section or ReportSection('general_section', '', [])
        self.general_section.name = self.general_section.name or 'general_section'

        if sections:
            self.sections = sections
            for section in sections:
                self.sections_by_name[section.name] = section

        elif sections_by_name:
            for name, section in sections_by_name.items():
                self.sections_by_name[name] = section
                self.sections.append(section)

        elif metrics_list:
            section = ReportSection('metrics_list', '', metrics_list)
            self.sections_by_name[''] = section
            self.sections.append(section)

    def add_metric(self, metric, section_name=''):
        self.sections_by_name[section_name].add_metric(metric)

    def get_metric(self, metric_name):
        return next((m for m in self.get_metrics() if m.name == metric_name), None)

    def get_metrics(self, sections=None):
        metrics = []
        for section in ((sections or self.sections) + [self.general_section]):
            if sections:
                pass
            for m in section.metrics:
                metrics.append(m)
        return metrics

    def __repr__(self, *args, **kwargs):
        return self.sections_by_name.__repr__(*args, **kwargs)

    @staticmethod
    def load(data):
        return MetricStorage(
            sections=[ReportSection.load(d) for d in data['sections']],
            general_section=ReportSection.load(data.get('general_section') or data.get('common_for_all_samples_section')))


def read_sample_names(sample_fpath):
    sample_names = []

    with open(sample_fpath) as f:
        for line in f:
            sample_name = line.strip()
            sample_names.append(sample_name)

    return sample_names


# def get_per_sample_fpaths_for_bcbio_final_dir(
#         bcbio_final_dir, sample_names, base_dir, ending, raw_ending=False):
#
#     single_report_fpaths = []
#
#     fixed_sample_names = []
#     for sample_name in sample_names:
#         report_name = ending if raw_ending else sample_name + ending
#         single_report_fpath = join(
#             bcbio_final_dir, sample_name, base_dir,
#             report_name)
#
#         info(single_report_fpath)
#
#         if not file_exists(single_report_fpath):
#             info('No ' + single_report_fpath + ', skipping.')
#             continue
#
#         if not verify_file(single_report_fpath):
#             critical(single_report_fpath + ' does not exist.')
#
#         single_report_fpaths.append(single_report_fpath)
#         fixed_sample_names.append(sample_name)
#
#     return single_report_fpaths, fixed_sample_names


# def summarize(cnf, get_sample_repot_fn, parse_report_fn, report_name):
#     """ Returns list of SampleReport objects:
#         [SampleReport(sample=Sample(name=), fpath=, records=[Record,...]),...]
#     """
#     return FullReport(
#         name=report_name,
#         sample_reports=[get_sample_repot_fn
#             SampleReport(sample_name, fpaths.html_fpath, parse_report_fn(cnf, fpaths.content_fpath))
#             for sample_name, fpaths in report_fpath_by_sample.items()])


def write_txt_rows(rows, output_dirpath, base_fname):
    if not rows:
        return None

    output_fpath = join(output_dirpath, base_fname + '.txt')

    col_widths = repeat(0)
    for row in rows:
        col_widths = [max(len(v), w) for v, w in izip(row, col_widths)]

    with open(output_fpath, 'w') as out:
        for row in rows:
            for val, w in izip(row, col_widths):
                out.write(val + (' ' * (w - len(val) + 2)))
            out.write('\n')

    return output_fpath


def write_tsv_rows(rows, output_dirpath, base_fname):
    if not rows:
        return None

    output_fpath = join(output_dirpath, base_fname + '.tsv')

    with open(output_fpath, 'w') as out:
        for row in rows:
            out.write('\t'.join([val for val in row]) + '\n')

    return output_fpath


def parse_tsv(tsv_fpath):
    report = []

    with open(tsv_fpath, 'r') as f:
        for line in f:
            report.append(line.split('\t'))

    if not report:
        critical('Data not found in ' + tsv_fpath)

    return report


def parse_value(string):
    val = string.replace(' ', '').replace(',', '')

    num_chars = []
    unit_chars = []

    i = 0
    while i < len(val) and (val[i].isdigit() or val[i] == '.'):
        num_chars += val[i]
        i += 1
    while i < len(val):
        unit_chars += val[i]
        i += 1

    val_num = ''.join(num_chars)
    val_unit = ''.join(unit_chars)

    try:
        val = int(val_num)
    except ValueError:
        try:
            val = float(val_num)
        except ValueError:
            val = val_num

    return val
