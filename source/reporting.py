from collections import OrderedDict, defaultdict
from itertools import repeat, izip, chain
from os.path import join, relpath
from json import load, dump, JSONEncoder, dumps
from source import verify_file
import datetime

from source.bcbio_structure import VariantCaller, BCBioSample
from source.logger import critical, info, err
from source.html_reporting.html_saver import write_html_report


def to_dict_by_name(objects):
    return OrderedDict((o.name, o) for o in objects)


class Record:
    def __init__(
            self,
            metric=None,
            value=None,
            meta=None,
            html_fpath=None,

            num = None,
            cell_contents = None,
            frac_width = None,
            right_shift = None,
            color = 'white',
            text_color = 'black'):  #TODO: get rid of those

        self.metric = metric
        self.set_value(value)
        self.meta = meta or dict()
        self.html_fpath = html_fpath

        self.num = None
        self.cell_contents = None
        self.frac_width = None
        self.right_shift = None
        self.color = 'white'
        self.text_color = 'black'

    def get_value(self):
        return self.value

    def set_value(self, value):
        if value is None:
            pass
        else:
            if isinstance(value, basestring):
                try:
                    value = int(value)
                except ValueError:
                    try:
                        value = float(value)
                    except ValueError:
                        pass
        self.value = value

    def del_value(self):
        del self.value

    value = property(get_value, set_value, del_value, 'value')

    def format(self, human_readable=True):
        return self.metric.format(self.value, human_readable=human_readable)

    def format_html(self):
        return self.metric.format(self.value, human_readable=True, is_html=True)

    def __repr__(self):
        return self.metric.name + ' ' + self.metric.format(self.value, human_readable=True)

    @staticmethod
    def load(data):
        data['metric'] = Metric.load(data['metric'])
        return Record(**data)


class Metric:
    def __init__(
            self,
            name=None,
            short_name=None,
            description=None,
            quality='More is better',  # "More is better", "Less is better", "Equal"
            unit='',
            common=False,

            values = None,
            min = None,
            max = None,
            med = None,
            low_outer_fence = None,
            low_inner_fence = None,
            top_inner_fence = None,
            top_outer_fence = None,
            all_values_equal = False):  #TODO: get read of those
        self.name = name
        self.short_name = short_name
        self.description = description
        self.quality = quality
        self.common = common
        self.unit = unit

        self.values = []
        self.min = None
        self.max = None
        self.med = None
        self.low_outer_fence = None
        self.low_inner_fence = None
        self.top_inner_fence = None
        self.top_outer_fence = None
        self.all_values_equal = False

    def get_display_name(self):
        return self.short_name if self.short_name is not None else self.name

    def format(self, value, human_readable=True, is_html=False):
        return Metric.format_value(value, unit=self.unit, human_readable=human_readable, is_html=is_html)

    @staticmethod
    def format_value(value, unit='', human_readable=True, is_html=False):
        if value is None:
            return '.'

        if unit and is_html:
            unit = '<span class=\'rhs\'>&nbsp;</span>' + unit

        if isinstance(value, basestring):
            if human_readable:
                return '{value}{unit}'.format(**locals())
            else:
                return value

        elif isinstance(value, int):
            if value == 0:
                return '0'
            if human_readable:
                if value <= 9999:
                    return str(value)
                else:
                    return '{value:,}{unit}'.format(**locals())
            else:
                return str(value)

        elif isinstance(value, float):
            if value == 0.0:
                return '0'
            if human_readable:
                if unit == '%':
                    value *= 100
                presision = 2
                for i in range(10, 2, -1):
                    if value < 1./(10**i):
                        presision = i + 1
                v = '{value:.{presision}f}'.format(**locals())
                if is_html:
                    v = v.replace('.', '<span class=\'hs\'></span>')
                v += unit
                return v
            else:
                return str(value)

        if isinstance(value, list):
            return ','.join(Metric.format_value(v, unit, human_readable, is_html) for v in list)

        return '.'

    def __repr__(self):
        return self.name

    @staticmethod
    def load(data):
        return Metric(**data)


# noinspection PyClassHasNoInit
class Report:
    def flatten(self, sections=None, human_readable=True):
        raise NotImplementedError()

    def save_txt(self, output_fpath, sections=None):
        fpath = write_txt_rows(self.flatten(sections, human_readable=True), output_fpath)
        self.txt_fpath = fpath
        return fpath

    def save_json(self, output_fpath, sections=None):
        self.dump(output_fpath)
        self.json_fpath = output_fpath
        return output_fpath

    def save_tsv(self, output_fpath, sections=None):
        fpath = write_tsv_rows(self.flatten(sections, human_readable=False), output_fpath)
        self.tsv_fpath = fpath
        return fpath

    def save_html(self, output_fpath, caption='', type_=None):
        # class Encoder(JSONEncoder):
        #     def default(self, o):
        #         if isinstance(o, (VariantCaller, BCBioSample)):
        #             return o.for_json()
        #         return o.__dict__

        # json = dumps(dict(
        #     date=datetime.datetime.now().strftime('%d %B %Y, %A, %H:%M:%S'),
        #     report=self,
        #     type_=type_,
        # ), separators=(',', ':'), cls=Encoder)
        fpath = write_html_report(self, type_, output_fpath, caption)
        self.html_fpath = fpath
        return fpath

    def dump(self, fpath):
        with open(fpath, 'w') as f:
            dump(self, f, default=lambda o: o.__dict__, indent=4)

    @staticmethod
    def find_record(records, metric_name):
        try:
            rec = next(r for r in records if r.metric.name == metric_name)
        except StopIteration:
            return None  # if no record for the metric
        else:
            return rec

    def get_common_records(self):
        raise NotImplementedError()


class SampleReport(Report):
    def __init__(self, sample=None, html_fpath=None,
                 records=None, metric_storage=None,
                 report_name='', plots=None, json_fpath=None,
                 display_name=None, caller_tag=None, project_tag=None,
                 **kwargs):
        self.sample = sample
        self.html_fpath = html_fpath
        self.records = records or []
        self.metric_storage = metric_storage
        self.report_name = report_name
        self.plots = plots or []  # TODO: make real JS plots, not just included PNG
        self.json_fpath = json_fpath

        self.display_name = display_name or (self.sample if isinstance(self.sample, basestring) else self.sample.name)
        self.caller_tag = None
        self.set_caller_tag(caller_tag)
        self.project_tag = None
        self.set_project_tag(project_tag)

        self.display_name = self.get_display_name()  # filled in only in save_html()!!!
        # self.display_name = sample if isinstance(sample, basestring) else sample.name
        # if self.sample.project_tag:
        #     self.display_name = self.sample.project_tag + ' ' + self.display_name

    # def set_display_name(self, name):
    #     self.display_name = name
    #     return self

    def get_common_records(self):
        common_records = []
        for record in self.records:
            if record.metric.common:
                common_records.append(record)
        return common_records

    def set_project_tag(self, tag):
        if not self.project_tag and tag:
            self.display_name = '[' + tag + ']' + ' ' + self.display_name
            self.project_tag = tag

    def set_caller_tag(self, tag):
        if not self.caller_tag and tag:
            self.display_name = self.display_name + ' ' + tag
            self.project_tag = tag

    def get_display_name(self):
        return self.display_name
    #     name = self.sample if isinstance(self.sample, basestring) else self.sample.name
    #     if self.project_tag:
    #         name = self.sample.project_tag + ' ' + name
    #     if self.caller_tag:
    #         name = name + ' ' + self.caller_tag
    #     return name

    def add_record(self, metric_name, value, meta=None):
        metric = self.metric_storage.get_metric(metric_name.strip())
        assert metric, metric_name
        rec = Record(metric, value, meta)
        self.records.append(rec)
        info(metric_name + ': ' + rec.format(human_readable=True))
        return rec

    def flatten(self, sections=None, human_readable=True):
        rows = []

        for m in self.metric_storage.general_section.metrics:
            r = Report.find_record(self.records, m.name)
            if r:
                if human_readable:
                    rows.append(['## ' + m.name + ': ' + r.format(human_readable=True)])
                else:
                    rows.append(['##' + m.name + '=' + r.format(human_readable=False)])

        rows.append(['Sample', self.get_display_name()])
        for metric in self.metric_storage.get_metrics(sections, skip_general_section=True):
            row = [metric.name]
            rec = Report.find_record(self.records, metric.name)
            row.append(rec.format(human_readable=human_readable) if rec is not None else '.')
            rows.append(row)
        return rows

    def save_html(self, output_fpath, caption='', type_=None):
        return Report.save_html(self, output_fpath, caption=caption, type_='SampleReport')

    def __repr__(self):
        return self.get_display_name() + (', ' + self.report_name if self.report_name else '')

    @staticmethod
    def load(data, sample=None, bcbio_structure=None):
        data['sample'] = sample or BCBioSample.load(data['sample'], bcbio_structure)
        data['records'] = [Record.load(d) for d in data['records']]
        data['metric_storage'] = MetricStorage.load(data['metric_storage'])

        return SampleReport(**data)


class PerRegionSampleReport(SampleReport):
    class Region:
        def __init__(self, parent_report):
            self.__parent_report = parent_report
            self.records = []

        def add_record(self, metric_name, value, meta=None):
            metric = self.__parent_report.metric_storage.get_metric(metric_name.strip())
            assert metric, metric_name
            rec = Record(metric, value, meta)
            self.records.append(rec)
            return rec

    def __init__(self, *args, **kwargs):
        SampleReport.__init__(self, *args, **kwargs)
        self.records = []
        self.regions = []

    def add_region(self):
        region = PerRegionSampleReport.Region(self)
        self.regions.append(region)
        return region

    def flatten(self, sections=None, human_readable=True):
        rows = []

        for m in self.metric_storage.general_section.metrics:
            rec = Report.find_record(self.records, m.name)
            if rec:
                if human_readable:
                    rows.append(['## ' + m.name + ': ' + rec.format(human_readable=True)])
                else:
                    rows.append(['##' + m.name + '=' + rec.format(human_readable=False)])

        header_row = []
        ms = self.metric_storage.get_metrics(sections, skip_general_section=True)
        for i, m in enumerate(ms):
            header_row.append(('#' if i == 0 else '') + m.get_display_name())
        rows.append(header_row)

        for reg in self.regions:
            row = []
            for m in self.metric_storage.get_metrics(sections, skip_general_section=True):
                rec = Report.find_record(reg.records, m.name)
                if rec:
                    row.append(rec.format(human_readable=human_readable))
                else:
                    pass

            rows.append(row)

        # next_list_value = next((r.value for r in self.records if isinstance(r.value, list)), None)
        # if next_list_value:
        #     for i in range(len(next_list_value)):
        #         header_row = []
        #         for m in self.metric_storage.get_metrics(sections):
        #             try:
        #                 r = next((r for r in self.records if r.metric.name == m.name), None)
        #                 if r:
        #                     if m.name in self.metric_storage.general_section.metrics_by_name:
        #                         val = r.value
        #                     elif not r.value:
        #                         val = None
        #                     else:
        #                         val = r.value[i]
        #                     header_row.append(r.metric.format(val))
        #             except StopIteration:
        #                 header_row.append('-')  # if no record for the metric
        #         rows.append(header_row)

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
    def __init__(self, name='', sample_reports=None, metric_storage=None, general_records=None, base_fname=None):
        self.name = name
        self.sample_reports = sample_reports or []
        self.metric_storage = metric_storage
        self._general_records = general_records
        self.base_fname = base_fname
        self.plots = []

        if metric_storage:
            for sample_report in sample_reports:
                sample_report.metric_storage = metric_storage

        elif sample_reports and sample_reports[0].metric_storage:
            self.metric_storage = sample_reports[0].metric_storage
            for sample_report in sample_reports:
                sample_report.metric_storage = metric_storage

    def get_common_records(self):
        common_records = []

        if self._general_records:
            common_records.extend(self._general_records)

        if self.sample_reports:
            sample_report = self.sample_reports[0]
            for record in sample_report.records:
                if record.metric and record.metric.common:  #TODO: why can record.metric be None?
                    common_records.append(record)
        return common_records

    def flatten(self, sections=None, human_readable=True):
        if len(self.sample_reports) == 0:
            err('No sample reports found: summary will not be produced.')
            return []

        rows = []

        some_rep = self.sample_reports[0]
        for m in self.metric_storage.general_section.metrics:
            rec = Report.find_record(some_rep.records, m.name)
            if rec:
                if human_readable:
                    rows.append(['## ' + m.name + ': ' + rec.format(human_readable=True)])
                else:
                    rows.append(['##' + m.name + '=' + rec.format(human_readable=False)])

        rows.append(['Sample'] + [rep.get_display_name() for rep in self.sample_reports])

        for metric in self.metric_storage.get_metrics(sections, skip_general_section=True):
            row = [metric.name]
            for sr in self.sample_reports:
                rec = Report.find_record(sr.records, metric.name)
                row.append(rec.format(human_readable=human_readable) if rec is not None else '.')
            rows.append(row)
        return rows

    @staticmethod
    def construct_from_sample_report_jsons(samples, bcbio_structure, output_dirpath,
            jsons_by_sample, htmls_by_sample):
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

    def save_into_files(self, base_path, caption, sections=None):
        return \
            self.save_txt(base_path + '.txt', sections), \
            self.save_tsv(base_path + '.tsv', sections), \
            self.save_html(base_path + '.html', caption)

    def save_html(self, output_fpath, caption='', type_=None, display_name=None):
        if len(self.sample_reports) == 0:
            err('No sample reports found: HTML summary will not be made.')
            return None

        return Report.save_html(self, output_fpath, caption=caption, type_='FullReport')

    def __repr__(self):
        return self.name


def parse_tsv(tsv_fpath):
    values = []

    with open(tsv_fpath, 'r') as f:
        for line in f:
            values.append(line.split('\t'))

    if not values:
        critical('Data not found in ' + tsv_fpath)

    return values


class ReportSection:
    def __init__(self, name='', title='', metrics=None, **kwargs):
        self.name = name
        self.title = title
        self.metrics = metrics or []
        self.metrics_by_name = dict((m.name, m) for m in metrics)

    def add_metric(self, metric):
        self.metrics.append(metric)
        self.metrics_by_name[metric.name] = metric

    def get_metrics(self):
        return self.metrics

    def get_metric(self, metric_name):
        return next((m for m in self.get_metrics() if m.name == metric_name), None)

    def remove_metric(self, metric_name):
        metric = self.metrics_by_name[metric_name]
        self.metrics.remove(metric)
        del self.metrics_by_name[metric_name]

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

    def remove_metric(self, metric_name):
        for section in self.sections:
            section.remove_metric(metric_name)

    def get_metrics(self, sections=None, skip_general_section=False):
        metrics = []
        for section in ((sections or self.sections) + ([] if skip_general_section else [self.general_section])):
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


def load_records(json_fpath):
    with open(json_fpath) as f:
        data = load(f, object_pairs_hook=OrderedDict)
        return [Record.load(d) for d in data['records']]


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


def write_txt_rows(rows, output_fpath):
    if not rows:
        return None

    # output_fpath = join(output_dirpath, base_fname + '.txt')

    col_widths = repeat(0)
    for row in rows:
        if not row[0].startswith('##'):
            col_widths = [max(len(v), w) for v, w in izip(row, col_widths)]

    with open(output_fpath, 'w') as out:
        for row in rows:
            if row[0].startswith('##'):
                out.write(row[0])
            else:
                for val, w in izip(row, col_widths):
                    out.write(val + (' ' * (w - len(val) + 2)))
            out.write('\n')

    return output_fpath


def write_tsv_rows(rows, output_fpath):
    if not rows:
        return None

    with open(output_fpath, 'w') as out:
        for row in rows:
            out.write('\t'.join([val for val in row]) + '\n')

    return output_fpath


# def parse_value(string):
#     val = string.replace(' ', '').replace(',', '')
#
#     num_chars = []
#     unit_chars = []
#
#     i = 0
#     while i < len(val) and (val[i].isdigit() or val[i] == '.'):
#         num_chars += val[i]
#         i += 1
#     while i < len(val):
#         unit_chars += val[i]
#         i += 1
#
#     val_num = ''.join(num_chars)
#     val_unit = ''.join(unit_chars)
#
#     try:
#         val = int(val_num)
#     except ValueError:
#         try:
#             val = float(val_num)
#         except ValueError:
#             val = val_num
#
#     return val
