# coding: utf-8

from collections import OrderedDict, defaultdict
from itertools import repeat, izip, chain
import os
from os.path import join, relpath, dirname, abspath, splitext, isdir, basename
from json import load, dump, JSONEncoder, dumps
from math import floor
import traceback
import datetime
from ext_modules.jsontemplate import jsontemplate
import shutil

from source.bcbio_structure import BCBioSample
from source.file_utils import file_transaction, file_exists, verify_file
from source.logger import critical, info, err, warn
from source.utils import mean


def get_real_path(path_in_html_saver):
    return join(dirname(abspath(__file__)), path_in_html_saver)


def to_dict_by_name(objects):
    return OrderedDict((o.name, o) for o in objects)


class Record:
    def __init__(
            self,
            metric=None,
            value=None,
            meta=None,
            html_fpath=None,
            parse=True,

            num = None,
            cell_contents = None,
            frac_width = None,
            right_shift = None,
            color = 'white',
            text_color = 'black'):  #TODO: get rid of those

        self.metric = metric
        if parse:
            self.set_value(value)
        else:
            self.value = value
        self.meta = meta or dict()
        self.html_fpath = html_fpath

        self.num = None
        self.cell_contents = None
        self.frac_width = None
        self.right_shift = None
        self.color = 'white'
        self.text_color = 'black'
        # self.color = lambda: self._color
        # self.text_color = lambda: self._text_color

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
        val = self.value
        # if self.html_fpath:
        #     if isinstance(self.html_fpath, dict):
        #         val = self.value + ': '
        #         val += ', '.join('<a href="' + v + '">' + k + '</a>'
        #                          for k, v in self.html_fpath.items())
        #     else:
        #         val = '<a href="' + self.html_fpath + '">' + self.value + '</a>'
        return self.metric.format(val, human_readable=True, is_html=True)

    def __repr__(self):
        return self.metric.name + ' ' + self.metric.format(self.value, human_readable=True)

    @staticmethod
    def load(data):
        data['metric'] = Metric.load(data['metric'])
        return Record(**data)


class Metric:
    skip_when_loading = [
        'numbers',
        'values',
        'min',
        'max',
        'med',
        'low_outer_fence',
        'low_inner_fence',
        'top_inner_fence',
        'top_outer_fence',
        'top_outer_fence',
        'all_values_equal'
    ]

    def __init__(
            self,
            name=None,
            short_name=None,
            description=None,
            quality='More is better',  # "More is better", "Less is better", "Equal"
            unit='',
            common=False,
            ok_threshold=None,
            bottom=None,
            is_hidden=False,
            with_heatmap=True,
            style='',
            class_='',
            align=None,  # TODO: legacy, remove
            width=None,
            max_width=None,
            min_width=None,

            numbers=None,
            values=None,
            min=None,
            max=None,
            med=None,
            low_outer_fence=None,
            low_inner_fence=None,
            top_inner_fence=None,
            top_outer_fence=None,
            all_values_equal=False):  # TODO: get read of those
        self.name = name
        self.short_name = short_name
        self.description = description
        self.quality = quality
        self.common = common
        self.unit = unit
        self.ok_threshold = ok_threshold
        self.bottom = bottom
        self.is_hidden = is_hidden
        self.with_heatmap = with_heatmap
        self.style = style
        self.class_ = class_
        self.max_width = max_width
        self.min_width = min_width

        self.numbers = []
        self.values = []
        self.min = None
        self.max = None
        self.med = med
        self.low_outer_fence = low_inner_fence
        self.low_inner_fence = low_inner_fence
        self.top_inner_fence = top_inner_fence
        self.top_outer_fence = top_outer_fence
        self.all_values_equal = False

    def get_display_name(self):
        return self.short_name if self.short_name is not None else self.name

    def format(self, value, human_readable=True, is_html=False):
        return Metric.format_value(value, unit=self.unit, human_readable=human_readable, is_html=is_html)

    @staticmethod
    def format_value(value, unit='', human_readable=True, is_html=False):
        if value is None:
            return '.'

        unit_str = unit
        if unit and is_html:
            unit_str = '<span class=\'rhs\'>&nbsp;</span>' + unit

        if isinstance(value, basestring):
            if human_readable:
                return '{value}{unit_str}'.format(**locals())
            else:
                return value

        elif isinstance(value, int):
            if value == 0:
                return '0'
            if human_readable:
                if value <= 9999:
                    return str(value)
                else:
                    v = '{value:,}{unit_str}'.format(**locals())
                    if is_html:
                        v = v.replace(',', '<span class=\'hs\'></span>')
                    return v
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
                return '{value:.{presision}f}{unit_str}'.format(**locals())
            else:
                return str(value)

        if isinstance(value, list):
            return ','.join(Metric.format_value(v, unit, human_readable, is_html) for v in list)

        return '.'

    def __repr__(self):
        return self.name

    @staticmethod
    def load(data):
        data = {k: v for k, v in data.items() if k not in Metric.skip_when_loading}
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

    def save_tsv(self, output_fpath, sections=None, human_readable=False):
        fpath = write_tsv_rows(self.flatten(sections, human_readable=human_readable), output_fpath)
        self.tsv_fpath = fpath
        return fpath

    def save_html(self, cnf, output_fpath, caption='', type_=None,
                  extra_js_fpaths=list(), extra_css_fpaths=list()):
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
        fpath = write_html_report(cnf, self, type_, output_fpath, caption=caption,
              extra_js_fpaths=extra_js_fpaths, extra_css_fpaths=extra_css_fpaths)
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
        except AttributeError:
            return None  # record with no metric
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
        self.records = [r for r in records if r.metric] if records else []
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

    def add_record(self, metric_name, value, meta=None, html_fpath=None, silent=False, parse=True):
        metric = self.metric_storage.find_metric(metric_name.strip())
        if not metric:
            err('Could not find metric ' + metric_name)
            return None
        rec = Record(metric, value, meta, html_fpath=html_fpath, parse=parse)
        self.records.append(rec)
        if not silent:
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

    def save_html(self, cnf, output_fpath, caption='', type_=None, extra_js_fpaths=list(), extra_css_fpaths=list()):
        return Report.save_html(self, cnf, output_fpath, caption=caption, type_='SampleReport',
                                extra_js_fpaths=list(), extra_css_fpaths=list())

    def __repr__(self):
        return self.get_display_name() + (', ' + self.report_name if self.report_name else '')

    @staticmethod
    def load(data, sample=None, bcbio_structure=None):
        data['sample'] = sample or BCBioSample.load(data['sample'], bcbio_structure)
        data['records'] = [Record.load(d) for d in data['records']]
        data['metric_storage'] = MetricStorage.load(data['metric_storage'])

        rep = SampleReport(**data)
        for rec in rep.records:
            rec.metric = rep.metric_storage.find_metric(rec.metric.name)
        return rep


class PerRegionSampleReport(SampleReport):
    class Region:
        def __init__(self, parent_report):
            self.__parent_report = parent_report
            self.records = []

        def add_record(self, metric_name, value, meta=None, parse=True):
            metric = self.__parent_report.metric_storage.find_metric(metric_name)
            assert metric, metric_name
            rec = Record(metric, value, meta, parse=parse)
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
        metrics = self.metric_storage.get_metrics(sections, skip_general_section=True)
        for i, m in enumerate(metrics):
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

    def save_html(self, cnf, output_fpath, caption='', type_=None,
                  extra_js_fpaths=list(), extra_css_fpaths=list()):
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
    def _sync_sections(dst_section, src_section):
        for metric in src_section.metrics:
            if not dst_section.find_metric(metric.name):
                dst_section.add_metric(metric)

    @staticmethod
    def construct_from_sample_report_jsons(samples, output_dirpath,
            jsons_by_sample, htmls_by_sample, bcbio_structure=None):
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
                    if metric_storage is None:
                        metric_storage = sample_report.metric_storage
                    else:  # Maximize metric storage
                        FullReport._sync_sections(metric_storage.general_section, sample_report.metric_storage.general_section)
                        for section in sample_report.metric_storage.sections:
                            if section.name not in metric_storage.sections_by_name:
                                metric_storage.add_section(section)
                            else:
                                FullReport._sync_sections(metric_storage.sections_by_name[section.name], section)

        for sample_report in full_report.sample_reports:
            sample_report.metric_storage = metric_storage
            for rec in sample_report.records:
                rec.metric = metric_storage.find_metric(rec.metric.name)
            sample_report.records = [r for r in sample_report.records if r.metric is not None]

        full_report.metric_storage = metric_storage

        return full_report

    def save_into_files(self, cnf, base_path, caption, sections=None):
        return \
            self.save_txt(base_path + '.txt', sections), \
            self.save_tsv(base_path + '.tsv', sections), \
            self.save_html(cnf, base_path + '.html', caption)

    def save_html(self, cnf, output_fpath, caption='', type_=None, display_name=None,
                  extra_js_fpaths=list(), extra_css_fpaths=list()):
        if len(self.sample_reports) == 0:
            err('No sample reports found: HTML summary will not be made.')
            return None

        return Report.save_html(self, cnf, output_fpath, caption=caption, type_='FullReport',
            extra_js_fpaths=extra_js_fpaths, extra_css_fpaths=extra_css_fpaths)

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

    def find_metric(self, metric_name):
        return self.metrics_by_name.get(metric_name, None)
        # return next((m for m in self.get_metrics() if m.name == metric_name), None)

    def remove_metric(self, metric_name):
        metric = self.metrics_by_name[metric_name]
        self.metrics.remove(metric)
        del self.metrics_by_name[metric_name]

    def __repr__(self):
        return self.name + ', ' + self.title

    @staticmethod
    def load(data):
        data['metrics'] = [Metric.load(m_data) for m_data in data['metrics']]
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

    def add_section(self, section):
        self.sections.append(section)
        self.sections_by_name[section.name] = section

    def add_metric(self, metric, section_name=''):
        self.sections_by_name[section_name].add_metric(metric)

    def find_metric(self, metric_name):
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


# HTML Saver

scripts_inserted = False

template_fpath = get_real_path(join('html_reporting', 'template.html'))
static_template_fpath = get_real_path(join('html_reporting', 'static_template.html'))

static_dirname = 'static'
static_dirpath = get_real_path(join('html_reporting', static_dirname))

aux_dirname = 'html_aux'

css_files = [
    'bootstrap/bootstrap.min.css',
    'common.css',
    'report.css',
    'table_sorter/style.css'
]
js_files = [
    'jquery-1.8.2.min.js',
    # 'flot/jquery.flot.min.js',
    # 'flot/jquery.flot.dashes.js',
    # 'scripts/hsvToRgb.js',
    'scripts/utils.js',
    # 'scripts/build_total_report.js',
    # 'scripts/build_report.js',
    # 'flot/excanvas.min.js',
    # 'ie_html5.js',
    'bootstrap/bootstrap.min.js',
    'bootstrap/bootstrap-tooltip-vlad.js',
    # 'bootstrap/bootstrap-tooltip-5px-lower.min.js',
    'dragtable.js',
    'table_sorter/tsort.js',
]
image_files = [
    # 'table_sorter/arrow_asc.png',
    # 'table_sorter/arrow_desc.png',
    # 'img/draggable.png',
]

def _build_report(report, type_):
    report_html = ''
    report_html += _build_common_records(report.get_common_records())

    for section in report.metric_storage.sections:
        column_names = [m.name for m in section.metrics]
        # columnOrder = (recoverOrderFromCookies section.name) or report.order or [0...columnNames.length]
        column_order = range(len(column_names))  # TODO: read column order

        report_html += _build_total_report(report, section, column_order)

    return report_html


def _build_common_records(general_records):
    for rec in general_records:
        _calc_record_cell_contents(rec)

    table = '<table cellspacing="0" class="common_table" id="common_table">'
    for rec in general_records:
        table += '\n<tr><td>'
        if rec.html_fpath:
            table += '<a href="' + rec.html_fpath + '>' + rec.cell_contents + '</a>'
        else:
            table += '<span class="metric_name">' + __get_metric_name_html(rec.metric, use_full_name=True) + ': </span>' + rec.cell_contents
        table += '</td></tr>'
    table += '\n</table>\n'
    table += '<div class="space_30px"></div>'
    return table


def __get_metric_name_html(metric, use_full_name=False):
    if metric.short_name and not use_full_name:
        metric_name = metric.short_name
        description = metric.description or metric.name
        return '<a class="metric_name" rel="tooltip" title="' + description + '">' + metric_name + '</a>'
    else:
        return metric.name


def make_cell_th(metric, pos=''):
    # return '<th class="' + metric.class_ + '" style="' + metric.style + '">' + metric.name + '</th>'
    html = ''

    if metric.is_hidden or not metric.values:
        return html

    style = metric.style
    if metric.max_width is not None:
        style += '; max-width: {w}px; width: {w}px;'.format(w=metric.max_width)
    if metric.min_width is not None:
        style += '; min-width: {w}px;'.format(w=metric.min_width)
    html += '\n<th style="' + style + '" class="second_through_last_col_headers_td"'

    if metric.numbers:
        sort_by = 'nosort' if metric.all_values_equal else 'numeric'
        direction = 'ascending' if metric.quality == 'Less is better' else 'descending'
        html += ' data-sortBy=' + sort_by + ' data-direction=' + direction + ' position=' + str(pos)

    # if metric.width is not None:
    #     html += ' width=' + str(metric.width)
    html += '><span class="metricName">' + __get_metric_name_html(metric) + '</span></th>'
    return html


def make_cell_td(rec, td_classes=''):
    html = ''

    if not rec.metric.values:
        return html

    if rec is None:
        html += "\n<td><span>-</span></td>"
        return html

    if not rec.metric:
        warn('rec.metric is None. (rec.value = ' + str(rec.value) + ')')
        return ''

    if rec.metric.is_hidden:
        return ''

    style = 'background-color: ' + rec.color + '; color: ' + rec.text_color + '; ' + rec.metric.style
    if rec.metric.max_width is not None:
        style += '; max-width: {w}px; width: {w}px;'.format(w=rec.metric.max_width)
    if rec.metric.min_width is not None:
        style += '; min-width: {w}px;'.format(w=rec.metric.min_width)

    html += '\n<td metric="' + rec.metric.name + '" style="' + style + '"'
    html += ' quality="' + str(rec.metric.quality) + '" class="td ' + td_classes + ' ' + rec.metric.class_ + ' '
    if rec.num:
        html += ' number" number="' + str(rec.value) + '" data-sortAs="' + str(rec.value)
    html += '"'
    # if rec.metric.width is not None:
    #     html += ' width=' + str(rec.metric.width)
    html += '>'

    if rec.right_shift:
        padding_style = 'margin-left: ' + str(rec.right_shift) + 'px; margin-right: -' + str(rec.right_shift) + 'px;'
    else:
        padding_style = ''

    if rec.html_fpath:
        if isinstance(rec.html_fpath, basestring):
            html += '<a href="' + rec.html_fpath + '">' + rec.cell_contents + '</a></td>'
        else:  # varQC -- several variant callers for one sample are possible
            if len(rec.html_fpath) == 0:
                rec.value = None
                _calc_record_cell_contents(rec)
                html += rec.cell_contents + '</td>'
            else:
                caller_links = ', '.join(
                    '<a href="' + html_fpath + '">' + caller + '</a>'
                    for caller, html_fpath in rec.html_fpath.items()
                    if html_fpath)
                html += rec.cell_contents + '(' + caller_links + ')</td>'

    else:
        html += '<span style="' + str(padding_style) + '" ' + \
                __get_meta_tag_contents(rec) + '>' + rec.cell_contents + '</span></td>'
    return html


def _build_total_report(report, section, column_order):
    html = ''

    if section.title:
        html += '\n<h3 class="table_name">' + section.title + '</h3>'
    calc_cell_contents(report, report.sample_reports if 'sample_reports' in report.__dict__ else [report], section)

    table = '\n<table cellspacing="0" class="report_table tableSorter fix-align-char" id="report_table_' + section.name + '">'
    table += '\n<thead>\n<tr class="top_row_tr">'
    table += '\n<th class="top_left_td left_column_td" data-sortBy="numeric"><span>Sample</span></th>'

    for col_num in range(len(section.metrics)):
        pos = column_order[col_num]
        metric = section.metrics[pos]
        table += make_cell_th(metric, pos)

    table += '\n</tr>\n</thead>\n<tbody>'

    i = 0
    if 'sample_reports' not in report.__dict__:
        sample_reports = [report]
    else:
        sample_reports = report.sample_reports
    sample_reports_length = len(sample_reports)

    for sample_report in sample_reports:
        line_caption = sample_report.display_name  # sample name
        max_sample_name_len = 50
        if len(line_caption) > max_sample_name_len:
            line_caption = '<span title="' + line_caption + '">' + line_caption[:max_sample_name_len] + '...</span>'

        second_row_tr = 'second_row_tr' if i == 0 else ''
        table += '\n<tr class="' + second_row_tr + '">\n<td class="left_column_td td" data-sortAs=' + str(sample_reports_length - i) + '>'
        if sample_reports_length == 1:
            table += '<span class="sample_name">' + line_caption + '</span>'
        else:
            if sample_report.html_fpath:
                if isinstance(sample_report.html_fpath, basestring):
                    table += '<a class="sample_name" href="' + sample_report.html_fpath + '">' + line_caption + '</a>'
                else:  # several links for one sample are possible multi-reports (e.g. TargQC)
                    if len(sample_report.html_fpath) == 0:
                        table += '<span class="sample_name">' + line_caption + '</span>'
                    else:
                        links = ', '.join('<a href="' + html_fpath + '">' + report_name + '</a>'
                                          for report_name, html_fpath in sample_report.html_fpath.items()
                                          if html_fpath)
                        table += '<span class="sample_name">' + line_caption + '(' + links + ')</span>'
            else:
                table += '<span class="sample_name">' + line_caption + '</span>'
        table += '\n</td>'

        for col_num in range(len(section.metrics)):
            pos = column_order[col_num]
            metric = section.metrics[pos]
            if not metric.values: continue
            if metric.is_hidden: continue
            rec = sample_report.find_record(sample_report.records, metric.name)
            if rec:
                table += make_cell_td(rec)

        table += '\n</tr>'
        i += 1
    table += '\n</tbody>\n</table>\n'

    html += table
    return html


def __get_meta_tag_contents(rec):
    # meta = rec.meta
    #
    # if meta? and (a for a of meta).length != 0
    #     if typeof meta is 'string'
    #         return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\" title=\"#{meta}\""
    #
    #     else  # qc
    #         (k for own k of meta).length isnt 0
    #         meta_table = '<table class=\'qc_meta_table\'>\n<tr><td></td>'
    #         for novelty, values of meta when novelty isnt 'all'
    #             meta_table += "<td>#{novelty}</td>"
    #         meta_table += '</tr>\n'
    #
    #         for novelty, val_by_db of meta
    #             dbs = (db for db, val of val_by_db when db isnt 'average')
    #             dbs.push 'average'
    #             break
    #
    #         short_table = true
    #         for novelty, val_by_db of meta
    #             if not check_all_values_equal (val for db, val of val_by_db when db isnt 'average')
    #                 short_table = false
    #
    #         if short_table  # Values are the same for each database
    #             meta_table += '<tr><td></td>'
    #             for novelty, val_by_db of meta when novelty isnt 'all'
    #                 pretty_str = toPrettyString(val_by_db[dbs[0]], rec.metric.unit)
    #                 meta_table += "<td>#{pretty_str}</td>"
    #             meta_table += '</tr>\n'
    #         else
    #             for db in dbs
    #                 meta_table += "<tr><td>#{db}</td>"
    #                 for novelty, val_by_db of meta when novelty isnt 'all'
    #                     meta_table += "<td>#{toPrettyString(val_by_db[db], rec.metric.unit)}</td>"
    #                 meta_table += '</tr>\n'
    #
    #         meta_table += '</table>\n'
    #
    #         return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\" title=\"#{meta_table}\""
    # else
    return ''
    # return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\""


def _calc_record_cell_contents(rec):
    rec.cell_contents = rec.format_html()
    if rec.value is not None and (isinstance(rec.value, int) or isinstance(rec.value, float)):
        rec.num = rec.value

    #TODO: intPartTextWidth
    return rec


# Color heat map
BLUE_HUE = 240
BLUE_OUTER_BRT = 55
BLUE_INNER_BRT = 65

GREEN_HUE = 120
GREEN_OUTER_BRT = 50
GREEN_INNER_BRT = 60

RED_HUE = 0
RED_OUTER_BRT = 50
RED_INNER_BRT = 60

MIN_NORMAL_BRT = 80
MEDIAN_BRT = 100  # just white.


def hsl2rgb(h, s, l):
    r, g, b = None, None, None

    if s == 0:
        r = g = b = l  # achromatic
    else:
        q = l * (1 + s) if l < 0.5 else l + s - l * s
        p = 2 * l - q
        r = hue2rgb(p, q, h + 1./3)
        g = hue2rgb(p, q, h)
        b = hue2rgb(p, q, h - 1./3)

    return map(int, [round(r * 255), round(g * 255), round(b * 255)])

def hue2rgb(p, q, t):
    if t < 0: t += 1
    if t > 1: t -= 1
    if t < 1./6: return p + (q - p) * 6 * t
    if t < 1./2: return q
    if t < 2./3: return p + (q - p) * (2./3 - t) * 6
    return p


def get_color(hue, lightness):
    lightness = lightness or 92
    # lightness = Math.round (Math.pow hue - 75, 2) / 350 + 35
    rgb = hsl2rgb(float(hue) / 360, 0.8, float(lightness) / 100)
    hex_rgb = [hex(c)[2:] for c in rgb]
    return '#' + ''.join(hex_rgb)


def calc_cell_contents(report, rows, section):
    max_frac_widths_by_metric = dict()

    # First round: calculatings max/min integral/fractional widths (for decimal alingment) and max/min values (for heatmaps)
    for row in rows:
        for rec in row.records:
            if rec.metric.name in section.metrics_by_name:
                _calc_record_cell_contents(rec)

        for rec in row.records:
            if rec.metric.name in section.metrics_by_name:
                if rec.metric.name not in max_frac_widths_by_metric or \
                                rec.frac_width > max_frac_widths_by_metric[rec.metric.name]:
                    max_frac_widths_by_metric[rec.metric.name] = rec.frac_width
                if rec.num:
                    rec.metric.numbers.append(rec.num)
                if rec.value is not None:
                    rec.metric.values.append(rec.value)
    # else:
    #     for rec in report.records:
    #         if rec.metric.name in section.metrics_by_name:
    #             if rec.num:
    #                 if not rec.metric.values:
    #                     rec.metric.values = []
    #                 rec.metric.values.append(rec.num)

    for metric in section.metrics:
        if metric.numbers and metric.with_heatmap:
            # For metrics where we know the "normal value" - we want to color everything above normal white,
            #   and everything below - red, starting from normal, finishing with bottom
            if metric.ok_threshold is not None:
                if isinstance(metric.ok_threshold, int) or isinstance(metric.ok_threshold, float):
                    metric.med = metric.ok_threshold
                    if metric.bottom is not None:
                        metric.low_outer_fence = metric.bottom

            def _cmp(a, b):  # None is always less than anything
                if a and b:
                    return cmp(a, b)
                elif a:
                    return 1
                else:
                    return -1

            numbers = sorted(
                [v for v in metric.numbers],
                cmp=_cmp)
            l = len(numbers)

            metric.min = numbers[0]
            metric.max = numbers[l - 1]
            metric.all_values_equal = metric.min == metric.max
            if metric.med is None:
                metric.med = numbers[(l - 1) / 2] if l % 2 != 0 else mean([numbers[l / 2], numbers[(l / 2) - 1]])
            q1 = numbers[int(floor((l - 1) / 4))]
            q3 = numbers[int(floor((l - 1) * 3 / 4))]

            d = q3 - q1
            metric.low_outer_fence = metric.low_outer_fence if metric.low_outer_fence is not None else q1 - 3   * d
            metric.low_inner_fence = metric.low_inner_fence if metric.low_inner_fence is not None else q1 - 1.5 * d
            metric.top_inner_fence = metric.top_inner_fence if metric.top_inner_fence is not None else q3 + 1.5 * d
            metric.top_outer_fence = metric.top_outer_fence if metric.top_outer_fence is not None else q3 + 3   * d

    # Second round: setting shift and color properties based on max/min widths and vals
    for row in rows:
        for rec in row.records:
            if rec.metric and rec.metric.name in section.metrics_by_name:
                # Padding based on frac width
                if rec.frac_width:
                    rec.right_shift = max_frac_widths_by_metric[rec.metric.name] - rec.frac_width

                metric = rec.metric

                # Color heatmap
                if rec.num and metric.with_heatmap:
                    [top_hue, inner_top_brt, outer_top_brt] = [BLUE_HUE, BLUE_INNER_BRT, BLUE_OUTER_BRT]
                    [low_hue, inner_low_brt, outer_low_brt] = [RED_HUE, RED_INNER_BRT, RED_OUTER_BRT]

                    if metric.quality == 'Less is better':  # then swap colors
                        [top_hue, low_hue] = [low_hue, top_hue]
                        [inner_top_brt, inner_low_brt] = [inner_low_brt, inner_top_brt]
                        [outer_top_brt, outer_low_brt] = [outer_low_brt, outer_top_brt]

                    if metric.ok_threshold is not None:
                        if isinstance(metric.ok_threshold, int) or isinstance(metric.ok_threshold, float):
                            if rec.num >= metric.ok_threshold:
                                continue  # white on blak

                            # rec_to_align_with = sample_report.find_record(sample_report.records, metric.threshold)
                            # if rec_to_align_with:
                            #     rec.text_color = lambda: rec_to_align_with.text_color()
                            #     continue

                    if not metric.all_values_equal:
                        rec.text_color = 'black'

                        # Low outliers
                        if rec.num < rec.metric.low_outer_fence:
                            rec.color = get_color(low_hue, outer_low_brt)
                            rec.text_color = 'white'

                        elif rec.num < rec.metric.low_inner_fence:
                            rec.color = get_color(low_hue, inner_low_brt)

                        # Normal values
                        elif rec.num < metric.med:
                            k = (MEDIAN_BRT - MIN_NORMAL_BRT) / (metric.med - rec.metric.low_inner_fence)
                            brt = round(MEDIAN_BRT - (metric.med - rec.num) * k)
                            rec.color = get_color(low_hue, brt)

                        # High outliers
                        elif rec.num > rec.metric.top_inner_fence:
                            rec.color = get_color(top_hue, inner_top_brt)

                        elif rec.num > rec.metric.top_outer_fence:
                            rec.color = get_color(top_hue, outer_top_brt)
                            rec.text_color = 'white'

                        elif rec.num > metric.med:
                            k = (MEDIAN_BRT - MIN_NORMAL_BRT) / (rec.metric.top_inner_fence - metric.med)
                            brt = round(MEDIAN_BRT - (rec.num - metric.med) * k)
                            rec.color = get_color(top_hue, brt)

        for rec in row.records:
            if rec.metric and rec.metric.name in section.metrics_by_name:
                metric = rec.metric

                if metric.ok_threshold is not None:
                    if isinstance(metric.ok_threshold, basestring):
                        rec_to_align_with = Report.find_record(row.records, metric.ok_threshold)
                        if rec_to_align_with:
                            rec.text_color = rec_to_align_with.text_color
                            rec.color = rec_to_align_with.color
    return report


def write_html_report(cnf, report, type_, html_fpath, caption='',
        extra_js_fpaths=list(), extra_css_fpaths=list(), image_by_key=dict()):

    with open(template_fpath) as f: html = f.read()

    html = _insert_into_html(html, caption, 'caption')
    html = _insert_into_html(html, datetime.datetime.now().strftime('%d %B %Y, %A, %H:%M:%S'), 'report_date')

    report_html = _build_report(report, type_)
    html = _insert_into_html(html, report_html, 'report')

    plots_html = ''.join('<img src="' + plot + '"/>' for plot in report.plots)
    html = _insert_into_html(html, plots_html, 'plots')

    return __write_html(cnf, html, html_fpath,
        extra_js_fpaths, extra_css_fpaths, image_by_key)


# def _copy_aux_files(results_dirpath, extra_files=list()):
#     aux_dirpath = join(results_dirpath, aux_dirname)
#
#     if isdir(aux_dirpath):
#         shutil.rmtree(aux_dirpath)
#
#     if not isdir(aux_dirpath):
#         try:
#             os.mkdir(aux_dirpath)
#         except OSError:
#             pass
#
#     def copy_aux_file(fname):
#         src_fpath = join(static_dirpath, fname)
#
#         if file_exists(src_fpath):
#             dst_fpath = join(aux_dirpath, fname)
#
#             if not file_exists(dirname(dst_fpath)):
#                 try:
#                     os.makedirs(dirname(dst_fpath))
#                 except OSError:
#                     pass
#             try:
#                 shutil.copyfile(src_fpath, dst_fpath)
#             except OSError:
#                 pass
#
#     for aux_f_relpath in js_files + css_files + image_files + extra_files:
#         if aux_f_relpath.endswith('.js'):
#             for ext in ['.js', '.coffee', '.map']:
#                 copy_aux_file(splitext(aux_f_relpath)[0] + ext)
#
#         elif aux_f_relpath.endswith('.css'):
#             for ext in ['.css', '.sass']:
#                 copy_aux_file(splitext(aux_f_relpath)[0] + ext)
#
#         else:
#             copy_aux_file(aux_f_relpath)


import base64

def _embed_images(html, report_dirpath, image_by_key, debug=False):
    ptrn = '<img src="{key}"'

    for key, fpath in image_by_key.items():
        info('Embedding image ' + fpath + '...', ending=' ')
        if not verify_file(fpath, silent=True):
            fpath = join(static_dirpath, join(*fpath.split('/')))
            if not verify_file(fpath):
                continue

        old_piece = ptrn.format(key=key)
        if debug:  # Not embedding, just adding links
            new_piece = old_piece.replace(key, relpath(fpath, report_dirpath))
        else:
            with open(fpath, 'rb') as f:
                encoded_string = base64.b64encode(f.read())
            binary_image = 'data:image/png;base64,' + encoded_string
            new_piece = old_piece.replace(key, binary_image)

        html = html.replace(old_piece, new_piece)
        info('Done.', print_date=False)

    return html


def _embed_css_and_scripts(html, report_dirpath,
       extra_js_fpaths=list(), extra_css_fpaths=list(), debug=False):
    js_line_tmpl = '<script type="text/javascript" src="{file_rel_path}"></script>'
    js_l_tag = '<script type="text/javascript" name="{name}">'
    js_r_tag = '    </script>'

    css_line_tmpl = '<link rel="stylesheet" type="text/css" href="{file_rel_path}" />'
    css_l_tag = '<style type="text/css" rel="stylesheet" name="{name}">'
    css_r_tag = '    </style>'

    for line_tmpl, files, l_tag, r_tag in [
            (js_line_tmpl,  extra_js_fpaths  + js_files,  js_l_tag,  js_r_tag),
            (css_line_tmpl, extra_css_fpaths + css_files, css_l_tag, css_r_tag),
        ]:
        for rel_fpath in files:
            info('Embedding ' + rel_fpath + '...', ending=' ')
            if verify_file(rel_fpath, silent=True):
                fpath = rel_fpath
                rel_fpath = basename(fpath)
            else:
                fpath = join(static_dirpath, join(*rel_fpath.split('/')))
                if not verify_file(fpath):
                    continue

            line = line_tmpl.format(file_rel_path=rel_fpath)
            l_tag_formatted = l_tag.format(name=rel_fpath)

            if debug:  # Not embedding, just adding links
                fpath = relpath(fpath, report_dirpath)
                line_formatted = line.replace(rel_fpath, fpath)
                html = html.replace(line, line_formatted)

            else:
                with open(fpath) as f:
                    contents = f.read()
                    contents = '\n'.join(' ' * 8 + l for l in contents.split('\n'))

                    # line_ascii = line.encode('ascii', 'ignore')
                    # html = html.decode('utf-8', 'ignore').encode('ascii')
                    # try:
                    #     html_ascii = html.decode('utf-8', 'ignore').encode('ascii')  #.encode('ascii', 'ignore')
                    # except:
                    #     pass
                    try:
                        html = ''.join(i for i in html if ord(i) < 128)
                        contents = ''.join(i for i in contents if ord(i) < 128)
                        html = html.replace(line, l_tag_formatted + '\n' + contents + '\n' + r_tag)
                    except:
                        err()
                        # print contents
                        err()
                        err(traceback.format_exc())
                        err()
                        err('Encoding problem embeding this file into html: ' + rel_fpath, print_date=False)
                        err()
                        continue

                # try:
                #     html = html_ascii.encode('ascii').replace(line, l_tag_fmt + '\n' + contents + '\n' + r_tag)
                # except:
                #     pass
                    # err('line - cannot replace')  # + str(l_tag_fmt + '\n' + contents + '\n' + r_tag))

            info('Done.', print_date=False)
    return html


def _insert_into_html(html, text, keyword):
    # substituting template text with json
    # html_text = re.sub('{{ ' + keyword + ' }}', text, html_text)
    html = ''.join(i for i in html if ord(i) < 128)
    html = html.replace('{{ ' + keyword + ' }}', text)

    return html


def write_static_html_report(cnf, data_dict, html_fpath, tmpl_fpath=None,
        extra_js_fpaths=list(), extra_css_fpaths=list(), image_by_key=dict()):

    tmpl_fpath = tmpl_fpath or static_template_fpath
    with open(tmpl_fpath) as f: html = f.read()

    html = jsontemplate.expand(html, data_dict)

    return __write_html(cnf, html, html_fpath,
        extra_js_fpaths, extra_css_fpaths, image_by_key)


def __write_html(cnf, html, html_fpath, extra_js_fpaths, extra_css_fpaths, image_by_key):
    html = _embed_css_and_scripts(html, dirname(html_fpath), extra_js_fpaths, extra_css_fpaths, cnf.debug)
    html = _embed_images(html, dirname(html_fpath), image_by_key, cnf.debug)

    with file_transaction(cnf.work_dir, html_fpath) as tx:
        with open(tx, 'w') as f:
            f.write(html)

    return html_fpath
