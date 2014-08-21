from collections import OrderedDict
from itertools import repeat, izip, chain
import json
from os.path import join
from source.file_utils import verify_file
from source.quast_reporting.html_saver import write_html_reports

from source.logger import critical, info
from source.file_utils import file_exists


class Record:
    def __init__(self,
                 metric=None,
                 value=None,
                 meta=None):
        self.metric = metric
        self.value = value
        self.meta = meta or dict()

    @staticmethod
    def dump_records(records, f):
        objects = [rec.__dict__ for rec in records]

        json.dump(objects, f,
                  default=lambda o: o.__dict__,
                  sort_keys=True,
                  indent=4)

    @staticmethod
    def load_records(fpath):
        with open(fpath) as f:
            records = []
            objects = json.load(f, object_pairs_hook=OrderedDict)
            for obj in objects:
                rec = Record()
                rec.__dict__ = dict(obj.items())
                m = Metric()
                m.__dict__.update(dict(rec.metric.items()))
                rec.metric = m
                records.append(rec)
            return records

    def format(self):
        return self.metric.format(self.value)


class Metric:
    def __init__(self,
                 name=None,
                 short_name=None,
                 description=None,
                 quality='More is better',  # More is better, Less is better, Equal
                 unit='',
                 common=False):
        self.name = name
        self.short_name = short_name
        self.description = description
        self.quality = quality
        self.common = common
        self.unit = unit

    @staticmethod
    def to_dict(metrics):
        return OrderedDict((m.name, m) for m in metrics)

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


class FullReport:
    def __init__(self, name='', sample_reports=list()):
        self.name = name
        self.sample_reports = sample_reports

    def copy(self):
        return FullReport(
            name=self.name,
            sample_reports=self.sample_reports[:])


class SampleReport:
    def __init__(self, sample=None, fpath=None, records=list(), name='', plots=list()):
        self.sample = sample
        self.fpath = fpath
        self.records = records
        self.name = name
        self.plots = plots  # TODO: make real JS plots, not just included PNG


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


def summarize(cnf, report_fpath_by_sample, parse_report_fn, report_name):
    """ Returns list of SampleReport objects:
        [SampleReport(sample=Sample(name=), fpath=, records=[Record,...]),...]
    """
    return FullReport(
        name=report_name,
        sample_reports=[
            SampleReport(sample_name, fpaths.html_fpath, parse_report_fn(cnf, fpaths.content_fpath))
            for sample_name, fpaths in report_fpath_by_sample.items()])


def write_summary_reports(output_dirpath, work_dirpath, full_reports, base_fname, caption):
    if not isinstance(full_reports, list):
        full_reports = [full_reports]

    return [fn(output_dirpath, work_dirpath, full_reports, base_fname, caption)
        for fn in [write_txt_reports,
                   write_tsv_reports,
                   write_html_reports]]


def _flatten_report(full_reports):
    # new_full_report = full_reports[0].copy()
    #
    # for full_report in full_reports[1:]:
    #     for new_sample_rep, sample_report in \
    #         zip(new_full_report.sample_reports,
    #             full_report.sample_reports):
    #
    #         new_sample_rep.extend(records)
    #
    #     new_full_report.sample_reports.extend(full_report.sample_reports)

    if isinstance(full_reports, FullReport):
        full_reports = [full_reports]
    else:
        if isinstance(full_reports, SampleReport):
            full_reports = [full_reports]

        if isinstance(full_reports, list) and isinstance(full_reports[0], SampleReport):
            full_reports = [FullReport(sample_reports=full_reports)]

    rows = [['Sample'] + [rep.sample.name for rep in full_reports[0].sample_reports]]

    for full_report in full_reports:
        for metric in (r.metric for r in full_report.sample_reports[0].records):
            row = [metric.name]
            for sample_report in full_report.sample_reports:
                row.append(next(
                    r.metric.format(r.value)
                    for r in sample_report.records
                    if r.metric.name == metric.name))
            rows.append(row)
    return rows


def write_txt_reports(output_dirpath, work_dirpath, full_reports, base_fname, caption=None):
    rows = _flatten_report(full_reports)
    return write_txt_rows(rows, output_dirpath, base_fname)


def write_txt_rows(rows, output_dirpath, base_fname):
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


def write_tsv_reports(output_dirpath, work_dirpath, full_reports, base_fname, caption=None):
    rows = _flatten_report(full_reports)
    return write_tsv_rows(rows, output_dirpath, base_fname)


def write_tsv_rows(rows, output_dirpath, base_fname):
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


def save_json(records, fpath):
    with open(fpath, 'w') as f:
        Record.dump_records(records, f)


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
