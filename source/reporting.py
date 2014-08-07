from collections import OrderedDict
from itertools import repeat, izip
import json
from os.path import join
from source.file_utils import verify_file
from source.quast_reporting.html_saver import write_html_report

from source.logger import critical, info
from source.file_utils import file_exists


class Record(object):
    def __init__(self,
                 metric=None,
                 value=None,
                 meta=None):
        self.metric = metric
        self.value = value
        self.meta = meta or dict()

    @staticmethod
    def dump_records(records, f):
        objects = {name: rec.__dict__ for name, rec in records.items()}

        json.dump(objects, f,
                  default=lambda o: o.__dict__,
                  sort_keys=True,
                  indent=4)

    @staticmethod
    def load_records(f):
        records = dict()
        objects = json.load(f)
        for rec_name, obj in objects.items():
            rec = Record()
            rec.__dict__ = obj
            m = Metric()
            m.__dict__ = rec.metric
            rec.metric = m
            records[rec_name] = rec
        return records


class Metric(object):
    def __init__(self,
                 name=None,
                 short_name=None,
                 description=None,
                 quality='More is better',  # More is better, Less is better
                 unit=''):
        self.name = name
        self.short_name = short_name,
        self.description = description,
        self.quality = quality
        self.unit = unit

    @staticmethod
    def to_dict(metrics):
        return {m.name: m for m in metrics}

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
                if value < 1/i:
                    presision = i
            return '{value:.{presision}f}{unit}'.format(**locals())


class SampleReport():
    def __init__(self, name=None, fpath=None, records=list()):
        self.name = name
        self.fpath = fpath
        self.records = records


def read_sample_names(sample_fpath):
    sample_names = []

    with open(sample_fpath) as f:
        for line in f:
            sample_name = line.strip()
            sample_names.append(sample_name)

    return sample_names


def get_sample_report_fpaths_for_bcbio_final_dir(
        bcbio_final_dir, sample_names, varqc_dir, ending):

    single_report_fpaths = []

    fixed_sample_names = []
    for sample_name in sample_names:
        single_report_fpath = join(
            bcbio_final_dir, sample_name, varqc_dir,
            sample_name + ending)

        info(single_report_fpath)

        if not file_exists(single_report_fpath):
            info('No ' + single_report_fpath + ', skipping.')
            continue

        if not verify_file(single_report_fpath):
            critical(single_report_fpath + ' does not exist.')

        single_report_fpaths.append(single_report_fpath)
        fixed_sample_names.append(sample_name)

    return single_report_fpaths, fixed_sample_names


def summarize(sample_names, report_fpaths, parse_report_fn):
    return [SampleReport(name, fpath, parse_report_fn(fpath).values())
            for name, fpath in zip(sample_names, report_fpaths)]


def write_summary_reports(output_dirpath, work_dirpath, report, base_fname, caption):

    return [fn(output_dirpath, work_dirpath, report, base_fname, caption)
        for fn in [write_txt_report,
                   write_tsv_report,
                   write_html_report]]


def _flatten_report(report):
    rows = [['Sample'] + [s.name for s in report]]

    for record in report[0].records:
        row = [record.metric.name]
        for sample in report:
            row.append(next(
                r.metric.format(r.value)
                for r in sample.records
                if r.metric.name == record.metric.name))
        rows.append(row)

    return rows


def write_txt_report(output_dirpath, work_dirpath, report, base_fname, caption=None):
    rows = _flatten_report(report)
    return write_txt(rows, output_dirpath, base_fname)


def write_txt(rows, output_dirpath, base_fname):
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


def write_tsv_report(output_dirpath, work_dirpath, report, base_fname, caption=None):
    rows = _flatten_report(report)
    return write_tsv(rows, output_dirpath, base_fname)


def write_tsv(rows, output_dirpath, base_fname):
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