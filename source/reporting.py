from collections import OrderedDict
from itertools import repeat, izip
from os.path import join
from source.file_utils import verify_file
from source.quast_reporting.html_saver import write_html_report

from source.logger import critical, info
from source.file_utils import file_exists


class Record(object):
    def __init__(self,
                 metric=None,
                 value=None,
                 meta=dict()):
        self.metric = metric
        self.value = value
        self.meta = meta


class Metric(object):
    def __init__(self,
                 name=None,
                 short_name=None,
                 description=None,
                 presision=0,  # number of decimal digits
                 quality='More is better',  # More is better, Less is better
                 unit=''):
        self.name = name
        self.short_name = short_name,
        self.description = description,
        self.presision = presision
        self.quality = quality
        self.unit = unit

    def format(self, value):
        if value is None:
            return '-'

        name = self.name
        unit = self.unit
        presision = self.presision

        if isinstance(value, basestring):
            return '{value}{unit}'.format(**locals())

        try:
            value = int(value)
        except ValueError:
            value = float(value)

        if self.presision == 0:
            return '{value:,}{unit}'.format(**locals())

        else:
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
            row.append(next(r.metric.format(r.value) for r in sample.records if r.metric.name == record.metric.name))
        rows.append(row)

    return rows


def write_txt_report(output_dirpath, work_dirpath, report, base_fname, caption=None):
    rows = _flatten_report(report)
    return write_txt(rows, output_dirpath, base_fname)


def write_txt(rows, output_dirpath, base_fname):
    output_fpath = join(output_dirpath, base_fname + '.txt')

    col_widths = repeat(0)
    col_widths = [max(len(v), w) for v, w in izip(rows[0], col_widths)]

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


