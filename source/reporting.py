from collections import OrderedDict
from itertools import repeat, izip
from os.path import join
from source.file_utils import verify_file
from source.quast_reporting.html_saver import write_html_report

from source.logger import critical, info
from source.utils_from_bcbio import file_exists


class Metric(object):
    def __init__(self):
        self.name = None
        self.quality = None
        self.value = None
        self.link = None
        self.meta = dict()


class SampleReport():
    def __init__(self, name=None, fpath=None, metrics=list()):
        self.name = name
        self.fpath = fpath
        self.metrics = metrics


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

    for metric in report[0].metrics:
        row = [metric.name]
        for sample in report:
            row.append(next(m.value for m in sample.metrics if m.name == metric.name))
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
                try:
                    val_int = int(val)
                except ValueError:
                    try:
                        val_double = float(val)
                    except ValueError:
                        pass
                    else:
                        val = '{0:.2f}'.format(val_double)
                else:
                    pass
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
            out.write('\t'.join(row) + '\n')

    return output_fpath


def parse_tsv(tsv_fpath):
    report = []

    with open(tsv_fpath, 'r') as f:
        for line in f:
            report.append(line.split('\t'))

    if not report:
        critical('Data not found in ' + tsv_fpath)

    return report


