from itertools import repeat, izip
from os.path import join
import sys
from source.quast_reporting.html_saver import _init_html, save_total_report

from source.utils import OrderedDefaultDict
from source.logger import critical


def summarize(sample_names, report_fpaths, get_pairs_fn):
    report = OrderedDefaultDict(list)  # metric -> [values]

    for sample_name, report_fpath in zip(sample_names, report_fpaths):
        for metric, value in get_pairs_fn(report_fpath):
            report[metric].append(value)

    return report


def write_summary_reports(report, sample_names, output_dirpath, base_fname):
    rows = [['Sample'] + sample_names]
    for metric, values in report.items():
        rows.append([metric] + values)

    return [fn(rows, sample_names, output_dirpath, base_fname)
        for fn in [write_txt_report,
                   write_tsv_report,
                   write_html_report]]


def write_html_report(rows, sample_names, output_dirpath, base_fname):
    group_metrics = []
    group = ['', group_metrics]
    report = [group]

    for row in rows[1:]:
        group_metrics.append(dict(metricName=row[0], values=row[1:],
                                  isMain=True, quality='More is better'))

    return save_total_report(output_dirpath, sample_names, base_fname, report)


def write_txt_report(rows, sample_names, output_dirpath, base_fname=None):
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


def write_tsv_report(rows, sample_names, output_dirpath, base_fname):
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


