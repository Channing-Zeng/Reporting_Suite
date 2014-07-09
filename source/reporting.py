from itertools import repeat, izip
from os.path import join, basename
from source.file_utils import verify_file
from source.quast_reporting.html_saver import write_html_report

from source.logger import critical, info
from source.utils import OrderedDefaultDict
from source.utils_from_bcbio import file_exists


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

    for sample_name in sample_names:
        single_report_fpath = join(
            bcbio_final_dir, sample_name, varqc_dir,
            sample_name + ending)

        info(single_report_fpath)

        if (not file_exists(single_report_fpath) and
                'mutect' in single_report_fpath and
                'vardict' in single_report_fpath):
            info('No ' + single_report_fpath + ', skipping.')
            continue

        if not verify_file(single_report_fpath):
            critical(single_report_fpath + ' does not exist.')
        single_report_fpaths.append(single_report_fpath)

    return single_report_fpaths


def summarize(sample_names, report_fpaths, get_rows_fn):
    metric_values = OrderedDefaultDict(list)
    metric_info = OrderedDefaultDict(dict)

    for sample_name, report_fpath in zip(sample_names, report_fpaths):
        for metric in get_rows_fn(report_fpath):
            metric_values[metric['metricName']].append(metric['value'])
            metric_info[metric['metricName']]['isMain'] = metric['isMain']
            metric_info[metric['metricName']]['quality'] = metric['quality']

    group = []
    report = [['Sample', group]]

    for (m_name1, metric_info), (m_name2, metric_values) \
            in zip(metric_info.items(), metric_values.items()):
        assert m_name1 == m_name2
        group.append(dict(
            metricName=m_name1, values=metric_values,
            isMain=metric_info['isMain'], quality=metric_info['quality']))

    return report


def write_summary_reports(output_dirpath, work_dirpath, report,
                          sample_names, base_fname, caption):

    return [fn(output_dirpath, work_dirpath, report,
               sample_names, base_fname, caption)
        for fn in [write_txt_report,
                   write_tsv_report,
                   write_html_report]]


def _flatten_report(report, sample_names):
    rows = [['Sample'] + sample_names]
    for group_name, group_metrics in report:
        for metric in group_metrics:
            if metric['isMain']:
                rows.append([metric['metricName']] + metric['values'])

    return rows


def write_txt_report(output_dirpath, work_dirpath, report,
                     sample_names, base_fname, caption=None):
    rows = _flatten_report(report, sample_names)
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


def write_tsv_report(output_dirpath, work_dirpath, report,
                     sample_names, base_fname, caption=None):
    rows = _flatten_report(report, sample_names)
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


