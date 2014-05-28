from collections import OrderedDict
from os.path import basename
from os.path import splitext

database = 'cosmic'
novelty = 'all'
metrics_header = 'Metric'
novelty_header = 'Novelty'
sample_header = 'Sample name:'


def summarize_qc(report_fpaths, output_summary_fpath, report_suffix=None):
    full_report = [['Sample']]
    for report_fpath in report_fpaths:
        _, report_dict = _parse_report_qc(report_fpath)
        sample_name = _get_sample_name(report_fpath, report_suffix)
        _add_to_full_report(full_report, sample_name, report_dict)
    _print_full_report(full_report, output_summary_fpath)


def summarize_cov(report_fpaths, output_summary_fpath, report_suffix=None):
    full_report = [['Sample']]
    for report_fpath in report_fpaths:
        report_dict = _parse_report_cov(report_fpath)
        sample_name = _get_sample_name(report_fpath, report_suffix)
        _add_to_full_report(full_report, sample_name, report_dict)
    _print_full_report(full_report, output_summary_fpath)


def _get_sample_name(report_fpath, report_suffix=None):
    if report_fpath.endswith(report_suffix):
        return basename(report_fpath)[:-len(report_suffix)]
    return splitext(basename(report_fpath))[0]


def _parse_report_qc(report_fpath):
    sample_name = ''
    report_dict = OrderedDict()
    with open(report_fpath, 'r') as f:
        # parsing Sample name and Database columns
        database_col_id = None
        novelty_col_id = None
        for line in f:
            if line.startswith(sample_header):
                sample_name = line[len(sample_header):].strip()
            elif line.startswith(metrics_header):
                if database in line:
                    database_col_id = line.split().index(database)
                if novelty_header in line:
                    novelty_col_id = line.split().index(novelty_header)
                break

        if database_col_id:
            # parsing rest of the report
            for line in f:
                if novelty_col_id and line.split()[novelty_col_id] != novelty:
                    continue
                cur_metric_name = line.split()[0]
                cur_value = line.split()[database_col_id]
                report_dict[cur_metric_name] = cur_value
    return sample_name, report_dict


def _parse_report_cov(report_fpath):
    report_dict = OrderedDict()
    with open(report_fpath, 'r') as f:
        for line in f:
            cur_value = line.split()[-1]
            cur_metric_name = line.strip()[:-len(cur_value)].strip()
            report_dict[cur_metric_name] = cur_value
    return report_dict


def _add_to_full_report(full_report, sample_name, report_dict):
    full_report[0].append(sample_name)
    if len(full_report) == 1:
        empty_report = True
    else:
        empty_report = False

    for id, (key, value) in enumerate(report_dict.items()):
        if empty_report:
            full_report.append([key])
        full_report[id + 1].append(value)


def _print_full_report(report, report_filename):
    col_widths = [0] * len(report[0])
    for row in report:
        for id, value in enumerate(row):
            col_widths[id] = max(len(value), col_widths[id])

    out = open(report_filename, 'w')
    for row in report:
        out.write('\t'.join(row) + '\n')
    out.close()


