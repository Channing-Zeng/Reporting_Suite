from collections import OrderedDict

database = 'cosmic'
novelty = 'all'
metrics_header = 'Metric'
novelty_header = 'Novelty'
sample_header = 'Sample name:'


def summarize_qc(input_reports, output_summary_fpath):
    full_report = [['Sample']]
    for report in input_reports:
        sample_name, report_dict = _parse_report(report)
        _add_to_full_report(full_report, sample_name, report_dict)
    _print_full_report(full_report, output_summary_fpath)


def _parse_report(report_filename):
    report_handler = open(report_filename, 'r')
    sample_name = ''
    report_dict = OrderedDict()
    # parsing Sample name and Database columns
    database_col_id = None
    novelty_col_id = None
    for line in report_handler:
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
        for line in report_handler:
            if novelty_col_id and line.split()[novelty_col_id] != novelty:
                continue
            cur_metric_name = line.split()[0]
            cur_value = line.split()[database_col_id]
            ### HARD_CODE!
            if cur_metric_name == 'variantRatePerBp':
                cur_metric_name = 'basesPerVariant'
                additional_metric_name = 'variantRate'
                additional_value = '%.8f' % (1.0 / float(cur_value))
                report_dict[additional_metric_name] = additional_value
            ###
            report_dict[cur_metric_name] = cur_value

    report_handler.close()
    return sample_name, report_dict


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
        out.write('  '.join('%-*s' % (col_width, value)
                            for col_width, value in zip(col_widths, row))
                  + "\r\n")
    out.close()
