from os.path import join

from source.logger import step_greetings
from source.utils import get_java_tool_cmdline, call_subprocess


def gatk_qc(cnf, qc_dir, vcf_fpath):
    step_greetings('Quality control reports')

    log = cnf['log']
    work_dir = cnf['work_dir']

    qc_cnf = cnf['quality_control']
    databases = qc_cnf.get('database_vcfs')
    novelty = qc_cnf.get('novelty')
    metrics = qc_cnf.get('metrics')

    executable = get_java_tool_cmdline(cnf, 'gatk')
    gatk_opts_line = ' '.join(cnf.get('gatk', {'options': []}).get('options', []))
    if 'threads' in cnf and ' -nt ' not in gatk_opts_line:
        gatk_opts_line += ' -nt ' + cnf.get('threads', '1')

    ref_fpath = cnf['genome']['seq']
    report_fpath = join(work_dir, cnf['name'] + '_gatk.report')

    cmdline = ('{executable} {gatk_opts_line} -R {ref_fpath} -T VariantEval'
               ' --eval:tmp {vcf_fpath} -o {report_fpath}').format(**locals())

    if 'dbsnp' in databases:
        cmdline += ' -D ' + databases['dbsnp']
    for db_name, db_path in databases.items():
        if not db_name == 'dbsnp':
            cmdline += ' -comp:' + db_name + ' ' + db_path

    call_subprocess(cnf, cmdline, None, report_fpath, stdout_to_outputfile=False,
         to_remove=[vcf_fpath + '.idx'])

    report = _parse_gatk_report(report_fpath, databases.keys(), novelty, metrics)

    final_report_fpath = join(qc_dir, cnf['name'] + '_qc.report')

    _make_final_report(report, final_report_fpath, cnf['name'],
                       databases.keys(), novelty, metrics)
    return final_report_fpath


def _parse_gatk_report(report_filename, databases, novelty, metrics):
    database_col_name = 'CompRod'
    database_col_id = None
    novelty_col_name = 'Novelty'
    novelty_col_id = None

    report = dict()
    comments_section = False
    cur_header = []
    cur_metrics_ids = []
    for line in open(report_filename):
        if not line.strip():
            continue
        if line.startswith('#'): # comment line
            comments_section = True
            continue
        elif comments_section:
            comments_section = False
            cur_header = line.split()
            cur_metrics_ids = []
            database_col_id = cur_header.index(database_col_name)
            novelty_col_id = cur_header.index(novelty_col_name)
            for metric in metrics:
                if metric in cur_header:
                    cur_metrics_ids.append(cur_header.index(metric))
                    if metric not in report:
                        report[metric] = dict()
        elif cur_metrics_ids:  # process lines only if there are metrics in current section
            values = line.split()
            cur_database = values[database_col_id]
            cur_novelty = values[novelty_col_id]
            if (cur_database not in databases) or (cur_novelty not in novelty):
                continue
            for metric_id in cur_metrics_ids:
                if cur_database not in report[cur_header[metric_id]]:
                    report[cur_header[metric_id]][cur_database] = dict()
                report[cur_header[metric_id]][cur_database][cur_novelty] = values[metric_id]
    return report


def _make_final_report(report_dict, report_filename, sample_name,
                       databases, novelty, metrics):
    header = ['Metric', 'Novelty'] + databases + ['Average']
    full_report = [header]
    for cur_metric in metrics:
        for cur_novelty in novelty:
            cur_row = [cur_metric, cur_novelty]
            sum = 0.0
            for cur_database in databases:
                if cur_metric == 'variantRatePerBp': # confusing name and value format
                    cur_row[0] = 'basesPerVariant'
                    cur_row.append("%.2f" % float(report_dict[cur_metric][cur_database][cur_novelty]))
                else:
                    cur_row.append(report_dict[cur_metric][cur_database][cur_novelty])
                sum += float(cur_row[-1])
            average = sum / len(databases)
            cur_row.append("%.2f" % average)
            full_report.append(cur_row)

    col_widths = [0] * len(header)
    for row in full_report:
        for id, value in enumerate(row):
            col_widths[id] = max(len(value), col_widths[id])

    out = open(report_filename, 'w')
    out.write('Sample name: ' + sample_name + '\n\n')
    for row in full_report:
        out.write('  '.join('%-*s' % (col_width, value) for col_width, value
                            in zip(col_widths, row)) + "\r\n")
    out.close()

