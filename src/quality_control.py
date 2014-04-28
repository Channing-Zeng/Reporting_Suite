from genericpath import isfile
import shutil
from os import mkdir
from os.path import basename, join

from src.utils import file_exists
from src.my_utils import info, err, verify_file, step_greetings, \
    get_tool_cmdline, get_java_tool_cmdline, call


def quality_control(cnf, qc_dir, vcf_fpath):
    if 'quality_control' not in cnf:
        return None, None

    if file_exists(qc_dir):
        shutil.rmtree(qc_dir)
    mkdir(qc_dir)

    qc_report_fpath = gatk_qc(cnf, qc_dir, vcf_fpath)
    qc_plots_fpaths = bcftools_qc(cnf, qc_dir, vcf_fpath)
    return qc_report_fpath, qc_plots_fpaths


def gatk_qc(cnf, qc_dir, vcf_fpath):
    step_greetings(cnf, 'Quality control reports')

    log = cnf['log']
    work_dir = cnf['work_dir']

    qc_cnf = cnf['quality_control']
    databases = qc_cnf.get('database_vcfs')
    novelty = qc_cnf.get('novelty')
    metrics = qc_cnf.get('metrics')

    executable = get_java_tool_cmdline(cnf, 'gatk')
    ref_fpath = cnf['genome']['seq']
    report_fpath = join(work_dir, cnf['name'] + '_gatk.report')

    cmdline = ('{executable} -nt 20 -R {ref_fpath} -T VariantEval'
               ' --eval:tmp {vcf_fpath} -o {report_fpath}').format(**locals())

    if 'dbsnp' in databases:
        cmdline += ' -D ' + databases['dbsnp']
        del databases['dbsnp']
    for db_name, db_path in databases.items():
        cmdline += ' -comp:' + db_name + ' ' + db_path

    call(cnf, cmdline, None, report_fpath, stdout_to_outputfile=False,
         to_remove=[vcf_fpath + '.idx'])

    report = _parse_gatk_report(report_fpath, databases.keys(), novelty, metrics)

    final_report_fpath = join(qc_dir, cnf['name'] + '.qc.report')

    _make_final_report(report, final_report_fpath, cnf['name'],
                       databases.keys(), novelty, metrics)
    return final_report_fpath


def _check_quality_control_config(cnf):
    qc_cnf = cnf.get('quality_control')
    if not qc_cnf:
        return

    if 'databases' not in qc_cnf:
        qc_cnf['databases'] = ['dbsnp']
        info('Warning: not databases for quality control, using [dbsnp]')

    if 'novelty' not in qc_cnf:
        qc_cnf['novelty'] = ['all', 'known', 'novel']
        info('Warning: no novelty specified for quality control, '
             'using default ' + ', '.join(qc_cnf['novelty']))

    if 'metrics' not in qc_cnf:
        qc_cnf['metircs'] = [
           'nEvalVariants', 'nSNPs', 'nInsertions', 'nDeletions',
           'nVariantsAtComp', 'compRate', 'nConcordant', 'concordantRate',
           'variantRate', 'variantRatePerBp', 'hetHomRatio', 'tiTvRatio']
        info('Warning: no metrics for quality control, using '
             'default ' + ', '.join(qc_cnf['metircs']))

    to_exit = False
    dbs_dict = {}
    for db in qc_cnf['databases']:
        if not db:
            err('Empty field for quality_control databases')
            to_exit = True
        elif file_exists(db):
            if not verify_file(db, 'Vcf'):
                to_exit = True
            dbs_dict[basename(db)] = db
        elif db not in cnf['genome']:
            to_exit = True
            err(cnf.get('log'), db + ' for variant qc is not found '
                                     'in genome resources in system config.')
        else:
            dbs_dict[db] = cnf['genome'][db]

    if to_exit:
        exit()

    qc_cnf['database_vcfs'] = dbs_dict


def bcftools_qc(cnf, qc_dir, vcf_fpath):
    step_greetings(cnf, 'Quality control plots')

    work_dir = cnf['work_dir']

    bgzip = get_tool_cmdline(cnf, 'bgzip')
    tabix = get_tool_cmdline(cnf, 'tabix')
    bcftools = get_tool_cmdline(cnf, 'bcftools')
    plot_vcfstats = get_tool_cmdline(cnf, 'plot_vcfstats')
    if not bcftools or not tabix or not bcftools or not plot_vcfstats:
        exit()

    gzipped_fpath = join(work_dir, basename(vcf_fpath) + '.gz')
    cmdline = '{bgzip} -c {vcf_fpath}'.format(**locals())
    call(cnf, cmdline, None, gzipped_fpath)

    tbi_fpath = gzipped_fpath + '.tbi'
    cmdline = '{tabix} -f -p vcf {gzipped_fpath}'.format(**locals())
    call(cnf, cmdline, None, tbi_fpath)

    text_report_fpath = join(work_dir, cnf['name'] + '_bcftools.report')
    cmdline = '{bcftools} stats {gzipped_fpath}'.format(**locals())
    call(cnf, cmdline, None, text_report_fpath, to_remove=[gzipped_fpath, tbi_fpath])

    viz_report_dir = join(work_dir, cnf['name'] + '_qc_plots/')
    if file_exists(viz_report_dir):
        shutil.rmtree(viz_report_dir)
    mkdir(viz_report_dir)
    cmdline = '{plot_vcfstats} -s {text_report_fpath} -p {viz_report_dir} ' \
              '--no-PDF'.format(**locals())
    call(cnf, cmdline, text_report_fpath, None, output_is_file=False)
    return _get_plots_from_bcftools(viz_report_dir, qc_dir)


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


def _get_plots_from_bcftools(bcftools_report_dir, output_dir):
    original_plots_names = ['indels.0.png', 'substitutions.0.png']
    final_plots_names = ['indels.png', 'substitution.png']
    for i, original_plot in enumerate(original_plots_names):
        plot_src_filename = join(bcftools_report_dir, original_plot)
        plot_dst_filename = join(output_dir, final_plots_names[i])
        if file_exists(plot_src_filename):
            shutil.copyfile(plot_src_filename, plot_dst_filename)
    return [join(output_dir, fname) for fname in final_plots_names]