from collections import OrderedDict
from source.logger import info
from source.reporting import summarize, write_summary_reports, get_sample_report_fpaths_for_bcbio_final_dir, Metric
from source.utils import OrderedDefaultDict

main_novelty = 'all'
metrics_header = 'Metric'
novelty_header = 'Novelty'
average_header = 'Average'


class VariantCaller:
    def __init__(self, suf):
        self.name = suf
        self.suf = suf
        self.single_qc_rep_fpaths = []
        self.summary_qc_report = None
        self.summary_qc_rep_fpaths = []


def make_summary_reports(cnf, sample_names):
    varqc_dir = cnf['base_name']

    vcf_sufs = cnf['vcf_suf'].split(',')
    callers = [VariantCaller(suf) for suf in vcf_sufs]

    for caller in callers:
        fpaths, sample_names = get_sample_report_fpaths_for_bcbio_final_dir(
            cnf['bcbio_final_dir'], sample_names, varqc_dir,
            '-' + caller.suf + '.varqc.txt')
        if fpaths:
            caller.single_qc_rep_fpaths = fpaths

    if len(callers) > 1:
        _make_for_multiple_variant_callers(callers, cnf, sample_names)

    else:
        _make_for_single_variant_caller(callers, cnf, sample_names)


def _make_for_single_variant_caller(callers, cnf, sample_names):
    full_report = summarize(sample_names, callers[0].single_qc_rep_fpaths, get_parse_qc_sample_report(cnf))

    full_summary_fpaths = write_summary_reports(
        cnf['output_dir'], cnf['work_dir'], full_report, 'varqc.summary', 'Variant QC')

    info()
    info('*' * 70)
    for fpath in full_summary_fpaths:
        info(fpath)


def _make_for_multiple_variant_callers(callers, cnf, sample_names):
    for caller in callers:
        caller.summary_qc_report = summarize(
            sample_names, caller.single_qc_rep_fpaths, get_parse_qc_sample_report(cnf))

        caller.summary_qc_rep_fpaths = write_summary_reports(
            cnf['output_dir'], cnf['work_dir'], caller.summary_qc_report,
            caller.suf + '.varqc.summary', 'Variant QC for ' + caller.name)

    all_single_reports = [r for c in callers for r in c.single_qc_rep_fpaths]
    all_sample_names = [sample_name + '-' + c.suf for sample_name in sample_names for c in callers]

    full_summary_report = summarize(all_sample_names, all_single_reports, get_parse_qc_sample_report(cnf))

    full_summary_fpaths = write_summary_reports(
        cnf['output_dir'], cnf['work_dir'], full_summary_report, 'varqc.summary', 'Variant QC')

    info()
    info('*' * 70)

    for caller in callers:
        info(caller.name)
        for fpath in caller.summary_qc_rep_fpaths:
            info('  ' + fpath)
        info()

    info('Total')
    for fpath in full_summary_fpaths:
        info('  ' + fpath)


def get_parse_qc_sample_report(cnf):
    def _parse_qc_sample_report(report_fpath):
        """ returns row_per_sample =
                dict(metricName=None, value=None,
                isMain=True, quality='More is better')
        """

        metrics = OrderedDefaultDict(Metric)
        rest_headers = []
                # metrics[metric_name]['meta']
        with open(report_fpath) as f:
            # parsing Sample name and Database columns
            main_value_col_id = None
            novelty_col_id = None
            for line in f:
                if not line.strip():
                    continue

                elif line.startswith(metrics_header):
                    if cnf.quality_control.db_for_summary in line:
                        main_value_col_id = line.split().index(cnf.quality_control.db_for_summary)

                    if novelty_header in line:
                        novelty_col_id = line.split().index(novelty_header)

                    rest_headers = line.split()[2:]

                elif novelty_col_id:
                    # parsing rest of the report
                    metric_name = line.split()[0]
                    novelty = line.split()[novelty_col_id]

                    metrics[metric_name].name = metric_name
                    metrics[metric_name].quality = 'More is better'
                    metrics[metric_name].meta[novelty] = dict(zip(rest_headers, line.split()[2:]))
                    if novelty == main_novelty:
                        metrics[metric_name].value = line.split()[main_value_col_id]

        return metrics

    return _parse_qc_sample_report

