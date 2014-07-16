from source.logger import info
from source.reporting import summarize, write_summary_reports, get_sample_report_fpaths_for_bcbio_final_dir
database = 'cosmic'
main_novelty = 'all'
metrics_header = 'Metric'
novelty_header = 'Novelty'
sample_header = 'Sample name:'


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
    full_report = summarize(sample_names, callers[0].single_qc_rep_fpaths, parse_qc_sample_report)

    full_summary_fpaths = write_summary_reports(
        cnf['output_dir'], cnf['work_dir'], full_report,
        sample_names, 'varqc.summary', 'Variant QC')

    info()
    info('*' * 70)
    for fpath in full_summary_fpaths:
        info(fpath)


def _make_for_multiple_variant_callers(callers, cnf, sample_names):
    for caller in callers:
        caller.summary_qc_report = summarize(
            sample_names, caller.single_qc_rep_fpaths, parse_qc_sample_report)

        caller.summary_qc_rep_fpaths = write_summary_reports(
            cnf['output_dir'], cnf['work_dir'], caller.summary_qc_report,
            sample_names, caller.suf + '.varqc.summary', 'Variant QC for ' + caller.name)

    all_single_reports = [r for c in callers for r in c.single_qc_rep_fpaths]
    all_sample_names = [sample_name + '-' + c.suf for sample_name in sample_names for c in callers]

    full_summary_report = summarize(all_sample_names, all_single_reports, parse_qc_sample_report)

    full_summary_fpaths = write_summary_reports(
        cnf['output_dir'], cnf['work_dir'], full_summary_report,
        all_sample_names, 'varqc.summary', 'Variant QC')

    info()
    info('*' * 70)

    for caller in callers:
        info(caller.name)
        for fpath in caller.summary_qc_rep_fpaths:
            info('  ' + fpath)
        info()

    for fpath in full_summary_fpaths:
        info('  ' + fpath)


def parse_qc_sample_report(report_fpath):
    """ returns row_per_sample =
            dict(metricName=None, value=None,
            isMain=True, quality='More is better')
    """
    row_per_sample = []

    with open(report_fpath) as f:
        # parsing Sample name and Database columns
        database_col_id = None
        novelty_col_id = None
        for line in f:
            if not line.strip():
                continue
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
                is_main = True
                novelty = line.split()[novelty_col_id]
                if novelty_col_id and novelty != main_novelty:
                    is_main = False

                cur_metric_name = line.split()[0]
                cur_value = line.split()[database_col_id]

                row_per_sample.append(dict(
                    metricName=cur_metric_name, value=cur_value,
                    isMain=is_main, quality='More is better'))

    return row_per_sample