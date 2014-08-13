from os.path import join
from source.bcbio_structure import BCBioStructure
from source.reporting import summarize, write_summary_reports, Record
from source.logger import step_greetings, info
from source.targetcov.cov import get_cov_metrics


def summary_reports(cnf, bcbio_structure):
    step_greetings('Coverage statistics for all samples')

    json_by_sample = bcbio_structure.get_targetcov_json_by_sample()
    sum_report = summarize(cnf, json_by_sample, _parse_targetseq_sample_report)

    final_summary_report_fpaths = write_summary_reports(
        cnf.output_dir,
        cnf.work_dir,
        sum_report,
        join(cnf.output_dir, BCBioStructure.targetseq_name),
        'Target coverage statistics')

    info()
    info('*' * 70)
    info('Summary:')
    for fpath in final_summary_report_fpaths:
        if fpath:
            info('  ' + fpath)


def _parse_targetseq_sample_report(cnf, json_fpath):
    records = Record.load_records(json_fpath)
    cov_metrics = get_cov_metrics(cnf.coverage_reports.depth_thresholds)
    for rec in records:
        if rec.metric.name in cov_metrics:
            rec.metric = cov_metrics[rec.metric.name]

    return records

