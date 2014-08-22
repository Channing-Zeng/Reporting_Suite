from os.path import join
from source.bcbio_structure import BCBioStructure
from source.reporting import summarize, write_summary_reports, Record
from source.logger import step_greetings, info
from source.targetcov.cov import get_basic_metrics, get_depth_metrics


def summary_reports(cnf, bcbio_structure):
    step_greetings('Coverage statistics for all samples')

    report_fpath_by_sample = bcbio_structure.get_targetcov_report_fpaths_by_sample()

    basic_metrics = get_basic_metrics()
    depth_metrics = get_depth_metrics(cnf.coverage_reports.depth_thresholds)

    basic_report = summarize(cnf, report_fpath_by_sample, _get_parse_reports(basic_metrics), '')
    depth_report = summarize(cnf, report_fpath_by_sample, _get_parse_reports(depth_metrics), 'Target coverage depth')

    final_summary_report_fpaths = write_summary_reports(
        cnf.output_dir,
        cnf.work_dir,
        [basic_report, depth_report],
        join(cnf.output_dir, cnf.name),
        'Target coverage statistics')

    info()
    info('*' * 70)
    info('Summary:')
    for fpath in final_summary_report_fpaths:
        if fpath:
            info('  ' + fpath)


def _get_parse_reports(metrics):
    def _parse_targetseq_report(cnf, json_fpath):
        _records = Record.load_records(json_fpath)
        records = []
        for rec in _records:
            if rec.metric.name in metrics:
                rec.metric = metrics[rec.metric.name]
                records.append(rec)
        return records

    return _parse_targetseq_report