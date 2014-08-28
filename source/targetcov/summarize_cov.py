from os.path import join
from source.bcbio_structure import BCBioStructure
from source.reporting import write_summary_report, Record, SampleReport, FullReport, Metric
from source.logger import step_greetings, info
from source.targetcov import cov


def summary_reports(cnf, bcbio_structure):
    step_greetings('Coverage statistics for all samples')

    jsons_by_sample = bcbio_structure.get_targetcov_report_fpaths_by_sample('json')
    htmls_by_sample = bcbio_structure.get_targetcov_report_fpaths_by_sample('html')

    metric_storage = cov.metric_storage  # TODO parse stroage from Json too
    for depth in cnf.coverage_reports.depth_thresholds:
        name = 'Part of target covered at least by ' + str(depth) + 'x'
        metric_storage.add_metric(
            Metric(name, short_name=str(depth) + 'x', description=name, unit='%'),
            'depth_metrics')

    full_report = FullReport(sample_reports=[
            SampleReport(
                sample,
                metric_storage=metric_storage,
                records=_parse_reports(jsons_by_sample[sample], metric_storage),
                html_fpath=htmls_by_sample[sample])
            for sample in bcbio_structure.samples
            if sample in jsons_by_sample and sample in htmls_by_sample])

    final_summary_report_fpaths = write_summary_report(
        cnf.output_dir, cnf.work_dir, full_report,
        cnf.name, 'Target coverage statistics')

    info()
    info('*' * 70)
    info('Summary:')
    for fpath in final_summary_report_fpaths:
        if fpath:
            info('  ' + fpath)


def _parse_reports(json_fpath, metric_storage):
    records = []
    for rec in Record.load_records(json_fpath):
        records.append(rec)
        rec.metric = metric_storage.get_metric(rec.metric.name)
    return records