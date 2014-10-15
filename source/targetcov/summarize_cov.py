from source.reporting import SampleReport, FullReport, Metric
from source.logger import step_greetings, info, send_email
from source.targetcov import cov


def summary_reports(cnf, bcbio_structure):
    step_greetings('Coverage statistics for all samples')

    metric_storage = cov.header_metric_storage  # TODO parse storage from Json too
    for depth in cnf.coverage_reports.depth_thresholds:
        name = 'Part of target covered at least by ' + str(depth) + 'x'
        metric_storage.add_metric(
            Metric(name, short_name=str(depth) + 'x', description=name, unit='%'),
            'depth_metrics')

    jsons_by_sample = bcbio_structure.find_targetcov_reports_by_sample('json')
    htmls_by_sample = bcbio_structure.find_targetcov_reports_by_sample('html')

    full_report = FullReport.construct_from_sample_report_jsons(
        bcbio_structure.samples, bcbio_structure,
        jsons_by_sample, htmls_by_sample, cnf.output_dir)

    final_summary_report_fpaths = full_report.save_into_files(
        cnf.output_dir, cnf.name, 'Target coverage statistics')

    info()
    info('*' * 70)
    info('Summary:')
    for fpath in final_summary_report_fpaths:
        if fpath:
            info('  ' + fpath)

    # send_email('TargetSeq summary: \n' + '\n'.join(final_summary_report_fpaths))
