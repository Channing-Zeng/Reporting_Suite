from os.path import join
from source.bcbio_structure import BCBioStructure
from source.reporting import write_summary_reports, Record, SampleReport, FullReport
from source.logger import step_greetings, info
from source.targetcov.cov import get_basic_metrics, get_depth_metrics



def summary_reports(cnf, bcbio_structure):
    step_greetings('Coverage statistics for all samples')

    jsons_by_sample = bcbio_structure.get_targetcov_report_fpaths_by_sample('json')
    htmls_by_sample = bcbio_structure.get_targetcov_report_fpaths_by_sample('html')

    basic_metrics = get_basic_metrics()
    depth_metrics = get_depth_metrics(cnf.coverage_reports.depth_thresholds)

    basic_report, depth_report = [FullReport(report_name, [
        SampleReport(sample,
                     records=_get_parse_reports(metrics, jsons_by_sample[sample]),
                     html_fpath=htmls_by_sample[sample])
            for sample in bcbio_structure.samples
            if sample in jsons_by_sample and sample in htmls_by_sample])
        for metrics, report_name in [(basic_metrics, cnf.name),
                                     (depth_metrics, 'Target coverage depth')]]

    final_summary_report_fpaths = write_summary_reports(
        cnf.output_dir, cnf.work_dir, [basic_report, depth_report],
        cnf.name, 'Target coverage statistics')

    info()
    info('*' * 70)
    info('Summary:')
    for fpath in final_summary_report_fpaths:
        if fpath:
            info('  ' + fpath)


def _get_parse_reports(metrics, json_fpath):
    records = []
    for rec in Record.load_records(json_fpath):
        if rec.metric.name in metrics:
            rec.metric = metrics[rec.metric.name]
            records.append(rec)
    return records