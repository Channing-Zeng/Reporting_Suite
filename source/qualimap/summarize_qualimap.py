from os.path import relpath
from source.reporting import FullReport, SampleReport
from source.logger import step_greetings, info, send_email
from source.qualimap import report_parser


def summary_reports(cnf, bcbio_structure):
    step_greetings('QualiMap statistics for all samples')

    htmls_by_sample = bcbio_structure.get_qualimap_report_fpaths_by_sample()
    if not htmls_by_sample:
        return []

    full_report = FullReport(cnf.name, [
        SampleReport(sample,
                     records=report_parser.parse_qualimap_sample_report(htmls_by_sample[sample.name]),
                     html_fpath=relpath(htmls_by_sample[sample.name], cnf.output_dir))
            for sample in bcbio_structure.samples
            if sample.name in htmls_by_sample],
        metric_storage=report_parser.metric_storage)

    final_summary_report_fpaths = full_report.save_into_files(
        cnf.output_dir, bcbio_structure.qualimap_name, 'QualiMap statistics')

    info()
    info('*' * 70)
    info('Summary:')
    for fpath in final_summary_report_fpaths:
        if fpath: info('  ' + fpath)
    # send_email('Qualimap summary: \n' + '\n'.join(final_summary_report_fpaths))

    return final_summary_report_fpaths
