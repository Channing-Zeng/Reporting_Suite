from collections import OrderedDict
from os.path import relpath
from source.reporting import SampleReport, FullReport, Metric, MetricStorage, ReportSection, Record, load_records
from source.logger import step_greetings, info, send_email, critical
from source.targetcov import cov
from source.qualimap import report_parser as qualimap_report_parser
from source.ngscat import report_parser as ngscat_report_parser


qualimap_to_targetcov_dict = {'Number of reads': cov.header_metric_storage.get_metric('Reads'),
                              'Mapped reads': cov.header_metric_storage.get_metric('Mapped reads'),
                              'Unmapped reads': cov.header_metric_storage.get_metric('Unmapped reads'),
                              'Mapped reads (on target)': cov.header_metric_storage.get_metric('Reads mapped on target'),
                              'Coverage Mean': cov.header_metric_storage.get_metric('Average target coverage depth'),
                              'Coverage Standard Deviation': cov.header_metric_storage.get_metric('Std. dev. of target coverage depth')}

ngscat_to_targetcov_dict = {'Number reads': cov.header_metric_storage.get_metric('Mapped reads'),
                            '% target bases with coverage >= 1x': cov.header_metric_storage.get_metric('Percentage of target covered by at least 1 read'),
                            '% reads on target': cov.header_metric_storage.get_metric('Reads mapped on target'),
                            'mean coverage': cov.header_metric_storage.get_metric('Average target coverage depth')}


def _get_targqc_metric(metric, report_type='targetcov'):  # report type is in ['targetcov', 'qualimap', 'ngscat']
    if report_type == 'targetcov':
        return metric
    elif report_type == 'qualimap':
        if metric.name in qualimap_to_targetcov_dict.keys():
            return qualimap_to_targetcov_dict[metric.name]
        return metric
    elif report_type == 'ngscat':
        if metric.name in ngscat_to_targetcov_dict.keys():
            return ngscat_to_targetcov_dict[metric.name]
        return metric
    critical('Incorrect usage of get_targqc_metric(), report_type is %s but should be one of the following: %s' %
             (report_type, ", ".join(['targetcov', 'qualimap', 'ngscat'])))
    return None


def _get_targqc_metric_storage(metric_storages_by_report_type):
    metric_list = []
    general_section_metric_list = []
    for report_type, metric_storage in metric_storages_by_report_type.items():
        # Note: skipping metrics with equivalents in TargetCov metrics
        general_section_metric_list += [metric for metric in metric_storage.general_section.metrics
                                        if metric == _get_targqc_metric(metric, report_type)]
        metric_list += [metric for metric in metric_storage.get_metrics(skip_general_section=True)
                        if metric == _get_targqc_metric(metric, report_type)]
    return MetricStorage(general_section=ReportSection('general_section', '', general_section_metric_list),
                         sections=[ReportSection('basic', '', metric_list)])


def _get_targqc_records(records_by_report_type):
    targqc_records = []
    filled_metric_names = []
    for report_type, records in records_by_report_type.items():
        for record in records:
            new_metric = _get_targqc_metric(record.metric, report_type)
            if new_metric.name not in filled_metric_names:
                filled_metric_names.append(new_metric.name)
                record.metric = new_metric
                targqc_records.append(record)
    return targqc_records


def summary_reports(cnf, bcbio_structure):
    step_greetings('Coverage statistics for all samples based on TargetSeq, ngsCAT, and Qualimap reports')

    targetcov_metric_storage = cov.header_metric_storage
    for depth in cnf.coverage_reports.depth_thresholds:
        name = 'Part of target covered at least by ' + str(depth) + 'x'
        targetcov_metric_storage.add_metric(
            Metric(name, short_name=str(depth) + 'x', description=name, unit='%'),
            'depth_metrics')
    targetcov_jsons_by_sample = bcbio_structure.find_targetcov_reports_by_sample('json')
    targetcov_htmls_by_sample = bcbio_structure.find_targetcov_reports_by_sample('html')
    ngscat_htmls_by_sample = bcbio_structure.get_ngscat_report_fpaths_by_sample()
    qualimap_htmls_by_sample = bcbio_structure.get_qualimap_report_fpaths_by_sample()

    all_htmls_by_sample = OrderedDict()
    for sample in bcbio_structure.samples:
        all_htmls_by_sample[sample.name] = OrderedDict()
        if sample.name in targetcov_htmls_by_sample:
            all_htmls_by_sample[sample.name]['targetcov'] = relpath(targetcov_htmls_by_sample[sample.name], cnf.output_dir)
        if sample.name in ngscat_htmls_by_sample:
            all_htmls_by_sample[sample.name]['ngscat'] = relpath(ngscat_htmls_by_sample[sample.name], cnf.output_dir)
        if sample.name in qualimap_htmls_by_sample:
            all_htmls_by_sample[sample.name]['qualimap'] = relpath(qualimap_htmls_by_sample[sample.name], cnf.output_dir)

    targqc_metric_storage = _get_targqc_metric_storage(OrderedDict(
        targetcov=targetcov_metric_storage,
        ngscat=ngscat_report_parser.metric_storage,
        qualimap=qualimap_report_parser.metric_storage))

    targqc_full_report = FullReport(cnf.name, [
        SampleReport(sample,
                     records=_get_targqc_records(OrderedDict(
                         targetcov=load_records(targetcov_jsons_by_sample[sample.name])
                         if sample.name in targetcov_jsons_by_sample else [],
                         ngscat=ngscat_report_parser.parse_ngscat_sample_report(ngscat_htmls_by_sample[sample.name])
                         if sample.name in ngscat_htmls_by_sample else [],
                         qualimap=qualimap_report_parser.parse_qualimap_sample_report(qualimap_htmls_by_sample[sample.name])
                         if sample.name in qualimap_htmls_by_sample else []
                     )),
                     html_fpath=all_htmls_by_sample[sample.name])
            for sample in bcbio_structure.samples
            if sample.name in set(targetcov_jsons_by_sample.keys() +
                                  ngscat_htmls_by_sample.keys() +
                                  qualimap_htmls_by_sample.keys())],
        metric_storage=targqc_metric_storage)

    final_summary_report_fpath = targqc_full_report.save_html(
        cnf.output_dir, bcbio_structure.targqc_name,
        'Coverage statistics for all samples based on TargetSeq, ngsCAT, and Qualimap reports')

    info()
    info('*' * 70)
    info('TargQC summary saved in: ')
    info('  ' + final_summary_report_fpath)
    send_email('TargQC summary: \n' + final_summary_report_fpath)
