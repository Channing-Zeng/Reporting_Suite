from os.path import join, relpath
from collections import OrderedDict
from source.bcbio_structure import BCBioStructure

from source.logger import info, step_greetings, send_email
from source.file_utils import verify_file
from source.reporting import Metric, Record, MetricStorage, ReportSection, SampleReport, FullReport


def make_project_level_report(cnf, bcbio_structure):
    step_greetings('Project-level report')

    general_section = ReportSection('general_section', '', [])
    general_records = _add_summary_reports(bcbio_structure, general_section)

    individual_reports_section = ReportSection('individual_reports', '', [])
    sample_reports_records = _add_per_sample_reports(bcbio_structure, general_records, individual_reports_section)

    metric_storage = MetricStorage(general_section=general_section, sections=[individual_reports_section])
    sample_reports = []
    for sample in bcbio_structure.samples:
        sample_reports.append(SampleReport(
            sample,
            records=sample_reports_records[sample.name],
            html_fpath=None,
            metric_storage=metric_storage))

    full_report = FullReport(cnf.name, sample_reports, metric_storage=metric_storage)
    final_summary_report_fpath = full_report.save_html(
        bcbio_structure.date_dirpath, bcbio_structure.project_name,
        'Project-level report for ' + bcbio_structure.project_name)

    info()
    info('*' * 70)
    info('Project-level report saved in: ')
    info('  ' + final_summary_report_fpath)
    send_email('Report for ' + bcbio_structure.project_name + ':'
               '\n  ' + final_summary_report_fpath)


def _add_summary_reports(bcbio_structure, general_section):
    general_records = []

    for (name, repr_name, summary_dir) in [
            (BCBioStructure.fastqc_name,      BCBioStructure.fastqc_repr,      BCBioStructure.fastqc_summary_dir),
            (BCBioStructure.targqc_name,      BCBioStructure.targqc_repr,      BCBioStructure.targqc_summary_dir),
            (BCBioStructure.varqc_name,       BCBioStructure.varqc_repr,       BCBioStructure.varqc_summary_dir),
            (BCBioStructure.varqc_after_name, BCBioStructure.varqc_after_repr, BCBioStructure.varqc_after_summary_dir)]:

        summary_report_fpath = join(bcbio_structure.date_dirpath, summary_dir, name + '.html')
        if verify_file(summary_report_fpath):
            cur_metric = Metric(repr_name + ' summary', common=True)
            general_section.add_metric(cur_metric)
            general_records.append(
                Record(metric=cur_metric,
                       value=cur_metric.name,
                       html_fpath=_convert_to_relpath(
                           summary_report_fpath,
                           bcbio_structure.date_dirpath)))
    return general_records


def _add_per_sample_reports(bcbio_structure, general_records, individual_reports_section):
    varqc_htmls_by_sample       = _add_varqc_reports(bcbio_structure, BCBioStructure.varqc_name, BCBioStructure.varqc_dir)
    varqc_after_htmls_by_sample = _add_varqc_reports(bcbio_structure, BCBioStructure.varqc_after_name, BCBioStructure.varqc_after_dir)
    targqc_htmls_by_sample      = _add_targqc_reports(bcbio_structure)
    fastqc_htmls_by_sample      = bcbio_structure.get_fastqc_report_fpaths_by_sample()

    sample_reports_records = dict()
    for sample in bcbio_structure.samples:
        sample_reports_records[sample.name] = list(general_records)

    for (repr_name, htmls_by_sample) in [
            (bcbio_structure.fastqc_repr,      fastqc_htmls_by_sample),
            (bcbio_structure.targqc_repr,      targqc_htmls_by_sample),
            (bcbio_structure.varqc_repr,       varqc_htmls_by_sample),
            (bcbio_structure.varqc_after_repr, varqc_after_htmls_by_sample)]:
        if htmls_by_sample:
            cur_metric = Metric(repr_name)
            individual_reports_section.add_metric(cur_metric)
            for sample in bcbio_structure.samples:
                if sample.name in htmls_by_sample:
                    sample_reports_records[sample.name].append(Record(
                        metric=cur_metric,
                        value=cur_metric.name,
                        html_fpath=_convert_to_relpath(
                            htmls_by_sample[sample.name],
                            bcbio_structure.date_dirpath)))
                else:
                    sample_reports_records[sample.name].append(
                        Record(metric=cur_metric, value=None, html_fpath=None))
    return sample_reports_records


def _add_varqc_reports(bcbio_structure, name, dir_name):
    callers = bcbio_structure.variant_callers.values()
    if len(callers) == 0:
        varqc_htmls_by_sample = None
    elif len(callers) == 1:
        varqc_htmls_by_sample = callers[0].find_fpaths_by_sample(
            dir_name, name, 'html')
    else:
        varqc_htmls_by_sample = OrderedDict()
        for sample in bcbio_structure.samples:
            varqc_htmls_by_sample[sample.name] = OrderedDict()
        for caller in callers:
            for sample, fpath in caller.find_fpaths_by_sample(
                    dir_name, name, 'html').items():
                varqc_htmls_by_sample[sample][caller.name] = fpath

    return varqc_htmls_by_sample


def _add_targqc_reports(bcbio_structure):
    targetcov_htmls_by_sample = bcbio_structure.find_targetcov_reports_by_sample('html')
    ngscat_htmls_by_sample = bcbio_structure.find_ngscat_reports_by_sample()
    qualimap_htmls_by_sample = bcbio_structure.find_qualimap_reports_by_sample()

    targqc_htmls_by_sample = OrderedDict()

    for sample in bcbio_structure.samples:
        targqc_htmls_by_sample[sample.name] = OrderedDict()

        if sample.name in targetcov_htmls_by_sample:
            targqc_htmls_by_sample[sample.name]['targetcov'] = targetcov_htmls_by_sample[sample.name]
        if sample.name in ngscat_htmls_by_sample:
            targqc_htmls_by_sample[sample.name]['ngscat'] = ngscat_htmls_by_sample[sample.name]
        if sample.name in qualimap_htmls_by_sample:
            targqc_htmls_by_sample[sample.name]['qualimap'] = qualimap_htmls_by_sample[sample.name]

    return targqc_htmls_by_sample


def _convert_to_relpath(value, base_dirpath):
    if isinstance(value, str):
        return relpath(value, base_dirpath)
    elif isinstance(value, dict):
        for k in value.keys():
            value[k] = relpath(value[k], base_dirpath)
        return value
    else:
        return value