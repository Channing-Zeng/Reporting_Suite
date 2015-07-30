import os
from os.path import join, relpath, dirname
from collections import OrderedDict
import getpass
from traceback import format_exc
import source

from source.bcbio_structure import BCBioStructure
from source.jira_utils import JiraCase
from source.logger import info, step_greetings, send_email, warn, err
from source.file_utils import verify_file, file_transaction, adjust_path, safe_mkdir, add_suffix
from source.reporting import Metric, Record, MetricStorage, ReportSection, SampleReport, FullReport
from source.html_reporting.html_saver import _write_static_html_report


def make_preproc_project_level_report(cnf, samples, project_name,
        project_dirpath, fastqc_dirpath, downsample_targqc_dirpath):
    step_greetings('Preproc project-level report')

    general_section = ReportSection('general_section', '', [])
    general_records = []

    # Summary
    for (name, repr_name, report_fpath) in [
            ('raw_fastqc_summary',            'Preproc FastQC',     join(fastqc_dirpath, 'FastQC.html')),
            ('downsampled_target_qc_summary', 'Seq QC downsamples', join(downsample_targqc_dirpath, 'targQC.html'))]:

        if verify_file(report_fpath):
            cur_metric = Metric(repr_name + ' summary', common=True)
            general_section.add_metric(cur_metric)
            general_records.append(
                Record(metric=cur_metric,
                       value=cur_metric.name,
                       html_fpath=_convert_to_relpath(
                           report_fpath,
                           project_dirpath)))

    # Individual reports
    individual_reports_section = ReportSection('individual_reports', '', [])
    sample_reports_records = dict()
    # _add_per_sample_reports(bcbio_structure, general_records, individual_reports_section)

    metric_storage = MetricStorage(general_section=general_section, sections=[individual_reports_section])
    sample_reports = []
    # for sample in samples:
    #     sample_reports.append(SampleReport(
    #         sample,
    #         records=sample_reports_records[sample.name],
    #         html_fpath=None,
    #         metric_storage=metric_storage))

    full_report = FullReport(cnf.project_name, sample_reports, metric_storage=metric_storage)
    html_fpath = join(project_dirpath, project_name + '.html')
    _save_static_html(cnf.work_dir, full_report, html_fpath, project_name)

    info()
    info('*' * 70)
    info('Project-level report saved in: ')
    info('  ' + html_fpath)
    return html_fpath


def make_postproc_project_level_report(cnf, bcbio_structure):
    step_greetings('Project-level report')

    general_section = ReportSection('general_section', '', [])
    general_records = []
    general_records = _add_summary_reports(bcbio_structure, general_section, general_records)
    general_records = _add_variants(bcbio_structure, general_section, general_records)

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

    full_report = FullReport(cnf.project_name, sample_reports, metric_storage=metric_storage)
    html_fpath = join(bcbio_structure.date_dirpath, bcbio_structure.project_name + '.html')
    _save_static_html(cnf.work_dir, full_report, html_fpath, bcbio_structure.project_name)

    info()
    info('*' * 70)
    info('Project-level report saved in: ')
    info('  ' + html_fpath)
    return html_fpath


def _add_variants(bcbio_structure, general_section, general_records):
    val = OrderedDict()
    single_val = OrderedDict()
    paired_val = OrderedDict()

    for caller in bcbio_structure.variant_callers.values():
        mut_fname = source.mut_fname_template.format(caller_name=caller.name)
        mut_fpath = join(bcbio_structure.date_dirpath, mut_fname)
        mut_pass_fpath = add_suffix(mut_fpath, source.mut_pass_suffix)
        if verify_file(mut_fpath, silent=True):
            val[caller.name] = mut_pass_fpath
        else:
            single_mut_fpath = add_suffix(mut_fpath, source.mut_single_suffix)
            paired_mut_fpath = add_suffix(mut_fpath, source.mut_paired_suffix)
            single_mut_pass_fpath = add_suffix(single_mut_fpath, source.mut_pass_suffix)
            paired_mut_pass_fpath = add_suffix(paired_mut_fpath, source.mut_pass_suffix)
            if verify_file(single_mut_pass_fpath, silent=True):
                single_val[caller.name] = single_mut_pass_fpath
                # _add_rec(single_mut_fpath, caller.name + ' mutations for separate samples')
            if verify_file(paired_mut_pass_fpath, silent=True):
                paired_val[caller.name] = paired_mut_pass_fpath
                # _add_rec(paired_mut_fpath, caller.name + ' mutations for paired samples')

    for val, metric_name in (
         (val, 'Mutations'),
         (single_val, 'Mutations for separate samples'),
         (paired_val, 'Mutations for paired samples')):
        if val:
            metric = Metric(metric_name, common=True)
            rec = Record(
                metric=metric,
                value=metric.name,
                html_fpath=_convert_to_relpath(val, bcbio_structure.date_dirpath))
            general_section.add_metric(metric)
            general_records.append(rec)

    return general_records


def _add_summary_reports(bcbio_structure, general_section, general_records):
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
    fastqc_htmls_by_sample      = dict([(s.name, verify_file(s.fastqc_html_fpath)) for s in bcbio_structure.samples])

    sample_reports_records = dict()
    for sample in bcbio_structure.samples:
        sample_reports_records[sample.name] = list(general_records)

    for (repr_name, links_by_sample) in [
            (bcbio_structure.fastqc_repr,      fastqc_htmls_by_sample),
            # ('BAM',                            bams_by_samples),
            (bcbio_structure.targqc_repr,      targqc_htmls_by_sample),
            (bcbio_structure.varqc_repr,       varqc_htmls_by_sample),
            (bcbio_structure.varqc_after_repr, varqc_after_htmls_by_sample)]:
        cur_metric = Metric(repr_name)
        individual_reports_section.add_metric(cur_metric)
        for sample in bcbio_structure.samples:
            if links_by_sample and links_by_sample.get(sample.name):
                sample_reports_records[sample.name].append(
                    Record(
                        metric=cur_metric,
                        value=cur_metric.name,
                        html_fpath=_convert_to_relpath(
                            links_by_sample[sample.name],
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
    targqc_htmls_by_sample = OrderedDict()

    for sample in bcbio_structure.samples:
        targqc_htmls_by_sample[sample.name] = OrderedDict()
        for report_name, report_html_fpath in [
            ('targetcov', sample.targetcov_html_fpath),
            ('ngscat', sample.ngscat_html_fpath),
            ('qualimap', sample.qualimap_html_fpath)]:
            if verify_file(report_html_fpath):
                targqc_htmls_by_sample[sample.name][report_name] = verify_file(report_html_fpath)

    return targqc_htmls_by_sample


def _convert_to_relpath(value, base_dirpath):
    if not value:
        return None
    if isinstance(value, str):
        return relpath(value, base_dirpath)
    elif isinstance(value, dict):
        for k in value.keys():
            if not value[k]:
                value[k] = None
            else:
                value[k] = relpath(value[k], base_dirpath)
        return value
    else:
        return value


def _save_static_html(work_dir, full_report, html_fpath, project_name):
    # metric name in FullReport --> metric name in Static HTML
    metric_names = OrderedDict([
        # (BCBioStructure.fastqc_repr + ' preproc', 'FastQC preproc'),
        (BCBioStructure.fastqc_repr, 'FastQC'),
        # ('BAM', 'BAM'),
        (BCBioStructure.targqc_repr, 'SeqQC'),
        # ('Downsampled ' + BCBioStructure.targqc_repr, 'Downsampled SeqQC'),
        # ('Mutations', 'Mutations'),
        # ('Mutations for separate samples', 'Mutations for separate samples'),
        # ('Mutations for paired samples', 'Mutations for paired samples'),
        (BCBioStructure.varqc_repr, 'VarQC'),
        (BCBioStructure.varqc_after_repr, 'VarQC after filtering')])

    def _process_record(record):
        new_html_fpath = []
        if isinstance(record["html_fpath"], basestring):
            new_html_fpath = [{"html_fpath_name": record["value"], "html_fpath_value": record["html_fpath"]}]
        elif isinstance(record["html_fpath"], dict):
            for k, v in record["html_fpath"].items():
                new_html_fpath.append({"html_fpath_name": k, "html_fpath_value": v})
        record["html_fpath"] = new_html_fpath
        return record

    def _get_summary_report_name(record):
        return record.value.lower().replace(' ', '_')

    # common records (summary reports)
    common_dict = dict()
    common_dict["project_name"] = project_name
    main_dict = dict()
    if full_report.sample_reports:
        sample_report = full_report.sample_reports[0]
        for record in sample_report.records:
            if record.metric.common:
                common_dict[_get_summary_report_name(record)] = record.__dict__
                if 'Mutations' in record.metric.name:
                    common_dict[_get_summary_report_name(record)]['contents'] = (
                        record.metric.name + ': ' + ', '.join('<a href={v}>{k}</a>'.format(k=k, v=v) for k, v in record.html_fpath.items()))

        # individual records
        main_dict = dict()
        main_dict["sample_reports"] = []
        main_dict["metric_names"] = [v for k, v in metric_names.items()
                                     if k in [m.name for m in full_report.metric_storage.get_metrics(skip_general_section=True)]]
        for sample_report in full_report.sample_reports:
            new_records = [_process_record(record.__dict__) for record in sample_report.records
                           if record.metric.name in metric_names.keys()]
            sample_report_dict = dict()
            sample_report_dict["records"] = new_records
            sample_report_dict["sample_name"] = sample_report.get_display_name()
            main_dict["sample_reports"].append(sample_report_dict)

    return _write_static_html_report(work_dir, {"common": common_dict, "main": main_dict}, html_fpath)