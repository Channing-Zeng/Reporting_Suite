import os
from os.path import join, relpath, dirname
from collections import OrderedDict, namedtuple
import getpass
from traceback import format_exc
from pysam.cvcf import defaultdict

import source
from source.preproc.dataset_structure import DatasetStructure
from source.bcbio_structure import BCBioStructure
from source.jira_utils import JiraCase
from source.logger import info, step_greetings, send_email, warn, err
from source.file_utils import verify_file, file_transaction, adjust_path, safe_mkdir, add_suffix
from source.reporting import Metric, Record, MetricStorage, ReportSection, SampleReport, FullReport
from source.html_reporting.html_saver import _write_static_html_report


def make_preproc_project_level_report(cnf, dataset_structure):
    step_greetings('Preproc project-level report')

    general_section = ReportSection('general_section', '', [])
    general_records = []
    general_records = _add_summary_reports(general_section, general_records, dataset_structure=dataset_structure)

    individual_reports_section = ReportSection('individual_reports', '', [])
    sample_reports_records = _add_per_sample_reports(general_records, individual_reports_section, dataset_structure=dataset_structure)

    metric_storage = MetricStorage(general_section=general_section, sections=[individual_reports_section])
    sample_reports = []
    for sample in dataset_structure.samples:
        sample_reports.append(SampleReport(
            sample,
            records=sample_reports_records[sample.name],
            html_fpath=None,
            metric_storage=metric_storage))

    full_report = FullReport(cnf.project_name, sample_reports, metric_storage=metric_storage, general_records=general_records)
    _save_static_html(cnf.work_dir, full_report, dataset_structure.project_report_html_fpath, dataset_structure.project_name)

    info()
    info('*' * 70)
    info('Project-level report saved in: ')
    info('  ' + dataset_structure.project_report_html_fpath)
    return dataset_structure.project_report_html_fpath


def make_postproc_project_level_report(cnf, bcbio_structure):
    step_greetings('Project-level report')

    general_section = ReportSection('general_section', '', [])
    general_records = []
    general_records = _add_summary_reports(general_section, general_records, bcbio_structure=bcbio_structure)
    general_records = _add_variants(bcbio_structure, general_section, general_records)

    individual_reports_section = ReportSection('individual_reports', '', [])
    sample_reports_records = _add_per_sample_reports(general_records, individual_reports_section, bcbio_structure=bcbio_structure)

    metric_storage = MetricStorage(general_section=general_section, sections=[individual_reports_section])
    sample_reports = []
    for sample in bcbio_structure.samples:
        sample_reports.append(SampleReport(
            sample,
            records=sample_reports_records[sample.name],
            html_fpath=None,
            metric_storage=metric_storage))

    full_report = FullReport(cnf.project_name, sample_reports, metric_storage=metric_storage)
    _save_static_html(cnf.work_dir, full_report, bcbio_structure.project_level_report_fpath, bcbio_structure.project_name)

    info()
    info('*' * 70)
    info('Project-level report saved in: ')
    info('  ' + bcbio_structure.project_level_report_fpath)
    return bcbio_structure.project_level_report_fpath


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


def _add_summary_reports(general_section, general_records, bcbio_structure=None, dataset_structure=None):
    """ We want links to be relative, so we make paths relative to the project-level-report parent directory.
        - If the bcbio_structure is set, project-level report is located at bcbio_structure.date_dirpath
        - If dataset_dirpath is set, project-level report is located right at dataset_dirpath
    """
    SummaryMetric = namedtuple('SummaryMetric', 'name report_fpath')
    summary_metrics = []
    base_dirpath = None

    if dataset_structure:
        summary_metrics.append(SummaryMetric(
            name=DatasetStructure.pre_fastqc_repr,
            report_fpath=dataset_structure.comb_fastqc_fpath))
        summary_metrics.append(SummaryMetric(
            name=DatasetStructure.downsample_targqc_repr,
            report_fpath=dataset_structure.downsample_targqc_report_fpath))
        base_dirpath = dirname(dataset_structure.project_report_html_fpath)
            # ('raw_fastqc_summary',            ,     ),
            # ('downsampled_target_qc_summary', , join(downsample_targqc_dirpath, 'targQC.html'))]:

    if bcbio_structure:
        for (name, repr_name, summary_dir) in [
            (BCBioStructure.fastqc_name,      BCBioStructure.fastqc_repr,      BCBioStructure.fastqc_summary_dir),
            (BCBioStructure.targqc_name,      BCBioStructure.targqc_repr,      BCBioStructure.targqc_summary_dir),
            (BCBioStructure.varqc_name,       BCBioStructure.varqc_repr,       BCBioStructure.varqc_summary_dir),
            (BCBioStructure.varqc_after_name, BCBioStructure.varqc_after_repr, BCBioStructure.varqc_after_summary_dir)]:

            summary_report_fpath = join(bcbio_structure.date_dirpath, summary_dir, name + '.html')
            summary_metrics.append(SummaryMetric(
                name=repr_name + ' summary',
                report_fpath=summary_report_fpath))
        base_dirpath = bcbio_structure.date_dirpath

    for sm in summary_metrics:
        if verify_file(sm.report_fpath):
            metric = Metric(sm.name, common=True)
            general_section.add_metric(metric)
            general_records.append(Record(metric=metric, value=metric.name, html_fpath=relpath(sm.report_fpath, base_dirpath)))

    return general_records


def _add_per_sample_reports(general_records, individual_reports_section, bcbio_structure=None, dataset_structure=None):
    to_add = []
    base_dirpath = None

    if dataset_structure:
        pre_fastqc_htmls_by_sample = dict([(s.name, verify_file(s.fastqc_html_fpath)) for s in dataset_structure.samples])
        targqc_htmls_by_sample     = _add_targqc_reports(dataset_structure.samples)
        to_add.extend([
            (DatasetStructure.pre_fastqc_repr,        pre_fastqc_htmls_by_sample),
            (DatasetStructure.downsample_targqc_repr, targqc_htmls_by_sample),
        ])
        base_dirpath = dirname(dataset_structure.project_report_html_fpath)

    if bcbio_structure:
        varqc_htmls_by_sample       = _add_varqc_reports(bcbio_structure, BCBioStructure.varqc_name, BCBioStructure.varqc_dir)
        varqc_after_htmls_by_sample = _add_varqc_reports(bcbio_structure, BCBioStructure.varqc_after_name, BCBioStructure.varqc_after_dir)
        targqc_htmls_by_sample      = _add_targqc_reports(bcbio_structure.samples)
        fastqc_htmls_by_sample      = dict([(s.name, verify_file(s.fastqc_html_fpath)) for s in bcbio_structure.samples])
        to_add.extend([
            (BCBioStructure.fastqc_repr,      fastqc_htmls_by_sample),
            # ('BAM',                            bams_by_samples),
            (BCBioStructure.targqc_repr,      targqc_htmls_by_sample),
            (BCBioStructure.varqc_repr,       varqc_htmls_by_sample),
            (BCBioStructure.varqc_after_repr, varqc_after_htmls_by_sample)
        ])
        base_dirpath = dirname(bcbio_structure.project_level_report_fpath)

    sample_reports_records = defaultdict(list)

    for (repr_name, links_by_sample) in to_add:
        cur_metric = Metric(repr_name)
        individual_reports_section.add_metric(cur_metric)

        samples = []
        if dataset_structure:
            samples.extend(dataset_structure.samples)
        if bcbio_structure:
            samples.extend(bcbio_structure.samples)

        for sample in samples:
            if links_by_sample and links_by_sample.get(sample.name):
                sample_reports_records[sample.name].append(
                    Record(
                        metric=cur_metric,
                        value=cur_metric.name,
                        html_fpath=_convert_to_relpath(links_by_sample[sample.name], base_dirpath)))
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


def _add_targqc_reports(samples):
    targqc_htmls_by_sample = OrderedDict()

    for sample in samples:
        targqc_htmls_by_sample[sample.name] = OrderedDict()
        for report_name, report_html_fpath in [
            ('targetcov', sample.targetcov_html_fpath),
            ('ngscat',    sample.ngscat_html_fpath),
            ('qualimap',  sample.qualimap_html_fpath)]:
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

    def _process_record(rec):
        new_html_fpath = []
        if isinstance(rec["html_fpath"], basestring):
            new_html_fpath = [{"html_fpath_name": rec["value"], "html_fpath_value": rec["html_fpath"]}]
        elif isinstance(rec["html_fpath"], dict):
            for k, v in rec["html_fpath"].items():
                new_html_fpath.append({"html_fpath_name": k, "html_fpath_value": v})
        rec["html_fpath"] = new_html_fpath
        return rec

    def _get_summary_report_name(rec):
        return rec.value.lower().replace(' ', '_')

    # common records (summary reports)
    common_dict = dict()
    common_dict["project_name"] = project_name
    common_records = full_report.get_common_records()
    for record in common_records:
        common_dict[_get_summary_report_name(record)] = record.__dict__
        if 'Mutations' in record.metric.name:
            common_dict[_get_summary_report_name(record)]['contents'] = (
                record.metric.name + ': ' + ', '.join('<a href={v}>{k}</a>'.format(k=k, v=v) for k, v in record.html_fpath.items()))

    main_dict = dict()
    if full_report.sample_reports:
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