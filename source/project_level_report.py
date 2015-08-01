import os
from os.path import join, relpath, dirname
from collections import OrderedDict, namedtuple
import getpass
from traceback import format_exc
from collections import defaultdict

import source
from source.preproc.dataset_structure import DatasetStructure
from source.bcbio_structure import BCBioStructure
from source.jira_utils import JiraCase
from source.logger import info, step_greetings, send_email, warn, err
from source.file_utils import verify_file, file_transaction, adjust_path, safe_mkdir, add_suffix
from source.reporting import Metric, Record, MetricStorage, ReportSection, SampleReport, FullReport
from source.html_reporting.html_saver import _write_static_html_report


FASTQC_NAME      = BCBioStructure.fastqc_repr
PRE_FASTQC_NAME  = 'Raw ' + FASTQC_NAME
SEQQC_NAME       = 'Seq QC'
PRE_SEQQC_NAME   = 'Downsampled ' + SEQQC_NAME
VARQC_NAME       = BCBioStructure.varqc_repr
VARQC_AFTER_NAME = BCBioStructure.varqc_after_repr
MUTATIONS_NAME   = 'Mutations'


metric_storage = MetricStorage(
    general_section=ReportSection(metrics=[
        Metric(PRE_FASTQC_NAME),
        Metric(FASTQC_NAME),
        Metric(PRE_SEQQC_NAME),
        Metric(SEQQC_NAME),
        Metric(VARQC_NAME),
        Metric(VARQC_AFTER_NAME)
    ]),
    sections=[ReportSection(metrics=[
        Metric(PRE_FASTQC_NAME),
        Metric(FASTQC_NAME),
        # Metric('BAM'),
        Metric(PRE_SEQQC_NAME),
        Metric(SEQQC_NAME),
        Metric(MUTATIONS_NAME),
        Metric(VARQC_NAME),
        Metric(VARQC_AFTER_NAME),
    ])])


def make_preproc_project_level_report(cnf, dataset_structure):
    step_greetings('Preproc project-level report')

    general_records = _add_summary_reports(metric_storage.general_section, dataset_structure=dataset_structure)
    sample_reports_records = _add_per_sample_reports(metric_storage.sections[0], dataset_structure=dataset_structure)

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

    general_records = _add_summary_reports(metric_storage.general_section, bcbio_structure=bcbio_structure)
    general_records.extend(_add_variants(metric_storage.general_section, bcbio_structure))

    sample_reports_records = _add_per_sample_reports(metric_storage.sections[0], bcbio_structure=bcbio_structure)

    sample_reports = []
    for sample in bcbio_structure.samples:
        sample_reports.append(SampleReport(
            sample,
            records=sample_reports_records[sample.name],
            html_fpath=None,
            metric_storage=metric_storage))

    full_report = FullReport(cnf.project_name, sample_reports, metric_storage=metric_storage, general_records=general_records)
    _save_static_html(cnf.work_dir, full_report, bcbio_structure.project_level_report_fpath, bcbio_structure.project_name)

    info()
    info('*' * 70)
    info('Project-level report saved in: ')
    info('  ' + bcbio_structure.project_level_report_fpath)
    return bcbio_structure.project_level_report_fpath


def _add_variants(general_section, bcbio_structure):
    records = []

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
            records.append(rec)

    return records


def _add_summary_reports(general_section, bcbio_structure=None, dataset_structure=None):
    """ We want links to be relative, so we make paths relative to the project-level-report parent directory.
        - If the bcbio_structure is set, project-level report is located at bcbio_structure.date_dirpath
        - If dataset_dirpath is set, project-level report is located right at dataset_dirpath
    """
    records = []

    if bcbio_structure:
        base_dirpath = bcbio_structure.date_dirpath
    else:
        base_dirpath = dirname(dataset_structure.project_report_html_fpath)

    if dataset_structure:
        for metric_name, report_fpath in [
            (PRE_FASTQC_NAME, dataset_structure.comb_fastqc_fpath),
            (PRE_SEQQC_NAME,  dataset_structure.downsample_targqc_report_fpath)]:

            if verify_file(report_fpath):
                metric = general_section.get_metric(metric_name)
                records.append(Record(metric=metric, value=metric.name, html_fpath=relpath(report_fpath, base_dirpath)))

    if bcbio_structure:
        for (report_file_name, metric_name, summary_dir) in [
            (BCBioStructure.fastqc_name,      FASTQC_NAME,      BCBioStructure.fastqc_summary_dir),
            (BCBioStructure.targqc_name,      SEQQC_NAME,       BCBioStructure.targqc_summary_dir),
            (BCBioStructure.varqc_name,       VARQC_NAME,       BCBioStructure.varqc_summary_dir),
            (BCBioStructure.varqc_after_name, VARQC_AFTER_NAME, BCBioStructure.varqc_after_summary_dir)]:

            report_fpath = join(bcbio_structure.date_dirpath, summary_dir, report_file_name + '.html')
            if verify_file(report_fpath):
                metric = general_section.get_metric(metric_name)
                if metric:
                    records.append(Record(metric=metric, value=metric.name, html_fpath=relpath(report_fpath, base_dirpath)))
            else:
                pass

    return records


def _add_per_sample_reports(individual_reports_section, bcbio_structure=None, dataset_structure=None):
    to_add = []
    base_dirpath = None

    if dataset_structure:
        pre_fastqc_htmls_by_sample = dict([(s.name, verify_file(s.fastqc_html_fpath)) for s in dataset_structure.samples])
        targqc_htmls_by_sample     = _add_targqc_reports(dataset_structure.samples)
        to_add.extend([
            (PRE_FASTQC_NAME,        pre_fastqc_htmls_by_sample),
            (PRE_SEQQC_NAME,         targqc_htmls_by_sample),
        ])
        base_dirpath = dirname(dataset_structure.project_report_html_fpath)

    if bcbio_structure:
        varqc_htmls_by_sample       = _add_varqc_reports(bcbio_structure, BCBioStructure.varqc_name, BCBioStructure.varqc_dir)
        varqc_after_htmls_by_sample = _add_varqc_reports(bcbio_structure, BCBioStructure.varqc_after_name, BCBioStructure.varqc_after_dir)
        targqc_htmls_by_sample      = _add_targqc_reports(bcbio_structure.samples)
        fastqc_htmls_by_sample      = dict([(s.name, verify_file(s.fastqc_html_fpath)) for s in bcbio_structure.samples])
        to_add.extend([
            (FASTQC_NAME,      fastqc_htmls_by_sample),
            # ('BAM',                            bams_by_samples),
            (SEQQC_NAME,       targqc_htmls_by_sample),
            (VARQC_NAME,       varqc_htmls_by_sample),
            (VARQC_AFTER_NAME, varqc_after_htmls_by_sample)
        ])
        base_dirpath = dirname(bcbio_structure.project_level_report_fpath)

    sample_reports_records = defaultdict(list)

    for (repr_name, links_by_sample) in to_add:
        cur_metric = Metric(repr_name)
        # individual_reports_section.add_metric(cur_metric)

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
    # metric_names = OrderedDict([
    #     (DatasetStructure.pre_fastqc_repr, DatasetStructure.pre_fastqc_repr),
    #     (BCBioStructure.fastqc_repr, 'FastQC'),
    #     # ('BAM', 'BAM'),
    #     (BCBioStructure.targqc_repr, 'SeqQC'),
    #     # ('Downsampled ' + BCBioStructure.targqc_repr, 'Downsampled SeqQC'),
    #     # ('Mutations', 'Mutations'),
    #     # ('Mutations for separate samples', 'Mutations for separate samples'),
    #     # ('Mutations for paired samples', 'Mutations for paired samples'),
    #     (BCBioStructure.varqc_repr, 'VarQC'),
    #     (BCBioStructure.varqc_after_repr, 'VarQC after filtering')])

    def __process_record(rec):
        new_html_fpath_value = []
        if isinstance(rec.html_fpath, basestring):
            new_html_fpath_value = [dict(html_fpath_name=rec.value, html_fpath_value=rec.html_fpath)]
        elif isinstance(rec.html_fpath, dict):
            for k, v in rec.html_fpath.items():
                new_html_fpath_value.append(dict(html_fpath_name=k, html_fpath_value=v))
        rec.html_fpath = new_html_fpath_value
        return rec.__dict__

    def _get_summary_report_name(rec):
        return rec.value.lower().replace(' ', '_')

    # common records (summary reports)
    common_dict = dict()
    common_dict["project_name"] = project_name
    common_records = full_report.get_common_records()
    for rec in common_records:
        rec.metric = rec.metric.__dict__
        rec_d = rec.__dict__
        common_dict[_get_summary_report_name(rec)] = rec_d
        if 'Mutations' in rec.metric['name']:
            common_dict[_get_summary_report_name(rec)]['contents'] = (
                rec.metric.name + ': ' + ', '.join('<a href={v}>{k}</a>'.format(k=k, v=v) for k, v in rec.html_fpath.items()))

    main_dict = dict()
    if full_report.sample_reports:
        # individual records
        main_dict = dict()
        main_dict["sample_reports"] = []
        main_dict["metric_names"] = [m.name for m in full_report.metric_storage.get_metrics(skip_general_section=True)]

        for sample_report in full_report.sample_reports:
            ready_records = []
            for m in metric_storage.get_metrics(skip_general_section=True):
                r = next((r for r in sample_report.records if r.metric.name == m.name), None)
                if r:
                    ready_records.append(__process_record(r))
                else:
                    ready_records.append(__process_record(Record(metric=m, value=None)))
            assert len(ready_records) == len(main_dict["metric_names"])

            sample_report_dict = dict()
            sample_report_dict["records"] = ready_records
            sample_report_dict["sample_name"] = sample_report.get_display_name()
            main_dict["sample_reports"].append(sample_report_dict)

    return _write_static_html_report(work_dir, {"common": common_dict, "main": main_dict}, html_fpath)