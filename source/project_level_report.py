import os
from os.path import join, relpath, dirname, basename, pardir, normpath, realpath
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


BASECALLS_NAME   = 'BaseCalls'
FASTQC_NAME      = BCBioStructure.fastqc_repr
PRE_FASTQC_NAME  = 'Raw ' + FASTQC_NAME
SEQQC_NAME       = 'Seq QC'
PRE_SEQQC_NAME   = 'Pre ' + SEQQC_NAME
VARQC_NAME       = BCBioStructure.varqc_repr
VARQC_AFTER_NAME = BCBioStructure.varqc_after_repr
MUTATIONS_NAME   = 'Mutations'


metric_storage = MetricStorage(
    general_section=ReportSection(metrics=[
        Metric(BASECALLS_NAME),
        Metric(PRE_FASTQC_NAME),
        Metric(FASTQC_NAME),
        Metric(PRE_SEQQC_NAME),
        Metric(SEQQC_NAME),
        Metric(VARQC_NAME),
        Metric(VARQC_AFTER_NAME),
        Metric(MUTATIONS_NAME)
    ]),
    sections=[ReportSection(metrics=[
        Metric(PRE_FASTQC_NAME),
        Metric(FASTQC_NAME),
        # Metric('BAM'),
        Metric(PRE_SEQQC_NAME),
        Metric(SEQQC_NAME),
        # Metric(MUTATIONS_NAME),
        Metric(VARQC_NAME),
        Metric(VARQC_AFTER_NAME),
    ])])


def make_project_level_report(cnf, dataset_structure=None, bcbio_structure=None):
    step_greetings('Making the preproc project-level report')

    if dataset_structure is None and bcbio_structure:
        analysis_dirpath = normpath(join(bcbio_structure.bcbio_project_dirpath, pardir))
        dataset_dirpath = realpath(join(analysis_dirpath, 'dataset'))
        dataset_structure = DatasetStructure.create(dataset_dirpath, bcbio_structure.project_name)

    general_records = _add_summary_reports(metric_storage.general_section, bcbio_structure, dataset_structure)
    sample_reports_records = _add_per_sample_reports(metric_storage.sections[0], bcbio_structure, dataset_structure)

    sample_reports = []
    for sample in (bcbio_structure or dataset_structure).samples:
        sample_reports.append(SampleReport(sample,
            records=sample_reports_records[sample.name],
            html_fpath=None,
            metric_storage=metric_storage))

    full_report = FullReport(cnf.project_name, sample_reports, metric_storage=metric_storage, general_records=general_records)
    project_report_html_fpath = dataset_structure.project_report_html_fpath
    project_name = dataset_structure.project_name
    if bcbio_structure:
        project_report_html_fpath = bcbio_structure.project_report_html_fpath
        project_name = bcbio_structure.project_name

    _save_static_html(cnf.work_dir, full_report, project_report_html_fpath, project_name)

    info()
    info('*' * 70)
    info('Project-level report saved in: ')
    info('  ' + project_report_html_fpath)
    return project_report_html_fpath


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
                html_fpath=_relpath_all(val, bcbio_structure.date_dirpath))
            general_section.add_metric(metric)
            records.append(rec)

    return records


def _make_path_record(html_fpath_value, metric, base_dirpath):
    if isinstance(html_fpath_value, basestring):
        html_fpath = relpath(html_fpath_value, base_dirpath) if verify_file(html_fpath_value) else None
        return Record(metric=metric, value=metric.name, html_fpath=html_fpath)
    elif isinstance(html_fpath_value, dict):
        html_fpath_value = OrderedDict([(k, relpath(html_fpath, base_dirpath)) for k, html_fpath in html_fpath_value.items() if verify_file(html_fpath)])
        return Record(metric=metric, value=metric.name, html_fpath=html_fpath_value)
    else:
        pass


def _add_summary_reports(general_section, bcbio_structure=None, dataset_structure=None):
    """ We want links to be relative, so we make paths relative to the project-level-report parent directory.
        - If the bcbio_structure is set, project-level report is located at bcbio_structure.date_dirpath
        - If dataset_dirpath is set, project-level report is located right at dataset_dirpath
    """
    if bcbio_structure:
        base_dirpath = bcbio_structure.date_dirpath
    else:
        base_dirpath = dirname(dataset_structure.project_report_html_fpath)

    recs = []

    if dataset_structure:
        if dataset_structure.basecall_stat_html_reports:
            val = OrderedDict([(basename(fpath), fpath) for fpath in dataset_structure.basecall_stat_html_reports])
            recs.append(_make_path_record(val, general_section.get_metric(BASECALLS_NAME), base_dirpath))
        recs.append(_make_path_record(dataset_structure.comb_fastqc_fpath,              general_section.get_metric(PRE_FASTQC_NAME), base_dirpath))
        recs.append(_make_path_record(dataset_structure.downsample_targqc_report_fpath, general_section.get_metric(PRE_SEQQC_NAME),  base_dirpath))

    if bcbio_structure:
        varqc_d = bcbio_structure.varqc_report_fpath_by_caller
        varqc_d['all'] = bcbio_structure.varqc_report_fpath

        varqc_after_d = bcbio_structure.varqc_after_report_fpath_by_caller
        varqc_after_d['all'] = bcbio_structure.varqc_after_report_fpath

        recs.append(_make_path_record(bcbio_structure.fastqc_summary_fpath, general_section.get_metric(FASTQC_NAME), base_dirpath))
        recs.append(_make_path_record(bcbio_structure.targqc_summary_fpath, general_section.get_metric(SEQQC_NAME),  base_dirpath))
        recs.append(_make_path_record(varqc_d,       general_section.get_metric(VARQC_NAME),       base_dirpath))
        recs.append(_make_path_record(varqc_after_d, general_section.get_metric(VARQC_AFTER_NAME), base_dirpath))

    return recs


def _add_per_sample_reports(individual_reports_section, bcbio_structure=None, dataset_structure=None):
    base_dirpath = None
    if dataset_structure:
        base_dirpath = dirname(dataset_structure.project_report_html_fpath)

    if bcbio_structure:
        base_dirpath = dirname(bcbio_structure.project_level_report_fpath)

    sample_reports_records = defaultdict(list)

    if dataset_structure:
        for s in dataset_structure.samples:
            sample_reports_records[s.name].extend([
                _make_path_record(
                    OrderedDict([('left', s.find_fastqc_html(s.l_fastqc_base_name)), ('right', s.find_fastqc_html(s.r_fastqc_base_name))]),
                    individual_reports_section.get_metric(PRE_FASTQC_NAME), base_dirpath),
                _make_path_record(
                    OrderedDict([('targqc', s.targetcov_html_fpath), ('ngscat', s.ngscat_html_fpath), ('qualimap', s.qualimap_html_fpath)]),
                    individual_reports_section.get_metric(PRE_SEQQC_NAME), base_dirpath)
            ])

    if bcbio_structure:
        for s in bcbio_structure.samples:
            targqc_d = OrderedDict([('targqc', s.targetcov_html_fpath), ('ngscat', s.ngscat_html_fpath), ('qualimap', s.qualimap_html_fpath)])
            varqc_d = OrderedDict([(k, s.get_varqc_fpath_by_callername(k)) for k in bcbio_structure.variant_callers.keys()])
            varqc_after_d = OrderedDict([(k, s.get_varqc_after_fpath_by_callername(k)) for k in bcbio_structure.variant_callers.keys()])

            sample_reports_records[s.name].extend([
                _make_path_record(s.fastqc_html_fpath, individual_reports_section.get_metric(FASTQC_NAME),      base_dirpath),
                _make_path_record(targqc_d,            individual_reports_section.get_metric(SEQQC_NAME),       base_dirpath),
                _make_path_record(varqc_d,             individual_reports_section.get_metric(VARQC_NAME),       base_dirpath),
                _make_path_record(varqc_after_d,       individual_reports_section.get_metric(VARQC_AFTER_NAME), base_dirpath)
            ])

    # for (repr_name, links_by_sample) in to_add:
    #     cur_metric = Metric(repr_name)
    #     # individual_reports_section.add_metric(cur_metric)
    #
    #     samples = []
    #     if dataset_structure:
    #         samples.extend(dataset_structure.samples)
    #     if bcbio_structure:
    #         samples.extend(bcbio_structure.samples)
    #
    #     for sample in samples:
    #         if links_by_sample and links_by_sample.get(sample.name):
    #             sample_reports_records[sample.name].append(
    #                 Record(
    #                     metric=cur_metric,
    #                     value=cur_metric.name,
    #                     html_fpath=_relpath_all(links_by_sample[sample.name], base_dirpath)))
    #         else:
    #             sample_reports_records[sample.name].append(
    #                 Record(metric=cur_metric, value=None, html_fpath=None))

    # fastqc_by_sample = dict()
    # targqc_by_sample = dict()
    # varqc_by_sample = dict()
    # varqc_after_by_sample = dict()
    #
    # if dataset_structure:
    #
    #     fastqc_by_sample.append(dict(PRE_FASTQC_NAME=pre_fastqc_htmls_by_sample))
    #     pre_fastqc_htmls_by_sample = dict([(s.name, verify_file(s.fastqc_html_fpath)) for s in dataset_structure.samples])
    #     targqc_htmls_by_sample     = _add_targqc_reports(dataset_structure.samples)
    #
    #     to_add.extend([
    #         (PRE_FASTQC_NAME,        pre_fastqc_htmls_by_sample),
    #         (PRE_SEQQC_NAME,         targqc_htmls_by_sample),
    #     ])
    #
    # if bcbio_structure:
    #     varqc_htmls_by_sample       = _add_varqc_reports(bcbio_structure, BCBioStructure.varqc_name, BCBioStructure.varqc_dir)
    #     varqc_after_htmls_by_sample = _add_varqc_reports(bcbio_structure, BCBioStructure.varqc_after_name, BCBioStructure.varqc_after_dir)
    #     targqc_htmls_by_sample      = _add_targqc_reports(bcbio_structure.samples)
    #     fastqc_htmls_by_sample      = dict([(s.name, verify_file(s.fastqc_html_fpath)) for s in bcbio_structure.samples])
    #     to_add.extend([
    #         (FASTQC_NAME,      fastqc_htmls_by_sample),
    #         # ('BAM',                            bams_by_samples),
    #         (SEQQC_NAME,       targqc_htmls_by_sample),
    #         (VARQC_NAME,       varqc_htmls_by_sample),
    #         (VARQC_AFTER_NAME, varqc_after_htmls_by_sample)
    #     ])

    return sample_reports_records


# def _add_varqc_reports(bcbio_structure, name, dir_name):
#     callers = bcbio_structure.variant_callers.values()
#     if len(callers) == 0:
#         varqc_htmls_by_sample = None
#     elif len(callers) == 1:
#         varqc_htmls_by_sample = callers[0].find_fpaths_by_sample(dir_name, name, 'html')
#     else:
#         varqc_htmls_by_sample = OrderedDict()
#         for sample in bcbio_structure.samples:
#             varqc_htmls_by_sample[sample.name] = OrderedDict()
#         for caller in callers:
#             for sample, fpath in caller.find_fpaths_by_sample(dir_name, name, 'html').items():
#                 varqc_htmls_by_sample[sample][caller.name] = fpath
#
#     return varqc_htmls_by_sample
#
#
# def _add_targqc_reports(samples):
#     targqc_htmls_by_sample = OrderedDict()
#
#     for sample in samples:
#         targqc_htmls_by_sample[sample.name] = OrderedDict()
#         for report_name, report_html_fpath in [
#             ('targetcov', sample.targetcov_html_fpath),
#             ('ngscat',    sample.ngscat_html_fpath),
#             ('qualimap',  sample.qualimap_html_fpath)]:
#             if verify_file(report_html_fpath):
#                 targqc_htmls_by_sample[sample.name][report_name] = verify_file(report_html_fpath)
#
#     return targqc_htmls_by_sample


def _relpath_all(value, base_dirpath):
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

    def __process_record(rec, short=False):
        d = rec.__dict__.copy()

        if isinstance(rec.html_fpath, basestring):
            d['contents'] = '<a href="' + rec.html_fpath + '">' + rec.value + '</a>'

        elif isinstance(rec.html_fpath, dict):
            d['contents'] = ', '.join('<a href="{v}">{k}</a>'.format(k=k, v=v) for k, v in rec.html_fpath.items()) if rec.html_fpath else '-'
            if not short:
                d['contents'] = rec.metric.name + ': ' + d['contents']

        else:
            d['contents'] = '-'

        d['metric'] = rec.metric.__dict__
        return d

    def _get_summary_report_name(rec):
        return rec.value.lower().replace(' ', '_')

    # common records (summary reports)
    common_dict = dict()
    common_dict["project_name"] = project_name
    common_records = full_report.get_common_records()
    for rec in common_records:
        if rec.value:
            common_dict[_get_summary_report_name(rec)] = __process_record(rec)  # rec_d

    main_dict = dict()
    if full_report.sample_reports:
        # individual records
        main_dict = dict()
        main_dict["sample_reports"] = []

        metrics = metric_storage.get_metrics(skip_general_section=True)
        metrics_with_values_set = set()
        for sample_report in full_report.sample_reports:
            for m in metric_storage.get_metrics(skip_general_section=True):
                r = next((r for r in sample_report.records if r.metric.name == m.name), None)
                if r:
                    metrics_with_values_set.add(m)

        metrics = [m for m in metrics if m in metrics_with_values_set]
        main_dict['metric_names'] = [m.name for m in metrics]

        for sample_report in full_report.sample_reports:
            ready_records = []
            for m in metrics:
                r = next((r for r in sample_report.records if r.metric.name == m.name), None)
                if r:
                    ready_records.append(__process_record(r, short=True))
                else:
                    ready_records.append(__process_record(Record(metric=m, value=None), short=True))
            assert len(ready_records) == len(main_dict["metric_names"])

            sample_report_dict = dict()
            sample_report_dict["records"] = ready_records
            sample_report_dict["sample_name"] = sample_report.get_display_name()
            main_dict["sample_reports"].append(sample_report_dict)

    return _write_static_html_report(work_dir, {"common": common_dict, "main": main_dict}, html_fpath)