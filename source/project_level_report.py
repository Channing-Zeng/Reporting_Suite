import getpass
from os.path import join, relpath
from collections import OrderedDict
import paramiko
import base64

from source.bcbio_structure import BCBioStructure
from source.logger import info, step_greetings, send_email, warn, err
from source.file_utils import verify_file
from source.reporting import Metric, Record, MetricStorage, ReportSection, SampleReport, FullReport
from source.html_reporting.html_saver import write_static_html_report


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
    # final_summary_report_fpath = full_report.save_html(
    #     bcbio_structure.date_dirpath, bcbio_structure.project_name,
    #     'Project-level report for ' + bcbio_structure.project_name)
    final_summary_report_fpath = _save_static_html(full_report, bcbio_structure.date_dirpath,
        report_base_name=bcbio_structure.project_name,
        project_name=bcbio_structure.project_name)

    copy_to_ngs_website(bcbio_structure.final_dirpath, final_summary_report_fpath, bcbio_structure.project_name)

    info()
    info('*' * 70)
    info('Project-level report saved in: ')
    info('  ' + final_summary_report_fpath)
    send_email('Report for ' + bcbio_structure.project_name + ':\n  ' + final_summary_report_fpath)


def copy_to_ngs_website(final_dirpath, html_report_fpath, project_name):
    server_url = 'ngs.usbod.astrazeneca.net'
    server_path = '/opt/lampp/htdocs/reports'
    username = 'klpf990'
    password = '123werasd'
    project_list_fpath = '/ngs/oncology/NGS.Project.csv'

    try:
        client = paramiko.SSHClient()
        client.connect(server_url, username=username, password=password)
    except:
        warn('Cannot connect to ' + server_url)
        pass
    else:
        client.exec_command('cd ' + server_path)
        client.exec_command('ln -s ' + final_dirpath + ' ' + project_name)
        client.close()

    if verify_file(project_list_fpath, 'Project list'):
        pids = set()
        with open(project_list_fpath) as f:
            fields = f.readline().strip().split(',')  # 'Updated By,PID,Name,JIRA URL,HTML report path,Why_IfNoReport,Data Hub,Analyses directory UK,Analyses directory US,Type,Division,Department,Sample Number,Reporter,Assignee,Description,IGV,Notes'
            index_of_pid = fields.index('PID')
            for l in f:
                values = l.strip().split(',')
                pids.add(values[index_of_pid])

        if project_name not in pids:
            with open(project_list_fpath, 'a') as f:
                values = {'PID': project_name, 'HTML report path': html_report_fpath, 'Updated By': getpass.getuser()}
                f.write('\n' + ','.join(values.get(f, '') for f in fields))


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
    fastqc_htmls_by_sample      = dict([(s.name, verify_file(s.fastqc_html_fpath)) for s in bcbio_structure.samples])

    sample_reports_records = dict()
    for sample in bcbio_structure.samples:
        sample_reports_records[sample.name] = list(general_records)

    for (repr_name, htmls_by_sample) in [
            (bcbio_structure.fastqc_repr,      fastqc_htmls_by_sample),
            (bcbio_structure.targqc_repr,      targqc_htmls_by_sample),
            (bcbio_structure.varqc_repr,       varqc_htmls_by_sample),
            (bcbio_structure.varqc_after_repr, varqc_after_htmls_by_sample)]:
        cur_metric = Metric(repr_name)
        individual_reports_section.add_metric(cur_metric)
        for sample in bcbio_structure.samples:
            if htmls_by_sample and htmls_by_sample.get(sample.name):
                sample_reports_records[sample.name].append(
                    Record(
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
    targqc_htmls_by_sample = OrderedDict()

    for sample in bcbio_structure.samples:
        targqc_htmls_by_sample[sample.name] = OrderedDict()
        for report_name, report_html_fpath in [('targetcov', sample.targetcov_html_fpath),
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


def _save_static_html(full_report, output_dirpath, report_base_name, project_name):
    # metric name in FullReport --> metric name in Static HTML
    metric_names = OrderedDict([
        ('FastQC', 'FastQC'),
        ('Target QC', 'SeqQC'),
        ('Var QC', 'VarQC'),
        ('Var QC after filtering', 'VarQC after filtering')])

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
    sample_report = full_report.sample_reports[0]
    for record in sample_report.records:
        if record.metric.common:
            common_dict[_get_summary_report_name(record)] = record.__dict__

    # individual records
    main_dict = dict()
    main_dict["sample_reports"] = []
    main_dict["metric_names"] = metric_names.values()
    for sample_report in full_report.sample_reports:
        new_records = [_process_record(record.__dict__) for record in sample_report.records
                       if record.metric.name in metric_names.keys()]
        sample_report_dict = dict()
        sample_report_dict["records"] = new_records
        sample_report_dict["sample_name"] = sample_report.display_name
        main_dict["sample_reports"].append(sample_report_dict)

    return write_static_html_report({"common": common_dict, "main": main_dict},
                                    output_dirpath, report_base_name)
