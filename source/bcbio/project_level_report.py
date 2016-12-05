import glob
import json
import os
import time
from genericpath import isdir, exists
from inspect import getsourcefile
from os import listdir
from os.path import join, relpath, dirname, basename, abspath, getmtime, isfile, pardir, islink
from collections import OrderedDict
from collections import defaultdict
import shutil

import variant_filtering
from ngs_reporting.oncoprints import create_oncoprints_link
from ngs_utils.file_utils import safe_symlink_to

from source.bcbio.bcbio_structure import BCBioStructure
from source.calling_process import call
from source.logger import info, step_greetings, warn, timestamp, critical, err
from source.file_utils import verify_file, add_suffix, verify_dir, file_transaction, safe_mkdir
from source.reporting.reporting import Metric, Record, MetricStorage, ReportSection, SampleReport, FullReport, \
    write_static_html_report
from source.tools_from_cnf import get_system_path, get_script_cmdline
from source.utils import get_ext_tools_dirname, get_version
from tools.prepare_data_for_exac import get_exac_us_url


FASTQC_NAME           = BCBioStructure.fastqc_repr
PRE_FASTQC_NAME       = 'Raw ' + FASTQC_NAME
EXAC_NAME             = 'ExAC browser'
CNV_NAME              = 'CNV'
MUTATIONS_NAME        = 'VarDict'
MUTATIONS_SINGLE_NAME = 'VarDict, single samples'
MUTATIONS_PAIRED_NAME = 'VarDict, paired samples'
GENDER                = 'Sex'
CLINICAL_NAME         = 'Oncology NGS report'
PHENOTYPE             = 'Phenotype'
NORM_MATCH            = 'Normal Match'
ABNORMAL_NAME         = 'Flagged regions'

## RNAseq reports
GENE_COUNTS_NAME      = 'Gene counts'
EXON_COUNTS_NAME      = 'Exon counts'
GENE_TPM_NAME         = 'Gene TPM'
ISOFORM_TPM_NAME      = 'Isoform TPM'

mutation_names = [
    MUTATIONS_NAME,
    MUTATIONS_SINGLE_NAME,
    MUTATIONS_PAIRED_NAME,
    CNV_NAME
]
expression_names = [
    GENE_COUNTS_NAME,
    EXON_COUNTS_NAME,
    GENE_TPM_NAME,
    ISOFORM_TPM_NAME,
]

metric_storage = MetricStorage(
    general_section=ReportSection(metrics=[
        Metric(PRE_FASTQC_NAME),
        Metric(FASTQC_NAME),
        Metric(EXAC_NAME),
        Metric(MUTATIONS_NAME),
        Metric(MUTATIONS_SINGLE_NAME),
        Metric(MUTATIONS_PAIRED_NAME),
        Metric(ABNORMAL_NAME),
        Metric(GENE_COUNTS_NAME),
        Metric(EXON_COUNTS_NAME),
        Metric(GENE_TPM_NAME),
        Metric(ISOFORM_TPM_NAME),
    ]),
    sections=[ReportSection(metrics=[
        Metric(PRE_FASTQC_NAME),
        Metric(FASTQC_NAME),
        # Metric('BAM'),
        # Metric(MUTATIONS_NAME),
        Metric(GENDER, description='If not defined, means that the target does not contain key male Y genes that we could check'),
        Metric(CLINICAL_NAME),
        Metric(CNV_NAME),
        Metric(GENE_COUNTS_NAME),
        Metric(PHENOTYPE),
        Metric(NORM_MATCH),
    ])])


def make_report_metadata(cnf, bcbio_structure, oncoprints_link=None):
    step_greetings('Making the %s project-level report' % ('preproc' if bcbio_structure is None else 'postproc'))

    # if dataset_structure is None and bcbio_structure:
    #     analysis_dirpath = normpath(join(bcbio_structure.bcbio_project_dirpath, pardir))
    #     dataset_dirpath = realpath(join(analysis_dirpath, 'dataset'))
    #     dataset_structure = DatasetStructure.create(dataset_dirpath, bcbio_structure.project_name)

    general_records = _add_summary_reports(cnf, metric_storage.general_section, bcbio_structure)
    # sample_reports_records = _add_per_sample_reports(cnf, metric_storage.sections[0], bcbio_structure)

    # sample_reports = []
    # if dataset_project:
    #     samples = dataset_project.sample_by_name.values()
    # if bcbio_structure:
    samples = bcbio_structure.samples
    # for sample in samples:
    #     sample_reports.append(SampleReport(sample,
    #         records=sample_reports_records[sample.name],
    #         html_fpath=None,
    #         metric_storage=metric_storage))

    full_report = FullReport(cnf.project_name, [], metric_storage=metric_storage, general_records=general_records)

    project_report_html_fpath = bcbio_structure.multiqc_fpath
    project_name = bcbio_structure.project_name

    additional_data = dict()
    normal_samples = [s for s in bcbio_structure.samples if s.phenotype == 'normal']
    if normal_samples:
        sample_match_on_hover_js = '<script type="text/javascript">\n'
        for s in bcbio_structure.samples:
            if s.phenotype != 'normal' and s.normal_match:
                sample_match_on_hover_js += ('' +
                    '\tdocument.getElementById("' + s.name + '_match").onmouseover = function() { document.getElementById("' + s.normal_match.name + '").style.backgroundColor = "#EEE"; };\n' +
                    '\tdocument.getElementById("' + s.name + '_match").onmouseleave = function() { document.getElementById("' + s.normal_match.name + '").style.backgroundColor = "white"; };\n'
                )
        sample_match_on_hover_js += '</script>\n'
        additional_data['sample_match_on_hover_js'] = sample_match_on_hover_js

    metadata = _report_to_multiqc_metadata(cnf, full_report,
        project_report_html_fpath, project_name, bcbio_structure,
        additional_data=additional_data, oncoprints_link=oncoprints_link)

    metadata_fpath = join(bcbio_structure.work_dir, 'az_multiqc_metadata.yaml')
    import yaml
    with open(metadata_fpath, 'w') as outfile:
        yaml.dump(metadata, outfile, default_flow_style=False)
    import json
    with open(metadata_fpath.replace('.yaml', '.json'), 'w') as outfile:
        json.dump(metadata, outfile)

    return metadata_fpath.replace('.yaml', '.json')


def _mutations_records(general_section, bcbio_structure, base_dirpath):
    records = []

    caller = bcbio_structure.variant_callers.get('vardict') or \
             bcbio_structure.variant_callers.get('vardict-java')

    _base_mut_fname = variant_filtering.mut_fname_template.format(caller_name=caller.name)
    _base_mut_fpath = join(bcbio_structure.date_dirpath, _base_mut_fname)
    mut_fpath = add_suffix(_base_mut_fpath, variant_filtering.mut_pass_suffix)
    single_mut_fpath = add_suffix(add_suffix(_base_mut_fpath, variant_filtering.mut_single_suffix), variant_filtering.mut_pass_suffix)
    paired_mut_fpath = add_suffix(add_suffix(_base_mut_fpath, variant_filtering.mut_paired_suffix), variant_filtering.mut_pass_suffix)
    mut_fpath = verify_file(mut_fpath, silent=True)
    single_mut_fpath = verify_file(single_mut_fpath, silent=True)
    paired_mut_fpath = verify_file(paired_mut_fpath, silent=True)

    for fpath, metric_name in (
         (mut_fpath, MUTATIONS_NAME),
         (single_mut_fpath, MUTATIONS_SINGLE_NAME),
         (paired_mut_fpath, MUTATIONS_PAIRED_NAME)):
        if fpath:
            metric = Metric(metric_name, common=True)
            rec = Record(
                metric=metric,
                value=basename(fpath),
                url=relpath(fpath, base_dirpath))
            general_section.add_metric(metric)
            records.append(rec)

    if isfile(bcbio_structure.seq2c_fpath):
        metric = Metric(CNV_NAME, common=True)
        fpath = bcbio_structure.seq2c_fpath
        rec = Record(
            metric=metric,
            value=basename(fpath),
            url=relpath(fpath, base_dirpath))
        general_section.add_metric(metric)
        records.append(rec)

    return records


def _make_url_record(html_fpath_value, metric, base_dirpath):
    # info('Adding paths to the report: ' + str(html_fpath_value))
    # if isinstance(html_fpath_value, dict):
    #     url = OrderedDict([(k, relpath(html_fpath, base_dirpath)) for k, html_fpath in html_fpath_value.items() if verify_file(html_fpath)])
    #     return Record(metric=metric, value=metric.name, url=url)
    # else:
    url = relpath(html_fpath_value, base_dirpath) if verify_file(html_fpath_value) else None
    return Record(metric=metric, value=metric.name, url=url)


def _add_summary_reports(cnf, general_section, bcbio_structure):
    """ We want links to be relative, so we make paths relative to the project-level-report parent directory.
        - If the bcbio_structure is set, project-level report is located at bcbio_structure.date_dirpath
        - If dataset_dirpath is set, project-level report is located right at dataset_dirpath
    """
    base_dirpath = dirname(bcbio_structure.multiqc_fpath)

    recs = []
    if bcbio_structure.is_rnaseq:
        recs = add_rna_summary_records(cnf, recs, general_section, bcbio_structure, base_dirpath)
    else:
        recs = add_dna_summary_records(cnf, recs, general_section, bcbio_structure, base_dirpath)
        # recs.append(_make_url_record(bcbio_structure.targqc_summary_fpath, general_section.find_metric(SEQQC_NAME), base_dirpath))
    # if verify_dir(bcbio_structure.flagged_regions_dirpath, is_critical=False):
    #     url_val = OrderedDict(
    #             [(region_type, join(bcbio_structure.flagged_regions_dirpath, 'flagged_' + region_type + '.html'))
    #                 for region_type in ['target', 'exons']])
    #     rec = _make_url_record(url_val, general_section.find_metric(ABNORMAL_NAME), base_dirpath)
    #     recs.append(rec)
    return recs


def add_rna_summary_records(cnf, recs, general_section, bcbio_structure, base_dirpath):
    recs.append(_make_url_record(
            join(bcbio_structure.expression_dirpath, 'html', 'counts.html'),
            general_section.find_metric(GENE_COUNTS_NAME), base_dirpath))
    recs.append(_make_url_record(
            join(bcbio_structure.expression_dirpath, 'html', 'dexseq.html'),
            general_section.find_metric(EXON_COUNTS_NAME), base_dirpath))
    recs.append(_make_url_record(
            join(bcbio_structure.expression_dirpath, 'html', 'gene.sf.tpm.html'),
            general_section.find_metric(GENE_TPM_NAME),    base_dirpath))
    recs.append(_make_url_record(
            join(bcbio_structure.expression_dirpath, 'html', 'isoform.sf.tpm.html'),
            general_section.find_metric(ISOFORM_TPM_NAME), base_dirpath))

    # rnaseq_html_fpath = join(bcbio_structure.date_dirpath, BCBioStructure.rnaseq_qc_report_name + '.html')
    # rnaseq_html_fpath = verify_file(rnaseq_html_fpath, is_critical=True)
    # recs.append(_make_url_record(rnaseq_html_fpath, general_section.find_metric(QC_REPORT_NAME), base_dirpath))

    return recs


def add_dna_summary_records(cnf, recs, general_section, bcbio_structure, base_dirpath):
    recs.extend(_mutations_records(general_section, bcbio_structure, base_dirpath))

    # m = general_section.find_metric(EXAC_NAME)
    # recs.append(Record(metric=m, value=m.name, url=get_exac_us_url(cnf.genome.name, bcbio_structure.project_name)))

    return recs


def make_multiqc_report(cnf, bcbio_structure, metadata_fpath=None):
    multiqc_bcbio_dirpath = join(bcbio_structure.date_dirpath, 'multiqc')
    new_multiqc_bcbio_dirpath = join(bcbio_structure.date_dirpath, 'log', 'multiqc_bcbio')
    if isdir(multiqc_bcbio_dirpath):
        if isdir(new_multiqc_bcbio_dirpath):
            try:
                shutil.rmtree(new_multiqc_bcbio_dirpath)
            except OSError:
                os.rename(new_multiqc_bcbio_dirpath, new_multiqc_bcbio_dirpath + '.' + timestamp().replace(':', '_').replace(' ', '_'))
        os.rename(multiqc_bcbio_dirpath, new_multiqc_bcbio_dirpath)

    multiqc_postproc_dirpath = safe_mkdir(dirname(bcbio_structure.multiqc_fpath))
    cmdl = 'multiqc -f -o ' + multiqc_postproc_dirpath + ' -t az ' + \
           '--filename ' + basename(bcbio_structure.multiqc_fpath) + ' --no-data-dir'
    if cnf.debug:
        cmdl += ' -v'
    if metadata_fpath:
        cmdl += ' --az-metadata ' + metadata_fpath

    input_list_fpath = join(bcbio_structure.date_dirpath, 'log', 'multiqc_list_files.txt')
    if not cnf.reuse_intermediate or not verify_file(input_list_fpath, silent=True):
        work_input_list_fpath = abspath(join(bcbio_structure.work_dir, pardir, 'list_files.txt'))
        if not verify_file(work_input_list_fpath, silent=True):
            critical('File list for MultiQC ' + work_input_list_fpath + ' not found. Please make sure you are using bcbio 1.0.0')
        with open(work_input_list_fpath) as inp, open(input_list_fpath, 'w') as out:
            qc_files_not_found = []
            for l in inp:
                fpath = l.strip()
                if '/work/' in fpath:
                    if fpath.endswith('target_info.yaml'):
                        correct_fpath = join(new_multiqc_bcbio_dirpath, basename(fpath))
                        if verify_file(fpath):
                            shutil.copy(fpath, correct_fpath)
                            out.write(correct_fpath + '\n')
                        else:
                            qc_files_not_found.append(fpath)
                            continue
                    else:
                        work_dirpath = fpath.split('/work/')[0] + '/work'
                        correct_fpath = fpath.replace(work_dirpath, bcbio_structure.final_dirpath)
                        for s in bcbio_structure.samples:
                            if '/qc/' + s.name + '/' in fpath:
                                correct_fpath = correct_fpath.replace('/qc/' + s.name + '/',
                                                                      '/' + s.name + '/qc/')
                                if '/qualimap/' in fpath:
                                    # QualiMap needs to be located in a directory named after the sample name:
                                    if not islink(join(s.dirpath, 'qc', 'qualimap', s.name)):
                                        os.symlink(join(s.dirpath, 'qc', 'qualimap'), join(s.dirpath, 'qc', 'qualimap', s.name))
                                elif '/qualimap_rnaseq/' in fpath:
                                    # QualiMap needs to be located in a directory named after the sample name:
                                    if not islink(join(s.dirpath, 'qc', 'qualimap_rnaseq', s.name)):
                                        os.symlink(join(s.dirpath, 'qc', 'qualimap_rnaseq'), join(s.dirpath, 'qc', 'qualimap_rnaseq', s.name))
                            else:
                                correct_fpath = correct_fpath.replace('/multiqc/report/metrics/' + s.name + '_bcbio.txt',
                                                                      '/' + s.name + '/qc/bcbio/' + s.name + '_bcbio.txt')
                        if not verify_file(correct_fpath):
                            qc_files_not_found.append(correct_fpath)
                            continue
                        out.write(correct_fpath + '\n')
                else:
                    if not verify_file(fpath):
                        qc_files_not_found.append(fpath)
                        continue
                    out.write(fpath + '\n')
            if qc_files_not_found:
                err('-')
                critical('Some QC files from list ' + work_input_list_fpath + ' were not found. Please check if work directory exists.')
            if bcbio_structure.is_rnaseq:
                pca_plot_fpath = create_rnaseq_pca_plot(cnf, bcbio_structure)
                info()
                if pca_plot_fpath and verify_file(pca_plot_fpath):
                    out.write(pca_plot_fpath + '\n')
            else:
                for s in bcbio_structure.samples:
                    snpeff_log_fpath = join(s.dirpath, BCBioStructure.varannotate_dir, s.name + '-vardict.snpEff_summary.csv')
                    if verify_file(snpeff_log_fpath, silent=True):
                        snpeff_log_multiqc_fpath = join(safe_mkdir(join(multiqc_postproc_dirpath, 'snpEff')), s.name)
                        shutil.copy(snpeff_log_fpath, snpeff_log_multiqc_fpath)
                        out.write(snpeff_log_multiqc_fpath + '\n')

    if verify_file(input_list_fpath, silent=True):
        cmdl += ' -l ' + input_list_fpath
    else:
        if not isfile(input_list_fpath):
            critical('Critical: MultiQC files list was not found in ' + input_list_fpath)
        # else:
            # targqc_dirpath = join(bcbio_structure.date_dirpath, BCBioStructure.targqc_dir)
            # if isdir(targqc_dirpath):
            #     info('Adding TargQC to MultiQC run')
            #     lines = open(input_list_fpath).readlines()
            #     lines.append(targqc_dirpath)
            #     with open(input_list_fpath, 'w') as f:
            #         f.writelines(lines)
            #     cmdl += ' -e qualimap'
            #     to_run = True
    if bcbio_structure.is_rnaseq:
        config_fname = 'multiqc_config_rna.yaml'
        cmdl += ' -e bcbio'
    else:
        config_fname = 'multiqc_config_dna.yaml'
        cmdl += ' -e qualimap -e bcbio'
    config_fpath = join(dirname(dirname(__file__)), config_fname)
    cmdl += ' -c ' + config_fpath
    res = call(cnf, cmdl, exit_on_error=False, silent=False)
    if not res:
        call(cnf, cmdl.replace('-e bcbio', ''), exit_on_error=True, silent=False)
    verify_file(bcbio_structure.multiqc_fpath, is_critical=True)
    return bcbio_structure.multiqc_fpath


def create_rnaseq_pca_plot(cnf, bcbio_structure):
    info('Making RNASeq PCA plot')

    csv_files_in_config_dir = [
        join(bcbio_structure.config_dir, fname)
        for fname in listdir(bcbio_structure.config_dir)
        if fname.endswith('.csv')]
    if not csv_files_in_config_dir:
        info('No CSV file found in config dir ' + bcbio_structure.config_dir)
        return None

    pca_r_script = get_script_cmdline(cnf, 'rscript', join('tools', 'pca.R'), is_critical=True)
    csv_fpath = join(bcbio_structure.config_dir, csv_files_in_config_dir[0])
    gene_counts_fpath = join(bcbio_structure.expression_dirpath, 'counts.tsv')
    gene_counts_fpath = verify_file(gene_counts_fpath, is_critical=True)
    output_fpath = join(bcbio_structure.work_dir, 'pca_data.txt')
    cmdl = pca_r_script + ' ' + csv_fpath + ' ' + output_fpath + ' ' + gene_counts_fpath
    call(cnf, cmdl, output_fpath=output_fpath, stdout_to_outputfile=False)
    return output_fpath


def _add_per_sample_reports(cnf, individual_reports_section, bcbio_structure):
    base_dirpath = dirname(bcbio_structure.multiqc_fpath)

    sample_reports_records = defaultdict(list)

    gender_record_by_sample = dict()
    for s in bcbio_structure.samples:
        gender_fpath = join(s.clinical_report_dirpath, 'gender.txt')
        if isfile(gender_fpath):
            gender = open(gender_fpath).read()
            rec = Record(individual_reports_section.find_metric(GENDER), gender)
            sample_reports_records[s.name].append(rec)

    normal_samples = [s for s in bcbio_structure.samples if s.phenotype == 'normal']
    for s in bcbio_structure.samples:
        if gender_record_by_sample.get(s.name):
            sample_reports_records[s.name].append(gender_record_by_sample.get(s.name))

        if normal_samples:
            rec = Record(individual_reports_section.find_metric(PHENOTYPE), s.phenotype)
            sample_reports_records[s.name].append(rec)
            if s.phenotype != 'normal' and s.normal_match:
                # if len(bcbio_structure.samples) > 1:
                rec = _make_relative_link_record(s.name, s.normal_match.name, individual_reports_section.find_metric(NORM_MATCH))
                # else:
                #     rec = Record(individual_reports_section.find_metric(NORM_MATCH), s.normal_match.name)
                sample_reports_records[s.name].append(rec)

        if bcbio_structure.is_rnaseq:
            pass
        else:
            recs = []
            rec = _make_url_record(verify_file(s.clinical_html, is_critical=True),
                                   individual_reports_section.find_metric(CLINICAL_NAME), base_dirpath)
            if rec and rec.value:
                recs.append(rec)
            sample_reports_records[s.name].extend(recs)

    return sample_reports_records


def _report_to_multiqc_metadata(cnf, full_report, html_fpath, project_name, bcbio_structure,
                                additional_data=None, oncoprints_link=None):
    def __process_record(rec, short=False):
        d = dict()
        d['metric_name'] = rec.metric.name
        if isinstance(rec.url, basestring):
            d['contents'] = '<a href="' + rec.url + '">' + rec.value + '</a>'
        else:
            d['contents'] = rec.metric.format_value(rec.value)
        if rec.metric.description:
            d['description'] = rec.metric.description
        return d

    def _get_summary_report_name(rec):
        return rec.value.lower().replace(' ', '_')

    # common records (summary reports)
    metadata_dict = dict(project_name=project_name)

    metadata_dict['run_section'] = get_run_info(cnf, bcbio_structure)

    mutations_links = []
    expression_links = []
    for rec in full_report.get_common_records():
        if rec.value:
            if rec.metric.name in mutation_names:
                mutations_links.append(__process_record(rec)['contents'])
            if rec.metric.name in expression_names:
                expression_links.append(__process_record(rec)['contents'])
    if oncoprints_link:
        mutations_links.append(('<a href="{oncoprints_link}" target="_blank">oncoprints</a> ' +
                             '(loading may take 5-10 seconds)').format(**locals()))
    if mutations_links:
        metadata_dict["mutations_links"] = mutations_links
    if expression_links:
        metadata_dict["expression_links"] = expression_links

    # SAMPLE_LEVEL DETAILS
    gender_by_sample = dict()
    for s in bcbio_structure.samples:
        gender_fpath = join(s.clinical_report_dirpath, 'gender.txt')
        if isfile(gender_fpath):
            gender = open(gender_fpath).read()
            gender_by_sample[s.name] = gender
    metadata_dict['gender_by_sample'] = gender_by_sample

    base_dirpath = dirname(bcbio_structure.multiqc_fpath)
    metadata_dict['ngs_report_by_sample'] = {s.name:
        relpath(s.clinical_html, base_dirpath)
        for s in bcbio_structure.samples
        if verify_file(s.clinical_html, silent=True)}

    # normal_samples = [s for s in bcbio_structure.samples if s.phenotype == 'normal']
    # for s in bcbio_structure.samples:
    #     if normal_samples:
    #         rec = Record(individual_reports_section.find_metric(PHENOTYPE), s.phenotype)
    #         sample_reports_records[s.name].append(rec)
    #         if s.phenotype != 'normal' and s.normal_match:
    #             # if len(bcbio_structure.samples) > 1:
    #             rec = _make_relative_link_record(s.name, s.normal_match.name, individual_reports_section.find_metric(NORM_MATCH))
    #             # else:
    #             #     rec = Record(individual_reports_section.find_metric(NORM_MATCH), s.normal_match.name)
    #             sample_reports_records[s.name].append(rec)

    if full_report.sample_reports:
        # individual records
        metadata_dict["sample_reports"] = []

        metrics = metric_storage.get_metrics(skip_general_section=True)
        metrics_with_values_set = set()
        for sample_report in full_report.sample_reports:
            for m in metric_storage.get_metrics(skip_general_section=True):
                r = next((r for r in sample_report.records if r.metric.name == m.name), None)
                if r:
                    metrics_with_values_set.add(m)

        metrics = [m for m in metrics if m in metrics_with_values_set]
        metadata_dict['metric_names'] = [m.name for m in metrics]

        for sample_report in full_report.sample_reports:
            ready_records = []
            for m in metrics:
                r = next((r for r in sample_report.records if r.metric.name == m.name), None)
                if r:
                    ready_records.append(__process_record(r, short=True))
                else:
                    ready_records.append(__process_record(Record(metric=m, value=None), short=True))
            assert len(ready_records) == len(metadata_dict["metric_names"])

            sample_report_dict = dict()
            sample_report_dict["records"] = ready_records
            sample_report_dict["sample_name"] = sample_report.display_name
            metadata_dict["sample_reports"].append(sample_report_dict)

    if additional_data:
        for k, v in additional_data.items():
            metadata_dict[k] = v

    return metadata_dict


def get_run_info(cnf, bcbio_structure):
    info('Getting run and codebase information...')
    run_info_dict = dict()
    cur_fpath = abspath(getsourcefile(lambda: 0))
    reporting_suite_dirpath = dirname(dirname(dirname(cur_fpath)))

    run_date = cnf.run_date if cnf.run_date else time.localtime()
    run_info_dict["run_date"] = time.strftime('%d %b %Y, %H:%M (GMT%z)', run_date)

    version = get_version()

    last_modified_datestamp = ''
    try:
        py_fpaths = set()
        for rootpath, dirnames, fnames in os.walk(reporting_suite_dirpath):
            for fn in fnames:
                if fn.endswith('.py'):
                    fpath = abspath(join(rootpath, fn))
                    if isfile(fpath):
                        py_fpaths.add(fpath)
        last_modification_time = max(getmtime(fpath) for fpath in py_fpaths)
    except:
        warn('Cannot get codebase files datestamps')
    else:
        last_modified_datestamp = time.strftime('%d %b %Y, %H:%M (GMT%z)', time.localtime(last_modification_time))

    if last_modified_datestamp or version:
        version_text = 'Reporting Suite '
        if version:
            version_text += 'v.' + version
        if version and last_modified_datestamp:
            version_text += ', '
        if last_modified_datestamp:
            version_text += 'last modified ' + last_modified_datestamp
        run_info_dict['suite_version'] = version_text

    # if bcbio_structure.is_rnaseq:
    base_dirpath = dirname(bcbio_structure.multiqc_fpath)
    program_versions = dict(l.strip().split(',') for l in open(bcbio_structure.program_versions_fpath).readlines())
    program_versions['reporting-suite'] = version
    with open(bcbio_structure.program_versions_fpath, 'w') as f:
        for p, v in sorted(program_versions.items(), key=lambda (k, _): k):
            f.write(p + ',' + v + '\n')
    programs_url = relpath(bcbio_structure.program_versions_fpath, base_dirpath) \
        if verify_file(bcbio_structure.program_versions_fpath) else None
    datas_url = relpath(bcbio_structure.data_versions_fpath, base_dirpath) \
        if verify_file(bcbio_structure.data_versions_fpath) else None

    run_info_dict['program_versions'] = '<a href="{programs_url}">program versions</a>'.format(**locals())
    run_info_dict['data_versions'] = '<a href="{datas_url}">data versions</a>'.format(**locals())
    run_info_dict['analysis_dir'] = bcbio_structure.final_dirpath

    # var_filtering_params = []
    # set_filtering_params(cnf, bcbio_structure=bcbio_structure)  # get variant filtering parameters from run_info yaml
    # dfts = defaults['variant_filtering']
    # for parameter in dfts:
    #     if dfts[parameter] != cnf.variant_filtering[parameter]:
    #         var_filtering_params.append('%s: %s' % (parameter, cnf.variant_filtering[parameter]))
    # if var_filtering_params:
    #     run_info_dict["filtering_params"] = ', '.join(var_filtering_params)
    # else:
    #     run_info_dict["filtering_params"] = 'default'
    return run_info_dict


def get_oncoprints_link(cnf, bcbio_structure, project_name):
    caller = bcbio_structure.variant_callers.get('vardict') or \
             bcbio_structure.variant_callers.get('vardict-java')

    _base_mut_fname = variant_filtering.mut_fname_template.format(caller_name=caller.name)
    _base_mut_fpath = join(bcbio_structure.date_dirpath, _base_mut_fname)
    mut_fpath = add_suffix(_base_mut_fpath, variant_filtering.mut_pass_suffix)
    single_mut_fpath = add_suffix(add_suffix(_base_mut_fpath, variant_filtering.mut_single_suffix), variant_filtering.mut_pass_suffix)
    paired_mut_fpath = add_suffix(add_suffix(_base_mut_fpath, variant_filtering.mut_paired_suffix), variant_filtering.mut_pass_suffix)
    mut_fpath = verify_file(mut_fpath, silent=True)
    single_mut_fpath = verify_file(single_mut_fpath, silent=True)
    paired_mut_fpath = verify_file(paired_mut_fpath, silent=True)
    mutations_fpaths = [f for f in [mut_fpath, single_mut_fpath, paired_mut_fpath] if f]

    return create_oncoprints_link(
        cnf.work_dir, cnf.genome.name, bcbio_structure.final_dirpath, project_name,
        mutations_fpaths, bcbio_structure.seq2c_fpath, bcbio_structure.bed)


def _make_relative_link_record(name, match_name, metric):
    value = '<a class="dotted-link" href="#{match_name}" id="{name}_match">{match_name}</a>'.format(name=name, match_name=match_name)
    return Record(metric=metric, value=value)
