import json
import os
import time
from genericpath import isdir, exists
from inspect import getsourcefile
from os import listdir
from os.path import join, relpath, dirname, basename, abspath, getmtime, isfile, pardir
from collections import OrderedDict
from collections import defaultdict

import shutil
import variant_filtering

import source
from ngs_reporting.oncoprints import create_oncoprints_link
from source.bcbio.bcbio_structure import BCBioStructure
from source.calling_process import call
from source.logger import info, step_greetings, warn, timestamp
from source.file_utils import verify_file, add_suffix, verify_dir, file_transaction
from source.reporting.reporting import Metric, Record, MetricStorage, ReportSection, SampleReport, FullReport, \
    write_static_html_report
from source.tools_from_cnf import get_system_path, get_script_cmdline
from source.utils import get_ext_tools_dirname, get_version
from tools.prepare_data_for_exac import get_exac_us_url


BASECALLS_NAME        = 'BaseCalls'
FASTQC_NAME           = BCBioStructure.fastqc_repr
PRE_FASTQC_NAME       = 'Raw ' + FASTQC_NAME
SEQQC_NAME            = 'Targ QC'
EXAC_NAME             = 'ExAC browser'
PRE_SEQQC_NAME        = 'Pre ' + SEQQC_NAME
VARQC_NAME            = 'Var QC'
VARQC_AFTER_NAME      = 'Var QC after filtering'
MUTATIONS_NAME        = 'Mutations'
MUTATIONS_SINGLE_NAME = 'Mutations for single samples'
MUTATIONS_PAIRED_NAME = 'Mutations for paired samples'
GENDER                = 'Sex'
CLINICAL_NAME         = 'Oncology NGS report'
PHENOTYPE             = 'Phenotype'
NORM_MATCH            = 'Normal Match'
ABNORMAL_NAME         = 'Flagged regions'
MULTIQC_NAME          = 'MultiQC'

## RNAseq reports
QC_REPORT_NAME        = 'QC report'
GENE_COUNTS_NAME      = 'Gene counts'
EXON_COUNTS_NAME      = 'Exon counts'
GENE_TPM_NAME         = 'Gene TPM'
ISOFORM_TPM_NAME      = 'Isoform TPM'
QUALIMAP_NAME         = 'Qualimap'

metric_storage = MetricStorage(
    general_section=ReportSection(metrics=[
        Metric(BASECALLS_NAME),
        Metric(PRE_FASTQC_NAME),
        Metric(FASTQC_NAME),
        Metric(PRE_SEQQC_NAME),
        Metric(SEQQC_NAME),
        Metric(EXAC_NAME),
        Metric(VARQC_NAME),
        Metric(VARQC_AFTER_NAME),
        Metric(MUTATIONS_NAME),
        Metric(MUTATIONS_SINGLE_NAME),
        Metric(MUTATIONS_PAIRED_NAME),
        Metric(ABNORMAL_NAME),
        Metric(QC_REPORT_NAME),
        Metric(GENE_COUNTS_NAME),
        Metric(EXON_COUNTS_NAME),
        Metric(GENE_TPM_NAME),
        Metric(ISOFORM_TPM_NAME),
        Metric(MULTIQC_NAME),
    ]),
    sections=[ReportSection(metrics=[
        Metric(PRE_FASTQC_NAME),
        Metric(FASTQC_NAME),
        # Metric('BAM'),
        Metric(PRE_SEQQC_NAME),
        Metric(SEQQC_NAME),
        Metric(QUALIMAP_NAME),
        # Metric(MUTATIONS_NAME),
        Metric(VARQC_NAME),
        Metric(VARQC_AFTER_NAME),
        Metric(GENDER, description='If not defined, means that the target does not contain key male Y genes that we could check'),
        Metric(CLINICAL_NAME),
        Metric(GENE_COUNTS_NAME),
        Metric(PHENOTYPE),
        Metric(NORM_MATCH),
    ])])


def make_project_level_report(cnf, dataset_structure=None, bcbio_structure=None, dataset_project=None,
                              oncoprints_link=None):

    step_greetings('Making the %s project-level report' % ('preproc' if bcbio_structure is None else 'postproc'))

    # if dataset_structure is None and bcbio_structure:
    #     analysis_dirpath = normpath(join(bcbio_structure.bcbio_project_dirpath, pardir))
    #     dataset_dirpath = realpath(join(analysis_dirpath, 'dataset'))
    #     dataset_structure = DatasetStructure.create(dataset_dirpath, bcbio_structure.project_name)

    bcbio_multiqc_available = bcbio_structure and isfile(bcbio_structure.multiqc_fpath)
    if bcbio_multiqc_available:
        info('Found MultiQC report ' + bcbio_structure.multiqc_fpath)

    general_records = _add_summary_reports(cnf, bcbio_multiqc_available, metric_storage.general_section, bcbio_structure, dataset_structure, dataset_project)
    sample_reports_records = _add_per_sample_reports(cnf, bcbio_multiqc_available, metric_storage.sections[0], bcbio_structure, dataset_structure, dataset_project)

    sample_reports = []
    samples = []
    if dataset_project:
        samples = dataset_project.sample_by_name.values()
    if bcbio_structure:
        samples = bcbio_structure.samples
    for sample in samples:
        sample_reports.append(SampleReport(sample,
            records=sample_reports_records[sample.name],
            html_fpath=None,
            metric_storage=metric_storage))

    full_report = FullReport(cnf.project_name, sample_reports,
                             metric_storage=metric_storage, general_records=general_records)

    project_report_html_fpath = None
    project_name = None
    if dataset_structure and dataset_project:
        project_report_html_fpath = dataset_project.project_report_html_fpath
        project_name = dataset_project.name
    if bcbio_structure:
        project_report_html_fpath = bcbio_structure.project_report_html_fpath
        project_name = bcbio_structure.project_name

    sample_match_on_hover_js = None
    if bcbio_structure:
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

    _save_static_html(cnf, full_report, project_report_html_fpath, project_name, bcbio_structure,
                      additional_data=dict(sample_match_on_hover_js=sample_match_on_hover_js),
                      oncoprints_link=oncoprints_link, dataset_project=dataset_project)

    info()
    info('*' * 70)
    info('Project-level report saved in: ')
    info('  ' + project_report_html_fpath)
    return project_report_html_fpath


def _mutations_records(general_section, bcbio_structure):
    records = []

    val = OrderedDict()
    single_val = OrderedDict()
    paired_val = OrderedDict()

    for caller in bcbio_structure.variant_callers.values():
        _base_mut_fname = variant_filtering.mut_fname_template.format(caller_name=caller.name)
        _base_mut_fpath = join(bcbio_structure.date_dirpath, _base_mut_fname)
        mut_fpath = add_suffix(_base_mut_fpath, variant_filtering.mut_pass_suffix)
        if verify_file(mut_fpath, silent=True):
            val[caller.name] = mut_fpath
        else:
            single_mut_fpath = add_suffix(add_suffix(_base_mut_fpath, variant_filtering.mut_single_suffix), variant_filtering.mut_pass_suffix)
            paired_mut_fpath = add_suffix(add_suffix(_base_mut_fpath, variant_filtering.mut_paired_suffix), variant_filtering.mut_pass_suffix)
            if verify_file(single_mut_fpath, silent=True):
                single_val[caller.name] = single_mut_fpath
                # _add_rec(single_mut_fpath, caller.name + ' mutations for separate samples')
            if verify_file(paired_mut_fpath, silent=True):
                paired_val[caller.name] = paired_mut_fpath
                # _add_rec(paired_mut_fpath, caller.name + ' mutations for paired samples')

    for val, metric_name in (
         (val, MUTATIONS_NAME),
         (single_val, MUTATIONS_SINGLE_NAME),
         (paired_val, MUTATIONS_PAIRED_NAME)):
        if val:
            metric = Metric(metric_name, common=True)
            rec = Record(
                metric=metric,
                value=metric.name,
                url=_relpath_all(val, bcbio_structure.date_dirpath))
            general_section.add_metric(metric)
            records.append(rec)

    return records


def _make_url_record(html_fpath_value, metric, base_dirpath):
    # info('Adding paths to the report: ' + str(html_fpath_value))
    if isinstance(html_fpath_value, dict):
        url = OrderedDict([(k, relpath(html_fpath, base_dirpath)) for k, html_fpath in html_fpath_value.items() if verify_file(html_fpath)])
        return Record(metric=metric, value=metric.name, url=url)
    else:
        url = relpath(html_fpath_value, base_dirpath) if verify_file(html_fpath_value) else None
        return Record(metric=metric, value=metric.name, url=url)


def get_base_dirpath(bcbio_structure, dataset_project):
    if bcbio_structure:
        return bcbio_structure.date_dirpath
    return dirname(dataset_project.project_report_html_fpath)

def _add_summary_reports(cnf, bcbio_multiqc_available, general_section, bcbio_structure=None, dataset_structure=None, dataset_project=None):
    """ We want links to be relative, so we make paths relative to the project-level-report parent directory.
        - If the bcbio_structure is set, project-level report is located at bcbio_structure.date_dirpath
        - If dataset_dirpath is set, project-level report is located right at dataset_dirpath
    """
    base_dirpath = get_base_dirpath(bcbio_structure, dataset_project)

    recs = []

    if dataset_structure and dataset_project:
        if dataset_structure.basecall_stat_html_reports:
            val = OrderedDict([(basename(fpath), fpath) for fpath in dataset_structure.basecall_stat_html_reports])
            recs.append(_make_url_record(val, general_section.find_metric(BASECALLS_NAME), base_dirpath))
        recs.append(_make_url_record(dataset_project.comb_fastqc_fpath,              general_section.find_metric(PRE_FASTQC_NAME), base_dirpath))
        recs.append(_make_url_record(dataset_project.downsample_targqc_report_fpath, general_section.find_metric(PRE_SEQQC_NAME),  base_dirpath))

    if bcbio_structure:
        if bcbio_multiqc_available:
            recs.append(_make_url_record(bcbio_structure.multiqc_fpath, general_section.find_metric(MULTIQC_NAME), base_dirpath))
        else:
            warn('MultiQC report ' + bcbio_structure.multiqc_fpath + ' not found')
            if not bcbio_structure.is_rnaseq:
                recs = add_dna_summary_records(cnf, recs, general_section, bcbio_structure, base_dirpath)
            # if isfile(bcbio_structure.fastqc_summary_fpath):
            #     recs.append(_make_url_record(bcbio_structure.fastqc_summary_fpath, general_section.find_metric(FASTQC_NAME), base_dirpath))

        if bcbio_structure.is_rnaseq:
            recs = add_rna_summary_records(cnf, recs, general_section, bcbio_structure, base_dirpath)
        else:
            recs.append(_make_url_record(bcbio_structure.targqc_summary_fpath, general_section.find_metric(SEQQC_NAME), base_dirpath))
        # if verify_dir(bcbio_structure.flagged_regions_dirpath, is_critical=False):
        #     url_val = OrderedDict(
        #             [(region_type, join(bcbio_structure.flagged_regions_dirpath, 'flagged_' + region_type + '.html'))
        #                 for region_type in ['target', 'exons']])
        #     rec = _make_url_record(url_val, general_section.find_metric(ABNORMAL_NAME), base_dirpath)
        #     recs.append(rec)

    return recs


def add_rna_summary_records(cnf, recs, general_section, bcbio_structure, base_dirpath):
    recs.append(_make_url_record(bcbio_structure.gene_counts_report_fpath, general_section.find_metric(GENE_COUNTS_NAME), base_dirpath))
    recs.append(_make_url_record(bcbio_structure.exon_counts_report_fpath, general_section.find_metric(EXON_COUNTS_NAME), base_dirpath))
    recs.append(_make_url_record(bcbio_structure.gene_tpm_report_fpath, general_section.find_metric(GENE_TPM_NAME), base_dirpath))
    recs.append(_make_url_record(bcbio_structure.isoform_tpm_report_fpath, general_section.find_metric(ISOFORM_TPM_NAME), base_dirpath))

    # rnaseq_html_fpath = join(bcbio_structure.date_dirpath, BCBioStructure.rnaseq_qc_report_name + '.html')
    # rnaseq_html_fpath = verify_file(rnaseq_html_fpath, is_critical=True)
    # recs.append(_make_url_record(rnaseq_html_fpath, general_section.find_metric(QC_REPORT_NAME), base_dirpath))

    return recs


def add_dna_summary_records(cnf, recs, general_section, bcbio_structure, base_dirpath):
    recs.extend(_mutations_records(general_section, bcbio_structure))

    m = general_section.find_metric(EXAC_NAME)
    recs.append(Record(metric=m, value=m.name, url=get_exac_us_url(cnf.genome.name, bcbio_structure.project_name)))

    # varqc_d = bcbio_structure.varqc_report_fpath_by_caller
    # varqc_d['all'] = bcbio_structure.varqc_report_fpath
    #
    # varqc_after_d = bcbio_structure.varqc_after_report_fpath_by_caller
    # varqc_after_d['all'] = bcbio_structure.varqc_after_report_fpath

    recs.append(_make_url_record(bcbio_structure.targqc_summary_fpath, general_section.find_metric(SEQQC_NAME),  base_dirpath))
    # recs.append(_make_url_record(varqc_d,       general_section.find_metric(VARQC_NAME),       base_dirpath))
    # recs.append(_make_url_record(varqc_after_d, general_section.find_metric(VARQC_AFTER_NAME), base_dirpath))
    return recs


def make_multiqc_report(cnf, bcbio_structure, oncoprints_link):
    multiqc_bcbio_dirpath = join(bcbio_structure.date_dirpath, 'multiqc')
    new_multiqc_bcbio_dirpath = join(bcbio_structure.date_dirpath, 'qc', 'multiqc_bcbio')
    if isdir(multiqc_bcbio_dirpath):
        if isdir(new_multiqc_bcbio_dirpath):
            try:
                shutil.rmtree(new_multiqc_bcbio_dirpath)
            except OSError:
                os.rename(new_multiqc_bcbio_dirpath, new_multiqc_bcbio_dirpath + '.' + timestamp().replace(':', '_').replace(' ', '_'))
        os.rename(multiqc_bcbio_dirpath, new_multiqc_bcbio_dirpath)
    multiqc_bcbio_dirpath = new_multiqc_bcbio_dirpath

    multiqc_postproc_dirpath = join(bcbio_structure.date_dirpath, 'multiqc_postproc')

    to_run = False
    cmdl = 'multiqc -f -o ' + multiqc_postproc_dirpath
    if cnf.debug:
        cmdl += ' -v'

    input_list_fpath = join(bcbio_structure.date_dirpath, 'qc', 'list_files.txt')
    if not verify_file(input_list_fpath, silent=True):
        work_input_list_fpath = join(bcbio_structure.work_dir, pardir, basename(input_list_fpath))
        if verify_file(work_input_list_fpath, silent=True):
            shutil.copy(work_input_list_fpath, input_list_fpath)

    if verify_file(input_list_fpath, silent=True):
        to_run = True
        cmdl += ' -l ' + input_list_fpath

    if bcbio_structure.is_rnaseq:
        if not to_run:
            for s in bcbio_structure.samples:
                for path in [
                    join(s.dirpath, 'qc', 'samtools'),
                    join(s.dirpath, 'qc', 'qualimap_rnaseq'),
                    join(s.dirpath, 'qc', 'fastqc'),
                    join(multiqc_bcbio_dirpath, 'report', 'metrics', s.name + '_bcbio.txt'),
                ]:
                    if exists(path):
                        cmdl += ' ' + path
                        to_run = True

        pca_plot_fpath = create_rnaseq_pca_plot(cnf, bcbio_structure)
        info()
        if pca_plot_fpath and verify_file(pca_plot_fpath):
            cmdl += ' ' + pca_plot_fpath
            to_run = True

    else:
        if not isfile(input_list_fpath):
            return join(multiqc_bcbio_dirpath, 'multiqc_report.html')
        else:
            targqc_dirpath = join(bcbio_structure.date_dirpath, BCBioStructure.targqc_dir)
            if isdir(targqc_dirpath):
                info('Adding TargQC to MultiQC run')
                lines = open(input_list_fpath).readlines()
                lines.append(targqc_dirpath)
                with open(input_list_fpath, 'w') as f:
                    f.writelines(lines)
                cmdl += ' -e qualimap'
                to_run = True

    if to_run:
        call(cnf, cmdl, exit_on_error=False)
        multiqc_fpath = join(multiqc_postproc_dirpath, 'multiqc_report.html')
        verify_file(multiqc_fpath, is_critical=True)
        return multiqc_fpath


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
    gene_counts_fpath = verify_file(bcbio_structure.gene_counts_fpath, is_critical=True, description='Gene counts')
    output_fpath = join(bcbio_structure.work_dir, 'pca_data.txt')
    cmdl = pca_r_script + ' ' + csv_fpath + ' ' + output_fpath + ' ' + gene_counts_fpath
    call(cnf, cmdl, output_fpath=output_fpath, stdout_to_outputfile=False)
    return output_fpath


def _add_per_sample_reports(cnf, bcbio_multiqc_available, individual_reports_section, bcbio_structure=None, dataset_structure=None, dataset_project=None):
    base_dirpath = get_base_dirpath(bcbio_structure, dataset_project)

    sample_reports_records = defaultdict(list)

    if dataset_project:
        for s in dataset_project.sample_by_name.values():
            sample_reports_records[s.name].append(_make_url_record(
                    OrderedDict([('left', s.find_fastqc_html(s.l_fastqc_base_name)), ('right', s.find_fastqc_html(s.r_fastqc_base_name))]),
                    individual_reports_section.find_metric(PRE_FASTQC_NAME), base_dirpath))
            # sample_reports_records[s.name].append(_make_url_record(
            #         OrderedDict([('targqc', s.targetcov_html_fpath), ('qualimap', s.qualimap_html_fpath)]),
            #         individual_reports_section.find_metric(PRE_SEQQC_NAME), base_dirpath))

    if bcbio_structure:
        gender_record_by_sample = dict()

        # for s in bcbio_structure.samples:
        #     if verify_file(s.targetcov_json_fpath):
        #         targqc_json = json.loads(open(s.targetcov_json_fpath).read(), object_pairs_hook=OrderedDict)
        #         sample_report = SampleReport.load(targqc_json, s, bcbio_structure)
        #         gender_rec = sample_report.find_record(sample_report.records, GENDER)
        #         if gender_rec and gender_rec.value:
        #             gender_record_by_sample[s.name] = gender_rec

        # if not gender_record_by_sample:
        #     individual_reports_section.

        normal_samples = [s for s in bcbio_structure.samples if s.phenotype == 'normal']
        for s in bcbio_structure.samples:
            if not bcbio_multiqc_available:
                if s.fastqc_html_fpath and isfile(s.fastqc_html_fpath):
                    sample_reports_records[s.name].append(_make_url_record(s.fastqc_html_fpath, individual_reports_section.find_metric(FASTQC_NAME), base_dirpath))

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
                sample_reports_records[s.name].extend(add_rna_sample_records(s, individual_reports_section, bcbio_structure, base_dirpath))
            else:
                sample_reports_records[s.name].extend(add_dna_sample_records(bcbio_multiqc_available, s, individual_reports_section, bcbio_structure, base_dirpath))

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


def add_rna_sample_records(s, individual_reports_section, bcbio_structure, base_dirpath):
    recs = []
    recs.append(_make_url_record(s.gene_counts, individual_reports_section.find_metric(GENE_COUNTS_NAME), base_dirpath))
    if verify_file(s.qualimap_html_fpath):
        recs.append(_make_url_record(s.qualimap_html_fpath, individual_reports_section.find_metric(QUALIMAP_NAME), base_dirpath))
    return recs


def add_dna_sample_records(bcbio_multiqc_available, s, individual_reports_section, bcbio_structure, base_dirpath):
    recs = []
    # if not bcbio_multiqc_available and not s.phenotype or s.phenotype != 'normal':
    #     varqc_d = OrderedDict([(k, s.get_varqc_fpath_by_callername(k)) for k in bcbio_structure.variant_callers.keys()])
    #     varqc_after_d = OrderedDict([(k, s.get_varqc_after_fpath_by_callername(k)) for k in bcbio_structure.variant_callers.keys()])
    #     recs.extend([
    #         _make_url_record(varqc_d,             individual_reports_section.find_metric(VARQC_NAME),       base_dirpath),
    #         _make_url_record(varqc_after_d,       individual_reports_section.find_metric(VARQC_AFTER_NAME), base_dirpath),
    #     ])

    if not verify_file(s.clinical_html, is_critical=False):
        clinical_html = join(dirname(dirname(s.clinical_html)), 'clinicalReport', basename(s.clinical_html))
        if verify_file(clinical_html):
            s.clinical_html = clinical_html
    rec = _make_url_record(s.clinical_html, individual_reports_section.find_metric(CLINICAL_NAME), base_dirpath)
    if rec and rec.value:
        recs.append(rec)
    return recs


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


def _save_static_html(cnf, full_report, html_fpath, project_name, bcbio_structure,
                      additional_data=None, oncoprints_link=None, dataset_project=None):
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

        if isinstance(rec.url, basestring):
            d['contents'] = '<a href="' + rec.url + '">' + rec.value + '</a>'

        elif isinstance(rec.url, dict):
            d['contents'] = ', '.join('<a href="{v}">{k}</a>'.format(k=k, v=v) for k, v in rec.url.items()) if rec.url else '-'
            if not short:
                d['contents'] = rec.metric.name + ': ' + d['contents']

        else:
            d['contents'] = rec.metric.format_value(rec.value)

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
    common_dict['run_section'] = get_run_info(cnf, bcbio_structure, dataset_project)

    if oncoprints_link:
        common_dict['oncoprints'] = {'oncoprints_link': '<a href="{oncoprints_link}" target="_blank">Oncoprints</a> ' \
                                                       '(loading may take 5-10 seconds)'.format(**locals())}

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
            sample_report_dict["sample_name"] = sample_report.display_name
            main_dict["sample_reports"].append(sample_report_dict)

    data = {"common": common_dict, "main": main_dict}
    if additional_data:
        for k, v in additional_data.items():
            data[k] = v

    return write_static_html_report(cnf, data, html_fpath)


def get_run_info(cnf, bcbio_structure, dataset_project):
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
    base_dirpath = get_base_dirpath(bcbio_structure, dataset_project)
    if bcbio_structure:
        programs_url = relpath(bcbio_structure.program_versions_fpath, base_dirpath) if verify_file(bcbio_structure.program_versions_fpath) else None
        run_info_dict['program_versions'] = '<a href="{programs_url}">Program versions</a>'.format(**locals())

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
    oncoprints_link = create_oncoprints_link(cnf, bcbio_structure, project_name)
    return oncoprints_link


def _make_relative_link_record(name, match_name, metric):
    value = '<a class="dotted-link" href="#{match_name}" id="{name}_match">{match_name}</a>'.format(name=name, match_name=match_name)
    return Record(metric=metric, value=value)
