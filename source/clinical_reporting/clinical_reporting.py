from collections import OrderedDict
from json import load, dump, JSONEncoder, dumps
import os
from os.path import join, islink, dirname, abspath, relpath

import source
from source import verify_file, info
from source.file_utils import add_suffix, verify_module
from source.file_utils import adjust_path
from source.html_reporting.html_saver import write_static_html_report
from source.reporting import MetricStorage, Metric, PerRegionSampleReport, ReportSection, SampleReport
from source.targetcov.flag_regions import get_depth_cutoff
from source.targetcov.summarize_targetcov import get_float_val, get_val
from tools import seq2c_plots


def make_key_gene_cov_report(cnf, sample, key_gene_names, ave_depth):
    info('Preparing coverage stats key gene tables')

    depth_cutoff = get_depth_cutoff(ave_depth, cnf.coverage_reports.depth_thresholds)

    stats_by_genename = dict()
    with open(sample.targetcov_detailed_tsv) as f_inp:
        for l in f_inp:
            if not l.startswith('#') and ('Whole-Gene' in l or 'Gene-Exon' in l):
                fs = l.split('\t')
                gene_name = get_val(fs[4])
                if gene_name in key_gene_names:
                    for t, field in zip(cnf.coverage_reports.depth_thresholds, fs[12:]):
                        if int(t) == depth_cutoff:
                            stats_by_genename[gene_name] = get_float_val(fs[9]), get_float_val(field)
                            continue

    clinical_cov_metrics = [
            Metric('Gene'),
            Metric('Mean coverage'),
            Metric('% covered at {}x'.format(depth_cutoff), unit='%')]
    seq2c_tsv = cnf.seq2c_tsv_fpath
    seq2c_data_by_genename = dict()
    if verify_file(seq2c_tsv, silent=True):
        with open(seq2c_tsv) as f_inp:
            for i, l in enumerate(f_inp):
                if i == 0:
                    continue
                fs = l.strip().split('\t')
                gene_name = fs[1]
                if fs[0] == sample.name and gene_name in key_gene_names and fs[9] in ['Del', 'Amp']:
                    seq2c_data_by_genename[gene_name] = fs[9] + ',' + fs[8]
        if seq2c_data_by_genename:
            clinical_cov_metrics.append(Metric('SNV'))

    clinical_cov_metric_storage = MetricStorage(
        sections=[ReportSection(metrics=clinical_cov_metrics)])
    key_genes_report = PerRegionSampleReport(sample=sample, metric_storage=clinical_cov_metric_storage)
    for gene_name in key_gene_names:
        ave_cov, depth_in_thresh = stats_by_genename.get(gene_name, (None, None))
        reg = key_genes_report.add_region()
        reg.add_record('Gene', gene_name)
        reg.add_record('Mean coverage', ave_cov)
        reg.add_record('% covered at {}x'.format(depth_cutoff), depth_in_thresh)
        if seq2c_data_by_genename:
            reg.add_record('SNV', seq2c_data_by_genename[gene_name] if gene_name in seq2c_data_by_genename else '')

    key_genes_report.save_tsv(sample.clinical_targqc_tsv, human_readable=True)
    info('Saved coverage report to ' + key_genes_report.tsv_fpath)
    info('-' * 70)
    info()
    return key_genes_report


def get_target_fraction(sample, targqc_json_fpath):
    with open(targqc_json_fpath) as f:
        data = load(f, object_pairs_hook=OrderedDict)
    sr = SampleReport.load(data, sample, None)
    r = sr.find_record(sr.records, 'Percentage of target covered by at least 1 read')
    if not r:
        r = sr.find_record(sr.records, 'Percentage of genome covered by at least 1 read')
    return r.value if r else None


def get_gender(sample, targqc_json_fpath):
    with open(targqc_json_fpath) as f:
        data = load(f, object_pairs_hook=OrderedDict)
    sr = SampleReport.load(data, sample, None)
    r = sr.find_record(sr.records, 'Gender')
    return r.value if r else None


def get_ave_coverage(sample, targqc_json_fpath):
    with open(targqc_json_fpath) as f:
        data = load(f, object_pairs_hook=OrderedDict)
    sr = SampleReport.load(data, sample, None)
    r = sr.find_record(sr.records, 'Average target coverage depth')
    if not r:
        r = sr.find_record(sr.records, 'Average genome coverage depth')
    return r.value if r else None


def get_total_variants_number(sample, varqc_json_fpath):
    with open(varqc_json_fpath) as f:
        data = load(f, object_pairs_hook=OrderedDict)
    sr = SampleReport.load(data, sample, None)
    r = sr.find_record(sr.records, 'Total variants')
    return r.value if r else None


def is_sample_presents_in_file(sample_name, mutations_fpath):
    with open(mutations_fpath) as f:
        for i, l in enumerate(f):
            if i == 0:
                continue
            fs = l.strip().split('\t')
            if fs[0] == sample_name:
                return True
    return False

def make_mutations_report(cnf, sample, key_gene_names, mutations_fpath):
    info('Preparing mutations stats for key gene tables')

    clinical_mut_metric_storage = MetricStorage(
        sections=[ReportSection(metrics=[
            Metric('Gene'),  # Gene & Transcript
            Metric('Transcript'),  # Gene & Transcript
            Metric('Variant'),            # c.244G>A, p.Glu82Lys
            Metric('Allele'),             # Het.
            Metric('Genomic Position'),       # hg19 chr11:g.47364249G>A
            Metric('Depth'),              # 658
            Metric('Frequency'),          # .19
            Metric('AA length'),          # 128
            Metric('ID'),                 # rs352343, COSM2123
            Metric('Type'),               # Frameshift
            Metric('Classification'),     # Likely Pathogenic
        ])])
    report = PerRegionSampleReport(sample=sample, metric_storage=clinical_mut_metric_storage)
    if not verify_file(mutations_fpath, silent=True):
        single_mutations_fpath = add_suffix(mutations_fpath, source.mut_single_suffix)
        paired_mutations_fpath = add_suffix(mutations_fpath, source.mut_paired_suffix)
        if verify_file(single_mutations_fpath, silent=True) and is_sample_presents_in_file(sample.name, single_mutations_fpath):
            mutations_fpath = single_mutations_fpath
        elif verify_file(paired_mutations_fpath, silent=True):
            mutations_fpath = paired_mutations_fpath

    info('Reading mutations from ' + mutations_fpath)
    with open(mutations_fpath) as f:
        for i, l in enumerate(f):
            if i == 0:
                continue
            fs = l.strip().split('\t')
            if len(fs) > 60:
                sample_name, chrom, start, id, ref, alt, type_, effect, func, codon_change, aa_change, cdna_change, \
                    aa_len, gene, transcr_biotype, coding, transcript, exon, cosmic_cds_change, cosmic_aa_change, \
                    cosmic_cnt, end, depth, af, bias, pmean, pstd, qual, qstd, sbf, gmaf, vd, clnsif, oddratio, hiaf, \
                    mq, sn, adjaf, nm, shift3, msi, dbsnpbuildid, vtype, status1, paired_pval, paired_oddratiom, \
                    m_depth, m_af, m_vd, m_bias, m_pmean, m_pstd, m_qual, m_qstd, m_hiaf, m_mq, m_sn, m_adjaf, m_nm, \
                    n_sample, n_var, pcnt_sample, ave_af, filter_, var_type, var_class, status = fs[:67]  # 67 of them
            else:
                sample_name, chrom, start, id, ref, alt, type_, effect, func, codon_change, aa_change, cdna_change, \
                    aa_len, gene, transcr_biotype, coding, transcript, exon, cosmic_cds_change, cosmic_aa_change, \
                    cosmic_cnt, end, depth, af, bias, pmean, pstd, qual, qstd, sbf, gmaf, vd, clnsif, oddratio, hiaf, \
                    mq, sn, adjaf, nm, shift3, msi, dbsnpbuildid, \
                    n_sample, n_var, pcnt_sample, ave_af, filter_, var_type, var_class, status = fs[:50]  # 50 of them

            if sample_name == sample.name and gene in key_gene_names:
                reg = report.add_region()
                reg.add_record('Gene', gene)
                reg.add_record('Transcript', transcript)
                reg.add_record('Variant', codon_change + (' p.' + aa_change if aa_change else ''))
                reg.add_record('Allele', None)
                reg.add_record('Genomic Position', str(cnf.genome) + ' ' + chrom + ':g.' +
                                   (Metric.format_value(int(start), human_readable=True) if start else '') +
                                   ' ' + ref + '>' + alt)
                reg.add_record('Depth', depth)
                reg.add_record('Frequency', af)
                reg.add_record('AA length', aa_len)
                reg.add_record('ID', ' '.join(id.split(';')))
                reg.add_record('Type', type_[0] + type_[1:].lower().replace('_', ' ') if type_ else type_)
                if status == 'likely':
                    status += ' pathogenic'
                reg.add_record('Classification', status[0].upper() + status[1:] if status else status)

    report.save_tsv(sample.clinical_mutation_tsv, human_readable=True)
    info('Saved mutations report to ' + report.tsv_fpath)
    info('-' * 70)
    info()
    return report


def make_clinical_html_report(cnf, sample, coverage_report, mutations_report,
                              ave_depth, target_frac, gender, total_variants, total_key_genes, key_gene_names):
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
    general_dict = dict()
    general_dict['project_name'] = cnf.project_name
    general_dict['sample'] = sample.name
    general_dict['sex'] = gender
    general_dict['project_name'] = cnf.project_name

    mutations_dict = dict()
    mutations_dict['first_col_header'] = mutations_report.metric_storage.get_metrics()[0].name
    mutations_dict['metric_names'] = [m.name for m in mutations_report.metric_storage.get_metrics()[1:]]
    mutations_dict['rows'] = [
        dict(first_col=region.records[0].value, records=[
                Metric.format_value(r.value) for r in region.records[1:]])
            for region in mutations_report.regions]

    coverage_dict = dict()
    coverage_dict['first_col_header'] = coverage_report.metric_storage.get_metrics()[0].name
    coverage_dict['metric_names'] = [m.name for m in coverage_report.metric_storage.get_metrics()[1:]]
    coverage_dict['rows'] = [
        dict(first_col=region.records[0].value, records=[
                Metric.format_value(r.value) for r in region.records[1:]])
            for region in coverage_report.regions]

    # if full_report.sample_reports:
    #     # individual records
    #     main_dict = dict()
    #     main_dict["sample_reports"] = []
    #
    #     metrics = metric_storage.get_metrics(skip_general_section=True)
    #     metrics_with_values_set = set()
    #     for sample_report in full_report.sample_reports:
    #         for m in metric_storage.get_metrics(skip_general_section=True):
    #             r = next((r for r in sample_report.records if r.metric.name == m.name), None)
    #             if r:
    #                 metrics_with_values_set.add(m)
    #
    #     metrics = [m for m in metrics if m in metrics_with_values_set]
    #     main_dict['metric_names'] = [m.name for m in metrics]
    #
    #     for sample_report in full_report.sample_reports:
    #         ready_records = []
    #         for m in metrics:
    #             r = next((r for r in sample_report.records if r.metric.name == m.name), None)
    #             if r:
    #                 ready_records.append(__process_record(r, short=True))
    #             else:
    #                 ready_records.append(__process_record(Record(metric=m, value=None), short=True))
    #         assert len(ready_records) == len(main_dict["metric_names"])
    #
    #         sample_report_dict = dict()
    #         sample_report_dict["records"] = ready_records
    #         sample_report_dict["sample_name"] = sample_report.get_display_name()
    #         main_dict["sample_reports"].append(sample_report_dict)
    #
    # regions_dict = dict()

    if cnf.seq2c_tsv_fpath:
        if verify_module('matplotlib'):
            seq2c_plot_fpath = seq2c_plots._draw_seq2c_plot(cnf, key_gene_names)
        else:
            seq2c_plot_fpath = None
        plots = [relpath(seq2c_plot_fpath, cnf.output_dir)]
    else:
        plots = None
    sample.clinical_html = write_static_html_report(cnf.work_dir, {
        'general': general_dict,
        'variants': mutations_dict,
        'coverage': coverage_dict,
    }, sample.clinical_html, tmpl_fpath=join(dirname(abspath(__file__)), 'report.html'), plots=plots)

    clin_rep_symlink = adjust_path(join(sample.dirpath, '..', sample.name + '.clinical_report.html'))
    # if islink(clin_rep_symlink):
    #     os.unlink(clin_rep_symlink)
    # os.symlink(sample.clinical_html, clin_rep_symlink)

    info('Saved clinical report to ' + clin_rep_symlink)
    info('-' * 70)
    info()
    return clin_rep_symlink


def tooltip_long(string, max_len=30):
    if len(string) < max_len:
        return string
    else:
        return '<a class="tooltip-link" rel="tooltip" title="' + string + '">' + string[:max_len - 2] + '...</a>'



