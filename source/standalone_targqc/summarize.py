from copy import deepcopy
import sys
import shutil
import os
from os import listdir
from os.path import relpath, join, exists, dirname, basename
from collections import OrderedDict, defaultdict
from joblib import Parallel, delayed

import source
from source.config import CallCnf
from source.reporting import SampleReport, FullReport, Metric, MetricStorage, ReportSection, write_tsv_rows, load_records, \
    Record, PerRegionSampleReport, Report
from source.logger import step_greetings, info, send_email, critical, warn, err
from source.targetcov import cov
from source.qualimap import report_parser as qualimap_report_parser
from source.ngscat import report_parser as ngscat_report_parser
from source.targetcov.bam_and_bed_utils import count_bed_cols, prepare_beds
from source.tools_from_cnf import get_system_path, get_qualimap_type
from source.calling_process import call
from source.file_utils import safe_mkdir, verify_file, verify_dir, intermediate_fname
from source.bcbio_structure import BCBioStructure
from source.variants.vcf_processing import bgzip_and_tabix
from source.utils import get_numeric_value
from sub_scripts.targetcov import Sample


def _run_multisample_qualimap(cnf, output_dir, samples, targqc_full_report):
    """ 1. Generates Qualimap2 plots and put into plots_dirpath
        2. Adds records to targqc_full_report.plots
    """
    plots_dirpath = join(output_dir, 'plots')

    # Qualimap2 run for multi-sample plots
    if len([s.qualimap_html_fpath for s in samples if s.qualimap_html_fpath]):
        qualimap = get_system_path(cnf, interpreter=None, name='qualimap')

        if qualimap is not None and get_qualimap_type(qualimap) == 'full':
            qualimap_output_dir = join(cnf.work_dir, 'qualimap_multi_bamqc')

            _correct_qualimap_genome_results(samples)
            _correct_qualimap_insert_size_histogram(samples)

            safe_mkdir(qualimap_output_dir)
            rows = []
            for sample in samples:
                if sample.qualimap_html_fpath:
                    rows += [[sample.name, sample.qualimap_html_fpath]]

            data_file = write_tsv_rows(rows, qualimap_output_dir, 'qualimap_results_by_sample')
            cmdline = '{qualimap} multi-bamqc --data {data_file} -outdir {qualimap_output_dir}'.format(**locals())
            ret_code = call(cnf, cmdline, exit_on_error=False, return_err_code=True, env_vars=dict(DISPLAY=None))

            targqc_full_report.plots = []
            qualimap_plots_dirpath = join(qualimap_output_dir, 'images_multisampleBamQcReport')
            if (ret_code is None or ret_code == 0) and verify_dir(qualimap_plots_dirpath):
                if exists(plots_dirpath):
                    shutil.rmtree(plots_dirpath)
                shutil.move(qualimap_plots_dirpath, plots_dirpath)
                for plot_fpath in listdir(plots_dirpath):
                    plot_fpath = join(plots_dirpath, plot_fpath)
                    if verify_file(plot_fpath) and plot_fpath.endswith('.png'):
                        targqc_full_report.plots.append(relpath(plot_fpath, output_dir))
            else:
                warn('Warning: Qualimap for multi-sample analysis failed to finish. TargQC will not contain plots.')
        else:
            warn('Warning: Qualimap for multi-sample analysis was not found. TargQC will not contain plots.')


def _make_targetcov_symlinks(samples):
    for sample in samples:
        new_link = join(
            dirname(dirname(sample.targetcov_detailed_txt)),
            basename(sample.targetcov_detailed_txt))
        if exists(new_link):
            os.unlink(new_link)
        os.symlink(sample.targetcov_detailed_txt, new_link)
        info('TargetCov TXT symlink saved to ' + new_link)


def _make_tarqc_html_report(cnf, output_dir, samples):
    header_storage = cov.get_header_metric_storage(cnf.coverage_reports.depth_thresholds)

    targqc_metric_storage = _get_targqc_metric_storage([
        ('targetcov', header_storage),
        ('ngscat', ngscat_report_parser.metric_storage),
        ('qualimap', qualimap_report_parser.metric_storage)])

    targqc_full_report = FullReport(cnf.name, [], metric_storage=targqc_metric_storage)

    for sample in samples:
        records_by_report_type = []
        if (verify_file(sample.targetcov_json_fpath, True) or
                verify_file(sample.ngscat_html_fpath, True) or
                verify_file(sample.qualimap_html_fpath, True)):
            records_by_report_type.append(('targetcov', load_records(sample.targetcov_json_fpath) if verify_file(
                sample.targetcov_json_fpath, silent=True) else []))
            records_by_report_type.append(('ngscat', ngscat_report_parser.parse_ngscat_sample_report(
                sample.ngscat_html_fpath) if verify_file(sample.ngscat_html_fpath, silent=True) else []))
            records_by_report_type.append(('qualimap', qualimap_report_parser.parse_qualimap_sample_report(
                sample.qualimap_html_fpath) if verify_file(sample.qualimap_html_fpath, silent=True) else []))

        targqc_full_report.sample_reports.append(
            SampleReport(
                sample,
                records=_get_targqc_records(records_by_report_type, header_storage),
                html_fpath=dict(
                    targetcov=relpath(sample.targetcov_html_fpath, output_dir) if sample.targetcov_html_fpath else None,
                    ngscat=relpath(sample.ngscat_html_fpath, output_dir) if sample.ngscat_html_fpath else None,
                    qualimap=relpath(sample.qualimap_html_fpath, output_dir) if sample.qualimap_html_fpath else None
                ),
                metric_storage=targqc_metric_storage
            )
        )

    _run_multisample_qualimap(cnf, output_dir, samples, targqc_full_report)

    txt_fpath = targqc_full_report.save_txt(output_dir, BCBioStructure.targqc_name)
    tsv_fpath = targqc_full_report.save_tsv(output_dir, BCBioStructure.targqc_name)
    html_fpath = targqc_full_report.save_html(output_dir, BCBioStructure.targqc_name,
        'Coverage statistics for all samples based on TargetSeq, ngsCAT, and Qualimap reports')

    return txt_fpath, tsv_fpath, html_fpath


def summarize_targqc(cnf, output_dir, samples, bed_fpath, exons_fpath, genes_fpath=None):
    step_greetings('Coverage statistics for all samples based on TargetSeq, ngsCAT, and Qualimap reports')

    for sample in samples:
        if not sample.targetcov_done():
            sys.exit(1)
        if not sample.ngscat_done():
            sample.ngscat_html_fpath = None
        if not sample.qualimap_done():
            sample.qualimap_html_fpath = None

    _make_targetcov_symlinks(samples)

    best_for_regions_fpath = _save_best_detailed_for_each_gene(cnf.coverage_reports.depth_thresholds, samples, output_dir)

    exons_bed, bed_fpath, _ = prepare_beds(cnf, exons_fpath, bed_fpath)

    # all_htmls_by_sample = OrderedDict()
    # for sample in samples:
    #     all_htmls_by_sample[sample.name] = OrderedDict()
    #     if sample.name in targetcov_htmls_by_sample:
    #         all_htmls_by_sample[sample.name]['targetcov'] = relpath(targetcov_htmls_by_sample[sample.name], output_dir)
    #     if sample.name in ngscat_htmls_by_sample:
    #         all_htmls_by_sample[sample.name]['ngscat'] =    relpath(ngscat_htmls_by_sample[sample.name], output_dir)
    #     if sample.name in qualimap_htmls_by_sample:
    #         all_htmls_by_sample[sample.name]['qualimap'] =  relpath(qualimap_htmls_by_sample[sample.name], output_dir)

    txt_fpath, tsv_fpath, html_fpath = _make_tarqc_html_report(cnf, output_dir, samples)

    norm_best_var_fpath, norm_comb_var_fpath = _report_normalize_coverage_for_variant_sites(cnf, output_dir, samples, 'oncomine', bed_fpath)

    info()
    info('*' * 70)
    info('TargQC summary saved in: ')
    for fpath in [txt_fpath, html_fpath]:
        if fpath: info('  ' + fpath)

    info()
    info('Best stats for regions saved in:')
    info('  ' + best_for_regions_fpath)

    if norm_best_var_fpath:
        info()
        info('Normalized depths for oncomine saved in:')
        info('        ' + norm_comb_var_fpath)
        info('  Best: ' + norm_best_var_fpath)


def get_ave_coverage(cnf, report_fpath):
    if verify_file(report_fpath):
        records = load_records(report_fpath)
        return next((r.value for r in records if r.metric.name == 'Average target coverage depth'), None)


def _clip_vcf_by_bed(cnf, vcf_fpath, bed_fpath):
    info('Clipping VCF ' + vcf_fpath + ' using BED ' + bed_fpath)

    bedtools = get_system_path(cnf, 'bedtools')

    clipped_vcf_fpath = intermediate_fname(cnf, vcf_fpath, 'clip')
    cmdline = '{bedtools} intersect -header	-a {vcf_fpath} -b {bed_fpath}'.format(**locals())
    res = call(cnf, cmdline, output_fpath=clipped_vcf_fpath)

    clipped_gz_vcf_fpath = bgzip_and_tabix(cnf, clipped_vcf_fpath)

    return clipped_gz_vcf_fpath


class Variant:
    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.depth = None


def _get_depth_for_each_variant(cnf_dict, samtools, bedtools, sample_name, bam_fpath, bed_fpath, vcf_fpath):
    info()
    info('Processing sample ' + sample_name)

    # http://www.1000genomes.org/faq/what-depth-coverage-your-phase1-variants
    # bedtools intersect -a oncomine.vcf -b Exons.az_key.bed -header > oncomine.az_key.vcf
    # /opt/az/local/tabix/tabix-0.2.6/bgzip oncomine.az_key.vcf
    # /opt/az/local/tabix/tabix-0.2.6/tabix -h -p vcf oncomine.az_key.vcf.gz
    # samtools view -b TRF004223.sorted.bam -L Exons.az_key.bed | bedtools genomecov -ibam stdin -bg > coverage.bg
    # bedtools intersect -a oncomine.az_key.vcf.gz -b coverage.bg -wa | cut -f1,2,4,5,8,11,12,13,14 > oncomine.az_key.depth_numbers.vcf

    cnf = CallCnf(cnf_dict)

    info()
    info('Depth of coverage for regions in BED ' + bed_fpath)
    cov_bg = join(cnf.work_dir, sample_name + '_coverage.bg')
    cmdline = '{samtools} view -L {bed_fpath} -b {bam_fpath} | {bedtools} genomecov -ibam stdin -bg'.format(**locals())
    call(cnf, cmdline, output_fpath=cov_bg)

    info()
    info('Intersection depth regions with VCF ' + vcf_fpath)
    vcf_depth_numbers_fpath = join(cnf.work_dir, sample_name + '_vcf_bg.intersect')
    cmdline = '{bedtools} intersect -a {vcf_fpath} -b {cov_bg} -wao'.format(**locals())
    res = call(cnf, cmdline, output_fpath=vcf_depth_numbers_fpath)
    # if res != oncomine_depth_numbers_fpath:
    #     info()
    #     info('Trying with uncompressed VCF')
    #     cmdline = 'gunzip {vcf_fpath} -c | {bedtools} intersect -a - -b {cov_bg} -wao | cut -f1,2,4,5,8,11,12,13,14,15'.format(**locals())
    #     call(cnf, cmdline, output_fpath=oncomine_depth_numbers_fpath)

    variants = []
    depths_per_var = defaultdict(list)
    with open(vcf_depth_numbers_fpath) as f:
        for l in f:
            # 1,2,4,5,8,11,12,13,14,15,16,17,18,19,20,21,22
            # c,p,r,a,f,ch,st,en,ge,ex,st,ft,bt,de,ov
            fs = l[:-1].split('\t')
            chrom, pos, _, ref, alt, _, _, info_fields = fs[:8]
            depth, overlap = fs[-2:]
            if any('om_MutClassPC=' + t in info_fields or 'om_MutClass=' + t
                    in info_fields for t in ['Hotspot', 'Deleterious']):
                var = Variant(chrom, pos, ref, alt)
                variants.append(var)
                if depth != '.':
                    depth, overlap = int(depth), int(overlap)
                    for i in range(overlap):
                        depths_per_var[var].append(depth)

    # getting avarage depth of coverage of each variant (exactly for those parts that were in BED)
    var_by_locus = dict()
    for var in variants:
        depths = depths_per_var[var]
        var.depth = (sum(depths) / len(depths)) if len(depths) != 0 else None
        var_by_locus[(var.chrom, var.pos, var.ref, var.alt)] = var

    info()
    info('Intersection BED with VCF to annotate variants with gene names')
    vcf_bed_intersect = join(cnf.work_dir, sample_name + '_vcf_bed.intersect')
    cmdline = '{bedtools} intersect -a {vcf_fpath} -b {bed_fpath} -wo'.format(**locals())
    res = call(cnf, cmdline, output_fpath=vcf_bed_intersect)

    bed_columns_num = count_bed_cols(bed_fpath)

    regions_in_order = []
    regions_set = set()
    vars_by_region = defaultdict(list)
    with open(vcf_bed_intersect) as f:
        for l in f:
            fs = l[:-1].split('\t')
            chrom, pos, _, ref, alt, _, _, info_fields = fs[:8]
            chrom_b, start_b, end_b, symbol, strand, feature, biotype = None, None, None, None, None, None, None
            if bed_columns_num == 8:
                chrom_b, start_b, end_b, symbol, _, strand, feature, biotype, _ = fs[-9:]
            if bed_columns_num == 4:
                chrom_b, start_b, end_b, symbol, _ = fs[-5:]
            assert chrom == chrom_b, l
            r = chrom, start_b, end_b, symbol, strand, feature, biotype
            if r not in regions_set:
                regions_set.add(r)
                regions_in_order.append(r)
            var = var_by_locus.get((chrom, pos, ref, alt))
            if var:
                vars_by_region[r].append(var)

    return sample_name, (vars_by_region, regions_in_order)


def _prep_comb_report(metric_storage, samples, shared_general_metrics, shared_metrics):
    comb_general_metrics = shared_general_metrics[:]
    comb_general_metrics.append(Metric('2 columns for each sample'))
    for s in samples:
        comb_general_metrics.append(Metric(s.name + ' ave depth'))

    comb_metrics = shared_metrics[:]
    for s in samples:
        comb_metrics.append(Metric(s.name + ' depth', short_name='DP'))
        comb_metrics.append(Metric(s.name + ' norm depth', short_name='Norm DP'))

    comb_report_metric_storage = MetricStorage(
        general_section=ReportSection('general_section', metrics=comb_general_metrics),
        sections=[ReportSection(metrics=comb_metrics)])

    report = PerRegionSampleReport(sample='Combined', metric_storage=comb_report_metric_storage)

    report.add_record('Sample', 'contains values from all samples: ' + ', '.join([s.name for s in samples]))
    report.add_record('2 columns for each sample', 'Depth, Normalized depth')

    m = metric_storage.get_metric('Average sample depth')
    for s in samples:
        val = Report.find_record(s.report.records, m.name).value
        report.add_record(s.name + ' ave depth', val)

    return report


def _prep_best_report(metric_storage, samples):
    report = PerRegionSampleReport(sample='Best', metric_storage=metric_storage)

    report.add_record('Sample', 'contains best values from all samples: ' + ', '.join([s.name for s in samples]))

    m = metric_storage.get_metric('Average sample depth')
    ave_sample_depth = max(Report.find_record(s.report.records, m.name).value for s in samples)
    report.add_record('Average sample depth', ave_sample_depth)

    return report


def _report_normalize_coverage_for_variant_sites(cnf, output_dir, samples, vcf_key, bed_fpath):
    step_greetings('Normalized coverage for ' + vcf_key + ' hotspots')
    vcf_fpath = cnf.genome.get(vcf_key)
    if not vcf_fpath:
        err('Error: no ' + vcf_key + ' for ' + cnf.genome.name + ' VCF fpath specified in ' + cnf.sys_cnf)
        return None

    ave_coverages_per_sample = {
        s.name: get_ave_coverage(cnf, s.targetcov_json_fpath)
        for s in samples if verify_file(s.targetcov_json_fpath)}

    vcf_fpath = _clip_vcf_by_bed(cnf, vcf_fpath, bed_fpath)
    samtools = get_system_path(cnf, 'samtools')
    bedtools = get_system_path(cnf, 'bedtools')

    vars_by_region_per_sample = OrderedDict(Parallel(n_jobs=cnf.threads)(delayed(_get_depth_for_each_variant)(
        cnf.__dict__, samtools, bedtools, s.name, s.bam, bed_fpath, vcf_fpath) for s in samples))

    ############################ For each sample ############################
    # Building metrics storages
    shared_general_metrics = [
        Metric('Sample', short_name='Sample', common=True)]

    shared_metrics = [
        Metric('Chr'),
        Metric('Start'),
        Metric('End'),
        Metric('Strand'),
        Metric('Feature'),
        Metric('Biotype'),
        Metric('Symbol'),
        Metric('Var pos'),
        Metric('Ref'),
        Metric('Alt')]

    metric_storage = MetricStorage(
        general_section=ReportSection('general_section', '',
            metrics=shared_general_metrics + [
                Metric('Average sample depth', short_name='Ave depth', common=True)
            ]),
        sections=[ReportSection(metrics=shared_metrics + [
            Metric('Depth'),
            Metric('Norm depth'),
        ])])

    for sample in samples:
        info()
        info('Saving report for sample ' + sample.name)
        report = PerRegionSampleReport(sample=sample, metric_storage=metric_storage)
        sample.report = report
        report.add_record('Sample', sample.name)

        ave_sample_depth = ave_coverages_per_sample[sample.name]
        report.add_record('Average sample depth', ave_sample_depth)

        total_variants = 0
        total_regions = 0
        vars_by_region, regions = vars_by_region_per_sample[sample.name]
        for r in regions:
            total_regions += 1

            (chrom, start, end, symbol, strand, feature, biotype) = r
            rep_region = report.add_region()
            rep_region.add_record('Chr', chrom)
            rep_region.add_record('Start', start)
            rep_region.add_record('End', end)
            rep_region.add_record('Symbol', symbol)
            rep_region.add_record('Strand', strand)
            rep_region.add_record('Feature', feature)
            rep_region.add_record('Biotype', biotype)
            rep_region.add_record('Var pos', None)
            rep_region.add_record('Ref', None)
            rep_region.add_record('Alt', None)
            rep_region.add_record('Depth', None)
            rep_region.add_record('Norm depth', None)

            for i, var in enumerate(vars_by_region[r]):
                total_variants += 1
                if total_variants % 10000 == 0:
                    info('Processed {0:,} variants, {0:,} regions.'.format(total_variants, total_regions))

                rep_region = report.add_region()
                rep_region.add_record('Chr', chrom)
                rep_region.add_record('Start', start)
                rep_region.add_record('End', end)
                rep_region.add_record('Symbol', symbol)
                rep_region.add_record('Strand', strand)
                rep_region.add_record('Feature', feature)
                rep_region.add_record('Biotype', biotype)
                    # rep_region.add_record('Chr', '')
                    # rep_region.add_record('Start', '')
                    # rep_region.add_record('End', '')
                    # rep_region.add_record('Symbol', symbol)
                    # rep_region.add_record('Strand', '')
                    # rep_region.add_record('Feature', '')
                    # rep_region.add_record('Biotype', '')

                rep_region.add_record('Var pos', var.pos)
                rep_region.add_record('Ref', var.ref)
                rep_region.add_record('Alt', var.alt)
                rep_region.add_record('Depth', var.depth)
                rep_region.add_record('Norm depth', var.depth / ave_sample_depth if var.depth and ave_sample_depth > 0 else None)

        best_report_basename = sample.name + '.' + source.targetseq_name  + '_' + vcf_key
        sample.targetcov_norm_depth_vcf_txt = report.save_txt(sample.dirpath, best_report_basename)
        sample.targetcov_norm_depth_vcf_tsv = report.save_tsv(sample.dirpath, best_report_basename)
        info('')
        info('Oncomine variants coverage report (total: {0:,} variants, {0:,} regions) saved into:'.format(total_variants, total_regions))
        info('  ' + sample.targetcov_norm_depth_vcf_txt)

    ############################ Best ############################
    info()
    info('*' * 70)
    info('Saving for all samples: combined and best values.')
    best_report = _prep_best_report(metric_storage, samples)
    comb_report = _prep_comb_report(metric_storage, samples, shared_general_metrics, shared_metrics)

    total_variants = 0
    nth_regions_from_each_sample = [s.report.get_regions() for s in samples]
    while True:
        nth_region_from_each_sample = [rs[total_variants] for rs in nth_regions_from_each_sample if total_variants < len(rs)]
        total_variants += 1
        if len(nth_region_from_each_sample) == 0:
            break
        assert len(nth_region_from_each_sample) == len(nth_regions_from_each_sample), 'Region files for samples are not euqal size'

        best_report_reg = best_report.add_region()
        comb_report_reg = comb_report.add_region()
        rand_line = nth_region_from_each_sample[0]
        for i in range(10):
            best_report_reg.records.append(rand_line.records[i])
            comb_report_reg.records.append(rand_line.records[i])

        best_depth = select_best(r.records[10].value for r in nth_region_from_each_sample)
        best_norm_depth = select_best(r.records[11].value for r in nth_region_from_each_sample)
        best_report_reg.add_record('Depth', best_depth)
        best_report_reg.add_record('Norm depth', best_norm_depth)

        for s, r in zip(samples, nth_region_from_each_sample):
            comb_report_reg.add_record(s.name + ' depth', r.records[10].value)
            comb_report_reg.add_record(s.name + ' norm depth', r.records[11].value)

    best_report_basename = 'Best.' + source.targetseq_name  + '_' + vcf_key
    comb_report_basename = 'Comb.' + source.targetseq_name  + '_' + vcf_key
    best_targetcov_norm_depth_vcf_txt = best_report.save_txt(output_dir, best_report_basename)
    best_targetcov_norm_depth_vcf_tsv = best_report.save_tsv(output_dir, best_report_basename)
    comb_targetcov_norm_depth_vcf_txt = comb_report.save_txt(output_dir, comb_report_basename)
    comb_targetcov_norm_depth_vcf_tsv = comb_report.save_tsv(output_dir, comb_report_basename)
    info('')
    info('Depths for Oncomine variants (total: {0:,} variants, {0:,} regions) saved into:'.format(total_variants))
    info('  Best:     ' + best_targetcov_norm_depth_vcf_txt)
    info('  Combined: ' + comb_targetcov_norm_depth_vcf_txt)

    return best_targetcov_norm_depth_vcf_txt, comb_targetcov_norm_depth_vcf_txt


def _get_targqc_metric(metric, header_metric_storage, report_type='targetcov'):  # report type is in ['targetcov', 'qualimap', 'ngscat']
    qualimap_to_targetcov_dict = {
        'Number of reads': header_metric_storage.get_metric('Reads'),
        'Mapped reads': header_metric_storage.get_metric('Mapped reads'),
        'Unmapped reads': header_metric_storage.get_metric('Unmapped reads'),
        'Mapped reads (on target)': header_metric_storage.get_metric('Reads mapped on target'),
        'Coverage Mean': header_metric_storage.get_metric('Average target coverage depth'),
        'Coverage Standard Deviation': header_metric_storage.get_metric('Std. dev. of target coverage depth')
    }

    ngscat_to_targetcov_dict = {
        'Number reads': header_metric_storage.get_metric('Mapped reads'),
        # '% target bases with coverage >= 1x': cov.header_metric_storage.get_metric('Percentage of target covered by at least 1 read'),
        '% reads on target': header_metric_storage.get_metric('Reads mapped on target'),
        'mean coverage': header_metric_storage.get_metric('Average target coverage depth')
    }

    if report_type == 'targetcov':
        return metric
    elif report_type == 'qualimap':
        if metric.name in qualimap_to_targetcov_dict:
            return qualimap_to_targetcov_dict[metric.name]
        return metric
    elif report_type == 'ngscat':
        if metric.name in ngscat_to_targetcov_dict:
            return ngscat_to_targetcov_dict[metric.name]
        return metric
    # critical('Incorrect usage of get_targqc_metric(), report_type is %s but should be one of the following: %s' %
    #          (report_type, ", ".join(['targetcov', 'qualimap', 'ngscat'])))
    return metric


def _get_targqc_metric_storage(metric_storages_by_report_type):
    class SectionId:
        def __init__(self, name, title):
            self.name = name
            self.title = title

        def __hash__(self):
            #return hash((self.name, self.title))
            return hash(self.name)  # use title from the first metric_storage

        def __eq__(self, other):
            #return (self.name, self.title) == (other.name, other.title)
            return self.name == other.name  # use title from the first metric_storage

    metrics_by_sections = OrderedDict()
    general_section_id = None
    general_section_metric_list = []

    for report_type, metric_storage in metric_storages_by_report_type:
        for section in metric_storage.sections:
            section_id = SectionId(section.name, section.title)
            if section_id not in metrics_by_sections.keys():
                metrics_by_sections[section_id] = []

            metrics_by_sections[section_id] += [metric
                for metric in metric_storage.get_metrics(sections=[section], skip_general_section=True)
                if metric == _get_targqc_metric(metric, dict(metric_storages_by_report_type)['targetcov'], report_type)]

        # specific behaviour for general section
        general_section_metric_list += [metric
            for metric in metric_storage.general_section.metrics
            if metric == _get_targqc_metric(metric, dict(metric_storages_by_report_type)['targetcov'], report_type)]
        if not general_section_id:
            general_section_id = SectionId(metric_storage.general_section.name, metric_storage.general_section.title)

    sections = []
    for section_id, metric_list in metrics_by_sections.items():
        sections.append(ReportSection(section_id.name, section_id.title, metric_list))

    return MetricStorage(
        general_section=ReportSection(
            general_section_id.name, general_section_id.title, general_section_metric_list),
        sections=sections)


def _get_targqc_records(records_by_report_type, header_storage):
    targqc_records = []
    filled_metric_names = []
    for report_type, records in records_by_report_type:
        for record in records:
            new_metric = _get_targqc_metric(record.metric, header_storage, report_type)
            if not new_metric or new_metric.name not in filled_metric_names:
                filled_metric_names.append(new_metric.name)
                record.metric = new_metric
                targqc_records.append(record)
    return targqc_records


def _correct_qualimap_genome_results(samples):
    """ fixing java.lang.Double.parseDouble error on entries like "6,082.49"
    """
    for s in samples:
        if verify_file(s.qualimap_genome_results_fpath):
            with open(s.qualimap_genome_results_fpath, 'r') as f:
                content = f.readlines()
            with open(s.qualimap_genome_results_fpath, 'w') as f:
                metrics_started = False
                for line in content:
                    if ">> Reference" in line:
                        metrics_started = True
                    if metrics_started:
                        line = line.replace(',', '')
                    f.write(line)


def _correct_qualimap_insert_size_histogram(samples):
    """ replacing Qualimap IS histogram with Picard one.
    """
    for s in samples:
        qualimap1_dirname = dirname(s.qualimap_ins_size_hist_fpath).replace('raw_data_qualimapReport', 'raw_data')
        qualimap2_dirname = dirname(s.qualimap_ins_size_hist_fpath)
        if exists(qualimap1_dirname):
            if not exists(qualimap2_dirname):
                shutil.move(qualimap1_dirname, qualimap2_dirname)
            else:
                shutil.rmtree(qualimap1_dirname)
        elif not exists(qualimap2_dirname):
            continue  # no data from both Qualimap v.1 and Qualimap v.2

        if verify_file(s.picard_ins_size_hist_fpath):
            with open(s.picard_ins_size_hist_fpath, 'r') as picard_f:
                one_line_to_stop = False
                for line in picard_f:
                    if one_line_to_stop:
                        break
                    if line.startswith('## HISTOGRAM'):
                        one_line_to_stop = True
                with open(s.qualimap_ins_size_hist_fpath, 'w') as qualimap_f:
                    for line in picard_f:
                        qualimap_f.write(line)
        elif verify_file(s.qualimap_ins_size_hist_fpath):
            os.remove(s.qualimap_ins_size_hist_fpath)


def select_best(values, fn=max):
    vs = [v for v in values if v is not None]
    return fn(vs) if len(vs) > 0 else None


def _save_best_detailed_for_each_gene(depth_threshs, samples, output_dir):
    metric_storage = cov.get_detailed_metric_storage(depth_threshs)

    report = PerRegionSampleReport(sample='Best', metric_storage=metric_storage)
    report.add_record('Sample', 'contains best values from all samples: ' + ', '.join([s.name for s in samples]))

    def get_int_val(v):
        v = _get_num(v)
        return int(v) if v else None

    def get_float_val(v):
        v = _get_num(v)
        return float(v) if v else None

    def _get_num(v):
        v = get_val(v)
        return ''.join(c for c in v if c.isdigit() or c == '.') if v else None

    def get_val(v):
        return v.strip() if v.strip() not in ['.', '-', ''] else None

    total_regions = 0
    open_tsv_files = [open(s.targetcov_detailed_tsv) for s in samples]

    first_col = 0
    while True:
        lines_for_each_sample = [next(f, None) for f in open_tsv_files]
        if not all(lines_for_each_sample):
            break
        l = lines_for_each_sample[0]
        if l.startswith('##'):
            continue
        elif l.startswith('#'):
            if l.startswith('#Sample'):
                first_col = 1
            break

    while True:
        lines_for_each_sample = [next(f, None) for f in open_tsv_files]
        if not all(lines_for_each_sample):
            break

        if all([not l.startswith('#') and 'Whole-Gene' in l for l in lines_for_each_sample]):
            shared_fields = lines_for_each_sample[0].split('\t')[first_col:first_col+8]
            reg = report.add_region()
            reg.add_record('Chr', get_val(shared_fields[0]))
            reg.add_record('Start', get_int_val(shared_fields[1]))
            reg.add_record('End', get_int_val(shared_fields[2]))
            reg.add_record('Size', get_int_val(shared_fields[3]))
            reg.add_record('Gene', get_val(shared_fields[4]))
            reg.add_record('Strand', get_val(shared_fields[5]))
            reg.add_record('Feature', get_val(shared_fields[6]))
            reg.add_record('Biotype', get_val(shared_fields[7]))

            min_depths, ave_depths, stddevs, withins = ([], [], [], [])
            percents_by_threshs = {t: [] for t in depth_threshs}

            for l in lines_for_each_sample:
                fs = l.split('\t')

                min_depths.append(get_int_val(fs[first_col+8]))
                ave_depths.append(get_float_val(fs[first_col+9]))
                stddevs.append(get_float_val(fs[first_col+10]))
                withins.append(get_float_val(fs[first_col+11]))
                for t, f in zip(depth_threshs, fs[first_col+12:]):
                    percents_by_threshs[t].append(get_float_val(f))

            # counting bests
            reg.add_record('Min depth', select_best(min_depths))
            reg.add_record('Ave depth', select_best(ave_depths))
            reg.add_record('Std dev', select_best(stddevs, max))
            reg.add_record('W/n 20% of ave depth', select_best(withins))
            for t in depth_threshs:
                reg.add_record('{}x'.format(t), select_best(percents_by_threshs[t]))

            total_regions += 1

    for f in open_tsv_files:
        f.close()

    gene_report_basename = 'Best.' + source.targetseq_name + source.targetcov.detail_gene_report_baseending
    txt_rep_fpath = report.save_txt(output_dir, gene_report_basename)
    tsv_rep_fpath = report.save_tsv(output_dir, gene_report_basename)
    info('')
    info('Best values for the regions (total ' + str(total_regions) + ') saved into:')
    info('  ' + txt_rep_fpath)

    return txt_rep_fpath