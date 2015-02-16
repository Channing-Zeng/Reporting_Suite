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
    Record, PerRegionSampleReport
from source.logger import step_greetings, info, send_email, critical, warn, err
from source.targetcov import cov
from source.qualimap import report_parser as qualimap_report_parser
from source.ngscat import report_parser as ngscat_report_parser
from source.tools_from_cnf import get_system_path, get_qualimap_type
from source.calling_process import call
from source.file_utils import safe_mkdir, verify_file, verify_dir, intermediate_fname
from source.bcbio_structure import BCBioStructure
from source.variants.vcf_processing import bgzip_and_tabix


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
    targetcov_metric_storage = cov.header_metric_storage
    for depth in cnf.coverage_reports.depth_thresholds:
        name = 'Part of target covered at least by ' + str(depth) + 'x'
        targetcov_metric_storage.add_metric(
            Metric(name, short_name=str(depth) + 'x', description=name, unit='%'),
            'depth_metrics')

    targqc_metric_storage = _get_targqc_metric_storage([
        ('targetcov', targetcov_metric_storage),
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
                records=_get_targqc_records(records_by_report_type),
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
    html_fpath = targqc_full_report.save_html(output_dir, BCBioStructure.targqc_name,
        'Coverage statistics for all samples based on TargetSeq, ngsCAT, and Qualimap reports')

    return txt_fpath, html_fpath


def summarize_targqc(cnf, output_dir, samples, bed_fpath):
    step_greetings('Coverage statistics for all samples based on TargetSeq, ngsCAT, and Qualimap reports')

    for sample in samples:
        if not sample.targetcov_done():
            sys.exit(1)
        if not sample.ngscat_done():
            sample.ngscat_html_fpath = None
        if not sample.qualimap_done():
            sample.qualimap_html_fpath = None

    _make_targetcov_symlinks(samples)

    best_for_regions_fpath = _save_best_detailed_for_each_gene(samples, output_dir)

    # all_htmls_by_sample = OrderedDict()
    # for sample in samples:
    #     all_htmls_by_sample[sample.name] = OrderedDict()
    #     if sample.name in targetcov_htmls_by_sample:
    #         all_htmls_by_sample[sample.name]['targetcov'] = relpath(targetcov_htmls_by_sample[sample.name], output_dir)
    #     if sample.name in ngscat_htmls_by_sample:
    #         all_htmls_by_sample[sample.name]['ngscat'] =    relpath(ngscat_htmls_by_sample[sample.name], output_dir)
    #     if sample.name in qualimap_htmls_by_sample:
    #         all_htmls_by_sample[sample.name]['qualimap'] =  relpath(qualimap_htmls_by_sample[sample.name], output_dir)

    txt_fpath, html_fpath = _make_tarqc_html_report(cnf, output_dir, samples)

    _report_normalize_coverage_and_hotspots(cnf, output_dir, samples, bed_fpath)

    info()
    info('*' * 70)
    info('TargQC summary saved in: ')
    for fpath in [txt_fpath, html_fpath]:
        if fpath: info('  ' + fpath)

    info()
    info('Best stats for regions saved in:')
    info('  ' + best_for_regions_fpath)

    # info()
    # info('Normalized cov for oncomine saved in:')
    # info('  ' + norm_oncomine_report_fpath)


def get_ave_coverage(cnf, report_fpath):
    if verify_file(report_fpath):
        records = load_records(report_fpath)
        return next((r.value for r in records if r.metric.name == 'Average target coverage depth'), None)


def _clip_oncomine_vcf_by_bed(cnf, bed_fpath):
    oncomine_vcf_fpath = cnf.genomes[cnf.genome].oncomine
    info('Clipping Oncomine VCF ' + oncomine_vcf_fpath)

    bedtools = get_system_path(cnf, 'bedtools')

    oncomine_clipped_vcf_fpath = intermediate_fname(cnf, oncomine_vcf_fpath, 'clip')
    cmdline = '{bedtools} intersect -a {oncomine_vcf_fpath} -b {bed_fpath}'.format(**locals())
    res = call(cnf, cmdline, output_fpath=oncomine_clipped_vcf_fpath)

    oncomine_clipped_gz_vcf_fpath = bgzip_and_tabix(cnf, oncomine_clipped_vcf_fpath)

    return oncomine_clipped_gz_vcf_fpath


class Variant:
    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.depth = None


def _get_depth_for_each_variant(cnf_dict, samtools, bedtools, sample_name, bam_fpath, bed_fpath, vcf_fpath):
    # http://www.1000genomes.org/faq/what-depth-coverage-your-phase1-variants
    # bedtools intersect -a oncomine.vcf -b Exons.az_key.bed -header > oncomine.az_key.vcf
    # /opt/az/local/tabix/tabix-0.2.6/bgzip oncomine.az_key.vcf
    # /opt/az/local/tabix/tabix-0.2.6/tabix -h -p vcf oncomine.az_key.vcf.gz
    # samtools view -b TRF004223.sorted.bam -L Exons.az_key.bed | bedtools genomecov -ibam stdin -bg > coverage.bg
    # bedtools intersect -a oncomine.az_key.vcf.gz -b coverage.bg -wa | cut -f1,2,4,5,8,11,12,13,14 > oncomine.az_key.depth_numbers.vcf

    cnf = CallCnf(cnf_dict)

    cov_bg = join(cnf.work_dir, sample_name + '_coverage.bg')
    cmdline = '{samtools} view -b {bam_fpath} -L {bed_fpath} | {bedtools} genomecov -ibam stdin -bg'.format(**locals())
    call(cnf, cmdline, output_fpath=cov_bg)

    oncomine_depth_numbers_fpath = intermediate_fname(cnf, vcf_fpath[:-3], 'depth_numbers')
    cmdline = '{bedtools} intersect -a {vcf_fpath} -b {cov_bg} -wao | cut -f1,2,4,5,8,11,12,13,14,15'.format(**locals())
    call(cnf, cmdline, output_fpath=oncomine_depth_numbers_fpath)

    variants = []
    depths_per_var = defaultdict(list())
    with open(oncomine_depth_numbers_fpath) as f:
        for l in f:
            # 1,2,4,5,8,11,12,13,14,15,16,17,18,19,20,21,22
            # c,p,r,a,f,ch,st,en,ge,ex,st,ft,bt,de,ov
            chrom, pos, ref, alt, fields, _, _, _, depth, overlap = l[:-1].split('\t')
            var = Variant(chrom, pos, ref, alt)
            variants.append(var)
            depth = int(depth)
            for i in range(int(overlap)):
                depths_per_var[var].append(depth)

    # getting avarage depth of coverage of each variant (exactly for those parts that were in BED)
    for var in variants:
        depths = depths_per_var[var]
        var.depth = sum(depths) / len(depths)

    return variants


def _report_normalize_coverage_and_hotspots(cnf, output_dir, samples, bed_fpath):
    ave_coverages_per_sample = {
        s.name: get_ave_coverage(cnf, s.targetcov_json_fpath)
        for s in samples if verify_file(s.targetcov_json_fpath)
    }

    oncomine_vcf_fpath = _clip_oncomine_vcf_by_bed(cnf, bed_fpath)
    samtools = get_system_path(cnf, 'samtools')
    bedtools = get_system_path(cnf, 'bedtools')

    def proc(s):
        info()
        info('Processing sample ' + s.name)
        return s.name, _get_depth_for_each_variant(cnf.__dict__, samtools, bedtools, s.name, s.bam, bed_fpath, oncomine_vcf_fpath)

    variants_per_sample = dict(Parallel(n_jobs=cnf.threads)(delayed(proc)(s) for s in samples))

    metric_storage = MetricStorage(
        general_section=ReportSection('general_section', '', [
            Metric('Sample name', short_name='Sample', common=True),
            Metric('Average sample depth', short_name='Ave depth', common=True),
        ]),
        sections=[ReportSection(metrics=[
            Metric('Sample'),
            Metric('Chr'),
            Metric('Pos'),
            Metric('Ref'),
            Metric('Alt'),
            Metric('Depth'),
            Metric('Normalized depth')
        ])])

    for sample in samples:
        info()
        info('Saving report for sample ' + sample.name)
        report = PerRegionSampleReport(sample=sample, metric_storage=metric_storage)

        ave_sample_depth = ave_coverages_per_sample[sample.name]
        report.add_record('Sample name', sample.name)
        report.add_record('Average sample depth', ave_sample_depth)

        i = 0
        for var in variants_per_sample[sample.name]:
            i += 1
            if i % 10000 == 0:
                info('Processed {0:,} regions.'.format(i))
            rep_region = report.add_region()
            rep_region.add_record('Chr', var.chrom)
            rep_region.add_record('Pos', var.pos)
            rep_region.add_record('Ref', var.ref)
            rep_region.add_record('Alt', var.alt)
            rep_region.add_record('Depth', var.depth)
            rep_region.add_record('Normalized depth', var.depth / ave_sample_depth if ave_sample_depth > 0 else None)

        report_basename = sample.name + '_oncomine'
        txt_rep_fpath = report.save_txt(output_dir, report_basename)
        tsv_rep_fpath = report.save_tsv(output_dir, report_basename)
        info('')
        info('Oncomine variants coverage report (total' + str(len(variants_per_sample)) + ') saved into:')
        info('  ' + txt_rep_fpath)


_qualimap_to_targetcov_dict = {
    'Number of reads': cov.header_metric_storage.get_metric('Reads'),
    'Mapped reads': cov.header_metric_storage.get_metric('Mapped reads'),
    'Unmapped reads': cov.header_metric_storage.get_metric('Unmapped reads'),
    'Mapped reads (on target)': cov.header_metric_storage.get_metric('Reads mapped on target'),
    'Coverage Mean': cov.header_metric_storage.get_metric('Average target coverage depth'),
    'Coverage Standard Deviation': cov.header_metric_storage.get_metric('Std. dev. of target coverage depth')}

_ngscat_to_targetcov_dict = {
    'Number reads': cov.header_metric_storage.get_metric('Mapped reads'),
    # '% target bases with coverage >= 1x': cov.header_metric_storage.get_metric('Percentage of target covered by at least 1 read'),
    '% reads on target': cov.header_metric_storage.get_metric('Reads mapped on target'),
    'mean coverage': cov.header_metric_storage.get_metric('Average target coverage depth')}


def _get_targqc_metric(metric, report_type='targetcov'):  # report type is in ['targetcov', 'qualimap', 'ngscat']
    if report_type == 'targetcov':
        return metric
    elif report_type == 'qualimap':
        if metric.name in _qualimap_to_targetcov_dict:
            return _qualimap_to_targetcov_dict[metric.name]
        return metric
    elif report_type == 'ngscat':
        if metric.name in _ngscat_to_targetcov_dict:
            return _ngscat_to_targetcov_dict[metric.name]
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
                if metric == _get_targqc_metric(metric, report_type)]

        # specific behaviour for general section
        general_section_metric_list += [metric
            for metric in metric_storage.general_section.metrics
            if metric == _get_targqc_metric(metric, report_type)]
        if not general_section_id:
            general_section_id = SectionId(metric_storage.general_section.name, metric_storage.general_section.title)

    sections = []
    for section_id, metric_list in metrics_by_sections.items():
        sections.append(ReportSection(section_id.name, section_id.title, metric_list))

    return MetricStorage(
        general_section=ReportSection(
            general_section_id.name, general_section_id.title, general_section_metric_list),
        sections=sections)


def _get_targqc_records(records_by_report_type):
    targqc_records = []
    filled_metric_names = []
    for report_type, records in records_by_report_type:
        for record in records:
            new_metric = _get_targqc_metric(record.metric, report_type)
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


def _save_best_detailed_for_each_gene(samples, output_dir):
    best_metrics_fpath = join(output_dir, 'Best.targetSeq.details.gene.tsv')
    with open(best_metrics_fpath, 'w') as best_f:
        with open(samples[0].targetcov_detailed_tsv) as f:
            header_fields = f.readline().split('\t')
            # 0        1      2      3    4     5     6       7        8        9          10         11        12                  13...
            # #Sample, Chrom, Start, End, Size, Gene, Strand, Feature, Biotype, Min Depth, Avg Depth, Std Dev., Within 20% of Mean, 1x...

        threshold_num = len(header_fields[13:])

        best_f.write('Chr\tStart\tEnd\tSize\tSymbol\tStrand\tFeature\tBiotype\t'
                     'Min depth\tAvg depth\tStd dev.\tW/n 20% of ave\t' + '\t'.join(header_fields[13:]))

        open_tsv_files = [open(s.targetcov_detailed_tsv) for s in samples]
        while True:
            lines_for_each_sample = [next(f, None) for f in open_tsv_files]
            if not all(lines_for_each_sample):
                break

            if all([not l.startswith('#') and 'Whole-Gene' in l for l in lines_for_each_sample]):
                shared_fields = lines_for_each_sample[0].split('\t')[1:9]
                best_f.write('\t'.join(shared_fields) + '\t')  # chrom, start, end, symbol, strand, feature, biotype, size

                min_depths, ave_depths, stddevs, withins = ([], [], [], [])
                percents_by_col_num = {t: [] for t in range(threshold_num)}

                for l in lines_for_each_sample:
                    fs = l.split('\t')

                    val = ''.join(c for c in fs[9] if c.isdigit())  # skipping decimal dots and spaces
                    if val not in ['', '.', '-']:
                        min_depths.append(int(val))

                    val = fs[10]
                    if val not in ['', '.', '-']:
                        ave_depths.append(float(val))

                    val = fs[11]
                    if val not in ['', '.', '-']:
                        stddevs.append(float(val))

                    val = fs[12][:-1]  # skipping "%" sign
                    if val not in ['', '.', '-']:
                        withins.append(float(val))

                    for col_num, f in enumerate(fs[13:]):
                        val = f.strip()[:-1]  # skipping "%" sign
                        if val not in ['', '.', '-']:
                            percents_by_col_num[col_num].append(float(val))

                # counting bests
                min_depth = max(min_depths) if min_depths else '.'
                ave_depth = max(ave_depths) if ave_depths else '.'
                stddev = min(stddevs) if stddevs else '.'
                within = max(withins) if withins else '.'
                percent_by_col_num = dict()
                for col_num, values in percents_by_col_num.iteritems():
                    percent_by_col_num[col_num] = max(values) if values else '.'

                best_f.write('{:,}\t'.format(min_depth) if min_depth != '.' else '.\t')
                best_f.write('{:.2f}\t'.format(ave_depth) if ave_depth != '.' else '.\t')
                best_f.write('{:.2f}\t'.format(stddev) if stddev != '.' else '.\t')
                best_f.write('{:.2f}%\t'.format(within) if within != '.' else '.\t')
                for col_num, val in percent_by_col_num.iteritems():
                    best_f.write('{:.2f}%\t'.format(val) if val != '.' else '.\t')
                best_f.write('\n')

        for f in open_tsv_files:
            f.close()

    return best_metrics_fpath