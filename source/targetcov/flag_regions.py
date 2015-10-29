# coding=utf-8

from collections import defaultdict
from os.path import join, splitext, basename

from source.calling_process import call
from source.file_utils import intermediate_fname
from source.logger import info, err, step_greetings
from source.reporting.reporting import Metric, MetricStorage, ReportSection, PerRegionSampleReport, load_records
import source
from source.targetcov.Region import save_regions_to_bed
from source.targetcov.bam_and_bed_utils import count_bed_cols
from source.targetcov.cov import get_detailed_metric_storage, make_flat_region_report
from source.tools_from_cnf import get_system_path
from source.variants.vcf_processing import bgzip_and_tabix

DEPTH_THRESH_FROM_AVE_COV = 0.5
MIN_DEPTH_PERCENT_AT_THRESH = 0.8


def get_depth_cutoff(ave_depth, depth_thresholds):
    depth_cutoff = ave_depth * DEPTH_THRESH_FROM_AVE_COV
    for thresh in depth_thresholds[::-1]:
        if thresh < depth_cutoff:
            return thresh
    return 1


def get_ave_sample_coverage(cnf, report_fpath):
    records = load_records(report_fpath)
    return next((r.value for r in records if r.metric.name == 'Average target coverage depth'), None)


def generate_flagged_regions_report(cnf, output_dir, sample, ave_depth, gene_by_name):
    depth_threshs = cnf.coverage_reports.depth_thresholds
    report = PerRegionSampleReport(sample=sample, metric_storage=get_detailed_metric_storage(depth_threshs))
    report.add_record('Sample', sample.name)

    ''' 1. Detect depth threshold (ave sample coverage * DEPTH_THRESH_FROM_AVE_COV)
        2. Select regions covered in less than MIN_DEPTH_PERCENT_AT_THRESH at threshold
        3. Sort by % at threshold
        4. Select those parts of those regions where % = 0, save to BED
        5. Find HotSpots at those regions
        6. Intersect HotSpots with tracks

        For each gene where are regions with parts % = 0:
            sort them by part where % = 0
    '''

    depth_cutoff = get_depth_cutoff(ave_depth, depth_threshs)

    # for gene in gene_by_name.values():
    #     total_size = 0
    #     total_bases_at_thresh = 0
    #     for r in gene.get_exons():
    #         total_size += r.get_size()
    #         total_bases_at_thresh += r.bases_within_threshs[depth_cutoff]
    #     if total_size:
    #         gene = float(total_bases_at_thresh) / total_size

    genes_sorted = sorted(gene_by_name.values(), key=lambda g: [g.rates_within_threshs[t] for t in depth_threshs])

    info('Selecting and saving low-covereg genes')
    low_cov_genes = [g for g in genes_sorted if
        any(e.rates_within_threshs[depth_cutoff] < MIN_DEPTH_PERCENT_AT_THRESH for e in g.get_exons()) or \
        any(a.rates_within_threshs[depth_cutoff] < MIN_DEPTH_PERCENT_AT_THRESH for a in g.get_amplicons())]

    final_regions = []
    for gene in low_cov_genes:
        final_regions.extend(gene.get_amplicons())
        final_regions.extend(gene.get_exons())
        final_regions.append(gene)

    selected_regions_bed_fpath = join(cnf.work_dir, 'selected_regions.bed')
    save_regions_to_bed(cnf, final_regions, selected_regions_bed_fpath)
    # TODO: extract only those subregions where cov is 0

    # Roport cov for Hotspots
    _report_normalize_coverage_for_variant_sites(
        cnf, output_dir, sample, ave_depth, 'oncomine', selected_regions_bed_fpath, depth_cutoff)

    report = make_flat_region_report(sample, final_regions, depth_threshs)

    report.save_tsv(sample.flagged_regions_tsv)
    info('')
    info('Selected regions (total ' + str(len(final_regions)) + ') saved into:')
    info('  ' + report.txt_fpath)

    return report


# sample_depths = [100, 120, 110, 130, 100, 90, 110, 200]
# med_sample_depth = median(sample_depths)
# mad = median(abs(med_sample_depth - v) for v in sample_depths)

# print med_sample_depth
# print mad

# z_scores = [(d - med_sample_depth) / mad for d in sample_depths]

# for d, z in zip(sample_depths, z_scores):
# 	print '{d}: {z}'.format(**locals())

##########################################
################# RORY ###################
##########################################
# import numpy
# import scipy.stats as stats

# def outlier_resistant(coverage):
#   median = float(numpy.median(coverage))
#   deviations = [abs(x - median) for x in coverage]
#   # median absolute deviation estimator of the standard deviation
#   mad = 1.4826 * float(numpy.median(deviations))
#   return int(median), int(mad)
#
# def modified_zscore(coverage):
#   median, mad = outlier_resistant(coverage)
#   z_scores = [float(x - median) / mad for x in coverage]
#   return z_scores
#
# example = [20,30,40,50,10,10,11,12,13,14,5,100,300,20,30,40]
#
# print modified_zscore(example)
#
# p_values = stats.norm.sf(modified_zscore(example))
#
# def keep_outliers(p_values, cutoff=0.05):
#     return [(x, y) for x, y in zip(example, p_values) if y < cutoff]
#
# print keep_outliers(p_values)

# you can see this flags only a couple of the samples as outliers with p-value < 0.05

# def FDR(x):
#     """
#     Copied from p.adjust function from R (ripped from
#     http://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python)
#     """
#     o = [i[0] for i in sorted(enumerate(x), key=lambda v:v[1],reverse=True)]
#     ro = [i[0] for i in sorted(enumerate(o), key=lambda v:v[1])]
#     q = sum([1.0/i for i in xrange(1, len(x)+1)])
#     l = [q*len(x)/i*x[j] for i, j in zip(reversed(xrange(1, len(x)+1)),o)]
#     l = [l[k] if l[k] < 1.0 else 1.0 for k in ro]
#     return l
#
# # the 50 one has a pretty marginal p-value though, 0.01. if we FDR adjust it and
# # filter on FDR < 0.05:
#
# fdr_values = FDR(p_values)
#
# print keep_outliers(fdr_values)
#
# def questionable_coverages(coverage, min_coverage=15, fdr_cutoff=0.05):
#     fdr_values = FDR(stats.norm.sf(modified_zscore(coverage)))
#     return [x for x, y in zip(coverage, fdr_values) if x < min_coverage or y < fdr_cutoff]
#
# print questionable_coverages(example, 10)
#
# def original_implementation(coverage):
#     coverage = sorted(coverage)
#     q1 = coverage[4]
#     q3 = coverage[12]
#     d = abs(q1 - q3)
#     bottom = q1 - 3 * d
#     top = q3 + 3 * d
#     return [x for x in coverage if x < bottom or x > top]
#
# print original_implementation(example)
##########################################
################# RORY ###################
##########################################


def _get_depth_for_each_variant(cnf, var_by_site, clipped_gz_vcf_fpath, bed_fpath, bam_fpath):
    # http://www.1000genomes.org/faq/what-depth-coverage-your-phase1-variants
    # bedtools intersect -a oncomine.vcf -b Exons.az_key.bed -header > oncomine.az_key.vcf
    # /opt/az/local/tabix/tabix-0.2.6/bgzip oncomine.az_key.vcf
    # /opt/az/local/tabix/tabix-0.2.6/tabix -h -p vcf oncomine.az_key.vcf.gz
    # samtools view -b TRF004223.sorted.bam -L Exons.az_key.bed | bedtools genomecov -ibam stdin -bg > coverage.bg
    # bedtools intersect -a oncomine.az_key.vcf.gz -b coverage.bg -wa | cut -f1,2,4,5,8,11,12,13,14 > oncomine.az_key.depth_numbers.vcf

    samtools = get_system_path(cnf, 'samtools')
    bedtools = get_system_path(cnf, 'bedtools')

    info()
    info('Depth of coverage for regions in BED ' + bed_fpath)
    cov_bg = join(cnf.work_dir, 'coverage.bg')
    cmdline = '{samtools} view -L {bed_fpath} -b {bam_fpath} | {bedtools} genomecov -ibam stdin -bg'.format(**locals())
    call(cnf, cmdline, output_fpath=cov_bg)

    info()
    info('Intersecting depth regions with VCF ' + clipped_gz_vcf_fpath)
    vcf_depth_numbers_fpath = join(cnf.work_dir, 'vcf_bg.intersect')
    cmdline = '{bedtools} intersect -a {clipped_gz_vcf_fpath} -b {cov_bg} -wao'.format(**locals())
    res = call(cnf, cmdline, output_fpath=vcf_depth_numbers_fpath)
    # if res != oncomine_depth_numbers_fpath:
    #     info()
    #     info('Trying with uncompressed VCF')
    #     cmdline = 'gunzip {vcf_fpath} -c | {bedtools} intersect -a - -b {cov_bg} -wao | cut -f1,2,4,5,8,11,12,13,14,15'.format(**locals())
    #     call(cnf, cmdline, output_fpath=oncomine_depth_numbers_fpath)

    depths_per_var = defaultdict(list)
    with open(vcf_depth_numbers_fpath) as f:
        for l in f:
            # 1,2,4,5,8,11,12,13,14,15,16,17,18,19,20
            # c,p,r,a,f,ch,st,en,ge,ex,st,ft,bt,de,ov
            fs = l[:-1].split('\t')
            chrom, pos, _, ref, alt = fs[:5]
            depth, overlap = fs[-2:]
            var = var_by_site.get((chrom, pos, ref, alt))
            if var and depth != '.':
                depth, overlap = int(depth), int(overlap)
                for i in range(overlap):
                    depths_per_var[(chrom, pos, ref, alt)].append(depth)

    # Getting avarage depth of coverage of each variant (exactly for those parts that were in BED)
    depth_by_var = {var: (sum(depths) / len(depths)) if len(depths) != 0 else None
                    for var, depths in depths_per_var.iteritems()}

    return depth_by_var


def _read_vcf_records_per_bed_region_and_clip_vcf(cnf, vcf_fpath, bed_fpath):
    info()
    info('Intersecting VCF ' + vcf_fpath + ' using BED ' + bed_fpath)

    vcf_columns_num = count_bed_cols(vcf_fpath)
    bed_columns_num = count_bed_cols(bed_fpath)

    vcf_bed_intersect = join(cnf.work_dir, splitext(basename(vcf_fpath))[0] + '_vcf_bed.intersect')
    bedtools = get_system_path(cnf, 'bedtools')
    cmdline = '{bedtools} intersect -header -a {vcf_fpath} -b {bed_fpath} -wo'.format(**locals())
    res = call(cnf, cmdline, output_fpath=vcf_bed_intersect)

    regions_in_order = []
    regions_set = set()
    vars_by_region = defaultdict(dict)
    var_by_site = dict()

    clipped_vcf_fpath = intermediate_fname(cnf, vcf_fpath, 'clip')

    with open(vcf_bed_intersect) as f, open(clipped_vcf_fpath, 'w') as clip_vcf:
        for l in f:
            l = l.strip()
            if not l or l.startswith('#'):
                clip_vcf.write(l + '\n')
                continue
            fs = l.split('\t')
            chrom, pos, id_, ref, alt, qual, filt, info_fields = fs[:8]
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

            cls = None
            if '=Hotspot' in info_fields: cls = 'Hotspot'
            if '=Deleterious' in info_fields: cls = 'Deleterious'
            if cls:
                var = Variant(chrom, pos, ref, alt, cls)
                vars_by_region[r][(chrom, pos, ref, alt)] = var
                var_by_site[(chrom, pos, ref, alt)] = var
                clip_vcf.write('\t'.join([chrom, pos, id_, ref, alt, qual, filt, info_fields]) + '\n')

    clipped_gz_vcf_fpath = bgzip_and_tabix(cnf, clipped_vcf_fpath)

    return clipped_gz_vcf_fpath, regions_in_order, vars_by_region, var_by_site


class Variant(object):
    def __init__(self, chrom, pos, ref, alt, cls):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.cls = cls

    def get_site(self):
        return self.chrom, self.pos, self.ref, self.alt

    def __repr__(self):
        fmt_pos = lambda pos: Metric.format_value(pos, human_readable=True)
        return '{pos} {var.ref}/{var.alt} {var.cls[0]}'.format(pos=fmt_pos(int(self.pos)), var=self)


class VariantsMetric(Metric):
    def format(self, value, human_readable=True):
        variants = value
        fmt_pos = lambda pos: Metric.format_value(int(pos), human_readable=True)
        return '  '.join('{pos}:{var.ref}/{var.alt}'.format(pos=fmt_pos(var.pos), var=var)
            for var in variants)


class DepthsMetric(Metric):
    def format(self, value, human_readable=True):
        depth_tuples = value
        fmt = lambda dp: Metric.format_value(dp, human_readable=True)
        return '  '.join('{depth}/{norm_depth}'.format(
            depth=fmt(int(depth) if depth is not None else None),
            norm_depth=fmt(float(norm_depth) if norm_depth is not None else None))
            for (depth, norm_depth) in depth_tuples)


shared_metrics = [
    Metric('Chr'),
    Metric('Start'),
    Metric('End'),
    Metric('Strand'),
    Metric('Feature'),
    Metric('Biotype'),
    Metric('Symbol'),
    Metric('Hotspots num', short_name='#HS'),
    VariantsMetric('Hotspots list', short_name='Hotspots')]

shared_general_metrics = [Metric('Sample', short_name='Sample', common=True)]

single_report_metric_storage = MetricStorage(
    general_section=ReportSection('general_section', '', metrics=shared_general_metrics + [
        Metric('Average sample depth', short_name='Ave depth', common=True)
    ]),
    sections=[ReportSection(metrics=shared_metrics + [
        DepthsMetric('Hotspots depths/norm depths', short_name='DP/Norm_DP')
    ])])


def _report_normalize_coverage_for_variant_sites(cnf, output_dir, sample, ave_sample_depth, vcf_key,
                                                 bed_fpath, depth_cutoff):
    step_greetings('Normalized coverage for ' + vcf_key + ' hotspots')
    vcf_fpath = cnf.genome.get(vcf_key)
    if not vcf_fpath:
        err('Error: no ' + vcf_key + ' for ' + cnf.genome.name + ' VCF fpath specified in ' + cnf.sys_cnf)
        return None

    clipped_gz_vcf_fpath, regions_in_order, vars_by_region, var_by_site = \
        _read_vcf_records_per_bed_region_and_clip_vcf(cnf, vcf_fpath, bed_fpath)

    depth_by_var = _get_depth_for_each_variant(cnf, var_by_site, clipped_gz_vcf_fpath, bed_fpath, sample.bam)

    info()
    info('Saving report for sample ' + sample.name)

    report = PerRegionSampleReport(sample=sample, metric_storage=single_report_metric_storage)
    report.add_record('Sample', sample.name)
    report.add_record('Average sample depth', ave_sample_depth)

    total_variants = 0
    total_variants_below_cutoff = 0
    total_regions = 0
    for r in regions_in_order:
        total_regions += 1

        (chrom, start, end, symbol, strand, feature, biotype) = r
        rep_region = sample.report.add_row()
        rep_region.add_record('Chr', chrom)
        rep_region.add_record('Start', start)
        rep_region.add_record('End', end)
        rep_region.add_record('Symbol', symbol)
        rep_region.add_record('Strand', strand)
        rep_region.add_record('Feature', feature)
        rep_region.add_record('Biotype', biotype)

        variants, depths = [], []
        for var in sorted(vars_by_region[r].values(), key=lambda v: v.pos):
            total_variants += 1

            depth = depth_by_var.get((var.get_site()))
            if depth_cutoff >= depth_cutoff:
                continue

            total_variants_below_cutoff += 1
            variants.append(var)

            norm_depth = depth / ave_sample_depth if depth and ave_sample_depth > 0 else None
            depths.append((depth, norm_depth))

            total_variants += 1
            if total_variants % 10000 == 0:
                info('Processed {0:,} variants, {0:,} regions.'.format(total_variants, total_regions))

        rep_region.add_record('Total hotspots num in regions', len(variants))
        rep_region.add_record('Hotspots num with depth less than ' + total_variants_below_cutoff)
        rep_region.add_record('Hotspots list', variants)
        rep_region.add_record('Hotspots depths/norm depths', depths)

    best_report_basename = sample.name + '.' + source.targetseq_name  + '_' + vcf_key
    report.save_txt(join(sample.dirpath, source.targetseq_name, best_report_basename + '.txt'))
    report.save_tsv(join(sample.dirpath, source.targetseq_name, best_report_basename + '.tsv'))
    info('')
    info('Oncomine variants coverage report (total: {0:,} variants, {0:,} variants below cut-off {0}, {0:,} regions) '
         'saved into:'.format(total_variants, total_variants_below_cutoff, depth_cutoff, total_regions))
    info('  ' + report.txt_fpath)

    return report


# def make_flagged_regions_reports(cnf, targetseq_dir, sample, filtered_vcf_by_callername=None):
#     if not filtered_vcf_by_callername:
#         info('No variant callset')
#         # return []
#
#     detail_gene_rep_fpath = join(
#         cnf.output_dir,
#         sample.name + '.' + targetseq_dir + source.detail_gene_report_baseending + '.tsv')
#
#     info('Reading regions from ' + detail_gene_rep_fpath)
#     regions = _read_regions(detail_gene_rep_fpath)
#
#     abnormal_regions_reports = _generate_flagged_regions_reports(
#         cnf, sample, regions, filtered_vcf_by_callername or dict())
#
#     return abnormal_regions_reports
#
#
# def _read_regions(gene_report_fpath):
#     regions = []
#     with open(gene_report_fpath) as f:
#         for line in f:
#             tokens = line.strip().split('\t')
#             if line.startswith('Sample') or line.startswith('#'):
#                 depth_threshs = [int(v[:-1]) for v in tokens[12:]]
#                 continue
#
#             sample_name, chrom, start, end, gene, exon_num, strand, feature, size, avg_depth, \
#                 std_dev, percent_within_normal = tokens[:12]
#             region = Region(
#                 sample_name=sample_name,
#                 chrom=chrom,
#                 start=int(''.join([c for c in start if c.isdigit()])),
#                 end=int(''.join([c for c in end if c.isdigit()])),
#                 gene_name=gene,
#                 exon_num=int(exon_num) if exon_num else None,
#                 strand=strand or None,
#                 feature=feature,
#                 size=int(''.join([c for c in size if c.isdigit()])),
#                 avg_depth=float(avg_depth),
#                 std_dev=float(std_dev),
#                 rate_within_normal=float(percent_within_normal[:-1])
#             )
#
#             percents = []
#             for token in tokens[12:]:
#                 try:
#                     v = float(token[:-1])
#                 except:
#                     v = '-'
#                 percents.append(v)
#
#             for thresh, percent in zip(depth_threshs, percents):
#                 region.percent_within_threshs[thresh] = percent
#             regions.append(region)
#     return regions
#
#
# def classify_based_on_min_and_max(cnf, sample, regions):
#     def simple_min_and_max():
#         cov_cnf = cnf.coverage_reports
#         _min_cov = min(cov_cnf.min_cov_factor * sample.median_cov, cov_cnf.min_cov)
#         _max_cov = cov_cnf.max_cov_factor * sample.median_cov
#
#         info('Assuming abnormal if below min(median*' +
#              str(cov_cnf.min_cov_factor) + ', ' + str(cov_cnf.min_cov) + ') = ' +
#              str(_min_cov) + ', or above median*' + str(cov_cnf.max_cov_factor) +
#              ' = ' + str(_max_cov))
#         return _min_cov, _max_cov
#
#     def min_and_max_based_on_outliers():
#         rs_by_depth = sorted(regions, key=lambda r: r.avg_depth)
#         l = len(rs_by_depth)
#         q1 = rs_by_depth[int((l - 1) / 4)].avg_depth
#         q3 = rs_by_depth[int((l - 1) * 3 / 4)].avg_depth
#         d = q3 - q1
#         _min_cov = q1 - 3 * d
#         _max_cov = q1 + 3 * d
#
#         info('Using outliers mechanism. l = {}, q1 = {}, q3 = {}, d = {},'
#              'min_cov = {}, max_cov = {}'.format(l, q1, q3, d, _min_cov, _max_cov))
#         return _min_cov, _max_cov
#
#     min_cov, max_cov = min_and_max_based_on_outliers()
#
#     low_regions, high_regions = [], []
#
#     def _classify_region(region, median_cov, minimal_cov, maximal_cov):
#         if region.avg_depth < minimal_cov:
#             region.cov_factor = region.avg_depth / median_cov
#             low_regions.append(region)
#
#         if region.avg_depth > maximal_cov:
#             region.cov_factor = region.avg_depth / median_cov
#             high_regions.append(region)
#
#     info('Classifying regions...')
#     proc_regions(regions, _classify_region, low_regions, high_regions,
#         sample.median_cov, min_cov, max_cov)
#
#     return low_regions, high_regions
#
#
# def classify_based_on_justin(cnf, sample, regions):
#     # Экзоны и ампликоны по отдельности. Если среднее покрытие 12k,
#     # выбрать колонку где покрытие хотя бы на 2k, отсортировать и показать
#     # те строки, где меньше чем на 80% покрытие. Вывести cosmic id где в этих
#     # регионах. Вывести HTML репорты со ссылками на бамы и vcf в IGV.
#     #
#     # Саммари репорт должен содержать те экзоны, которые оказались флаггед
#     # хотя бы в 20% проектов.
#     low_regions, high_regions = [], []
#
#
#     return low_regions, high_regions
#
#
# def _generate_flagged_regions_reports(cnf, sample, regions, filtered_vcf_by_callername):
#     info('Calculation median coverage...')
#
#     median_cov = median((r.avg_depth for r in regions))
#     sample.median_cov = median_cov
#     if median_cov == 0:
#         err('Median coverage is 0')
#         return None, None
#     info('Median: ' + str(median_cov))
#     info()
#
#     info('Extracting abnormally covered regions.')
#     # low_regions, high_regions = classify_based_on_min_and_max(cnf, sample, regions)
#     low_regions, high_regions = classify_based_on_justin(cnf, sample, regions)
#
#     abnormal_regions_report_fpaths = []
#     for caller_name, vcf_fpath in filtered_vcf_by_callername:
#         vcf_dbs = [
#             # VCFDataBase('dbsnp', 'DBSNP', cnf.genome),
#             VCFDataBase('cosmic', 'Cosmic', cnf.genome),
#             VCFDataBase('oncomine', 'Oncomine', cnf.genome),
#         ]
#         for kind, regions, f_basename in zip(
#                 ['low', 'high'],
#                 [low_regions, high_regions],
#                 [BCBioStructure.detail_lowcov_gene_report_baseending,
#                  BCBioStructure.detail_highcov_gene_report_baseending]):
#             if len(regions) == 0:
#                 err('No flagged regions with too ' + kind + ' coverage, skipping counting missed variants.')
#                 continue
#
#             report_base_name = cnf.name + '_' + caller_name + '.' + BCBioStructure.targetseq_name + f_basename
#
#             abnormal_regions_bed_fpath = save_regions_to_bed(cnf, regions, report_base_name)
#             for vcf_db in vcf_dbs:
#                 vcf_db.missed_vcf_fpath, vcf_db.annotated_bed = _find_missed_variants(cnf, vcf_db,
#                     caller_name, vcf_fpath, report_base_name, abnormal_regions_bed_fpath)
#
#             report = _make_flagged_region_report(cnf, sample, regions, vcf_fpath, caller_name, vcf_dbs)
#
#             regions_html_rep_fpath = report.save_html(cnf.output_dir, report_base_name,
#                 caption='Regions with ' + kind + ' coverage for ' + caller_name)
#
#             regions_txt_rep_fpath = report.save_txt(cnf.output_dir, report_base_name,
#                 [report.metric_storage.sections_by_name[kind + '_cov']])
#
#             abnormal_regions_report_fpaths.append(regions_txt_rep_fpath)
#             info('Too ' + kind + ' covered regions (total ' + str(len(regions)) + ') saved into:')
#             info('  HTML: ' + str(regions_html_rep_fpath))
#             info('  TXT:  ' + regions_txt_rep_fpath)
#
#     return abnormal_regions_report_fpaths
#
#
# class VCFDataBase():
#     def __init__(self, name, descriptive_name, genome_cnf):
#         self.name = name
#         self.descriptive_name = descriptive_name
#         self.vcf_fpath = genome_cnf.get(name)
#         self.annotated_bed = None
#         self.missed_vars_by_region_info = OrderedDict()
#
#
# def _make_flagged_region_report(cnf, sample, regions, filtered_vcf_fpath, caller_name, vcf_dbs):
#     regions_metrics = [
#         Metric('Sample'),
#         Metric('Chr'),
#         Metric('Start'),
#         Metric('End'),
#         Metric('Gene'),
#         Metric('Feature'),
#         Metric('Size'),
#         Metric('AvgDepth'),
#         Metric('StdDev', description='Coverage depth standard deviation'),
#         Metric('Wn 20% of Mean', unit='%', description='Persentage of the region that lies within 20% of an avarage depth.'),
#         Metric('Depth/Median', description='Average depth of coveraege of the region devided by median coverage of the sample'),
#         Metric('Total variants'),
#     ]
#     for vcf_db in vcf_dbs:
#         regions_metrics.append(Metric(vcf_db.descriptive_name + ' missed'))
#
#     for thres in cnf.coverage_reports.depth_thresholds:
#         regions_metrics.append(Metric('{}x'.format(thres), unit='%',
#             description='Bases covered by at least {} reads'.format(thres)))
#
#     general_metrics = [
#         Metric('Median depth'),
#         Metric('Variants'),
#     ]
#     for vcf_db in vcf_dbs:
#         general_metrics.append(
#             Metric(vcf_db.descriptive_name + ' missed variants'))
#
#     region_metric_storage = MetricStorage(
#         general_section=ReportSection(metrics=general_metrics),
#         sections_by_name=OrderedDict(
#             low_cov=ReportSection('Low covered regions', 'Low covered regions', regions_metrics[:]),
#             high_cov=ReportSection('Regions', 'Amplcons and exons coverage depth statistics', regions_metrics[:],
#         )))
#
#     report = PerRegionSampleReport(sample, metric_storage=region_metric_storage)
#
#     for vcf_db in vcf_dbs:
#         if not vcf_db.annotated_bed:
#             report.add_record(vcf_db.descriptive_name + ' missed variants', None)
#             continue
#
#         with open(vcf_db.annotated_bed) as f:
#             for line in f:
#                 if not line.startswith('#'):
#                     tokens = line.strip().split('\t')
#                     chrom, start, end, gene, feature, missed_vars = tokens[:9]
#                     vcf_db.missed_vars_by_region_info[(chrom, int(start), int(end), feature)] = int(missed_vars)
#
#         info('Missed variants from ' + vcf_db.descriptive_name + ': ' +
#              str(sum(vcf_db.missed_vars_by_region_info.values())))
#         report.add_record(vcf_db.descriptive_name + ' missed variants', vcf_db.missed_vcf_fpath)
#
#         info()
#
#     # cosmic_missed_vcf_fpath = add_suffix(filtered_vcf_fpath, 'cosmic')
#     # shutil.copy(filtered_vcf_fpath, cosmic_missed_vcf_fpath)
#
#     report.add_record('Median depth', sample.median_cov)
#     report.add_record('Variants', filtered_vcf_fpath)
#
#     for r in regions:
#         for vcf_db in vcf_dbs:
#             if vcf_db.descriptive_name not in r.missed_by_db:
#                 r.missed_by_db[vcf_db.descriptive_name] = 0
#             r.missed_by_db[vcf_db.descriptive_name] += \
#                 vcf_db.missed_vars_by_region_info.get((r.chrom, r.start, r.end, r.feature)) or 0
#
#     rec_by_name = {metric.name: Record(metric, []) for metric in regions_metrics}
#     for rec in rec_by_name.values():
#         report.records.append(rec)
#
#     proc_regions(regions, _fill_in_record_info, rec_by_name, sample.median_cov)
#
#     return report
#
#
# def _add_vcf_header(source_vcf_fpath, dest_vcf_fpath):
#     with open(source_vcf_fpath) as inp, open(dest_vcf_fpath, 'w') as out:
#         for line in inp:
#             if line.startswith('#'):
#                 out.write(line)
#     with open(dest_vcf_fpath + '_tmp') as inp, open(dest_vcf_fpath, 'a') as out:
#         for line in inp:
#             out.write(line)
#     os.remove(dest_vcf_fpath + '_tmp')
#     info(dest_vcf_fpath)
#
#
# def _find_missed_variants(cnf, vcf_db, caller_name, input_vcf_fpath, report_base_name, bed_fpath):
#     if not vcf_db.vcf_fpath and not verify_file(vcf_db.vcf_fpath, vcf_db.descriptive_name):
#         info('Skipping counting variants from ' + vcf_db.descriptive_name)
#         return None, None
#
#     info('Counting missed variants from ' + vcf_db.descriptive_name)
#     db_in_roi_vcf_fpath = join(cnf.work_dir, report_base_name + '_' + vcf_db.name + '_roi.vcf')
#     missed_vcf_fpath = join(cnf.output_dir, report_base_name + '_' + caller_name + '_' + vcf_db.name + '_missed.vcf')
#     annotated_regions_bed_fpath = join(cnf.work_dir, report_base_name + '_' + caller_name + '_' + vcf_db.name + '.bed')
#
#     bedtools = get_system_path(cnf, 'bedtools')
#
#     db_vcf_fpath = vcf_db.vcf_fpath
#     # Take variants from DB VCF that lie in BED
#     # (-wa means to take whole regions from a, not the parts that overlap)
#     cmdline = '{bedtools} intersect -wa -a {db_vcf_fpath} -b {bed_fpath}'.format(**locals())
#     call(cnf, cmdline, output_fpath=db_in_roi_vcf_fpath + '_tmp')
#     _add_vcf_header(db_vcf_fpath, db_in_roi_vcf_fpath)
#     info()
#
#     # Take variants from VCF missed in DB (but only in the regions from BED)
#     cmdline = '{bedtools} intersect -v -a {input_vcf_fpath} -b {db_in_roi_vcf_fpath}'.format(**locals())
#     call(cnf, cmdline, output_fpath=missed_vcf_fpath + '_tmp')
#     _add_vcf_header(input_vcf_fpath, missed_vcf_fpath)
#     info()
#
#     # Take variants from VCF missed in DB (but only in the regions from BED)
#     cmdline = '{bedtools} intersect -c -a {bed_fpath} -b {missed_vcf_fpath}'.format(**locals())
#     call(cnf, cmdline, output_fpath=annotated_regions_bed_fpath)
#     info()
#
#     return missed_vcf_fpath, annotated_regions_bed_fpath
#
#
# # def _count_variants(region, vcf_fpath):
# #     num = 0
# #     with open(vcf_fpath) as inp:
# #         reader = vcf_parser.Reader(inp)
# #         for rec in reader:
# #             if region.start <= rec.POS <= region.end:
# #                 num += 1
# #     return num
#
#
# def _fill_in_record_info(region, rec_by_name, median_depth):
#     rec_by_name['Sample'].value.append(         region.sample_name)
#     rec_by_name['Chr'].value.append(            region.chrom)
#     rec_by_name['Start'].value.append(          region.start)
#     rec_by_name['End'].value.append(            region.end)
#     rec_by_name['Gene'].value.append(           region.gene_name)
#     rec_by_name['Feature'].value.append(        region.feature)
#     rec_by_name['Size'].value.append(           region.get_size())
#     rec_by_name['AvgDepth'].value.append(      region.avg_depth)
#     rec_by_name['StdDev'].value.append(         region.std_dev)
#     rec_by_name['Wn 20% of Mean'].value.append( region.percent_within_normal)
#     rec_by_name['Depth/Median'].value.append(   region.avg_depth / median_depth if median_depth != 0 else None)
#     rec_by_name['Total variants'].value.append( region.var_num)
#
#     for db_descriptive_name, num in region.missed_by_db.items():
#         rec_by_name[db_descriptive_name + ' missed'].value.append(num)
#
#     for depth_thres, percent in region.percent_within_threshs.items():
#         rec_by_name['{}x'.format(depth_thres)].value.append(percent)

