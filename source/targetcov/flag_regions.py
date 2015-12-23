# coding=utf-8

from collections import defaultdict
import os
from os.path import join, splitext, basename, dirname, abspath

from source.calling_process import call
from source.file_utils import intermediate_fname, verify_file, add_suffix, safe_mkdir
from source.logger import info, err, step_greetings
from source.reporting.reporting import Metric, MetricStorage, ReportSection, PerRegionSampleReport, load_records
import source
from source.targetcov.Region import save_regions_to_bed
from source.targetcov.bam_and_bed_utils import count_bed_cols
from source.targetcov.cov import get_detailed_metric_storage, make_flat_region_report
from source.tools_from_cnf import get_system_path
from source.variants.vcf_processing import bgzip_and_tabix

DEPTH_THRESH_FROM_AVE_COV = 0.8
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
    safe_mkdir(sample.flagged_regions_dirpath)
    ''' 1. Detect depth threshold (ave sample coverage * DEPTH_THRESH_FROM_AVE_COV)
        2. Select regions covered in less than MIN_DEPTH_PERCENT_AT_THRESH at threshold
        3. Sort by % at threshold
        4. Select those parts of those regions where % = 0, save to BED
        5. Find HotSpots at those regions
        6. Intersect HotSpots with tracks

        For each gene where are regions with parts % = 0:
            sort them by part where % = 0
    '''
    #vcf_dbs = ['oncomine', 'dbsnp', 'cosmic']
    vcf_dbs = ['oncomine']

    from source.clinical_reporting.clinical_parser import get_key_or_target_bed_genes
    key_genes, _ = get_key_or_target_bed_genes(cnf.bed, cnf.key_genes)
    depth_cutoff = get_depth_cutoff(ave_depth, depth_threshs)
    genes_sorted = sorted(gene_by_name.values())

    for coverage_type in ['low', 'high']:
        info('Selecting and saving ' + coverage_type + ' covered genes')
        if coverage_type == 'low':
            selected_genes = [g for g in genes_sorted if g.gene_name in key_genes and (
                any(e.rates_within_threshs[depth_cutoff] < MIN_DEPTH_PERCENT_AT_THRESH for e in g.get_exons()) or
                any(a.rates_within_threshs[depth_cutoff] < MIN_DEPTH_PERCENT_AT_THRESH for a in g.get_amplicons()))]
        else:
            min_cov, max_cov = min_and_max_based_on_outliers(genes_sorted)
            selected_genes = [g for g in genes_sorted if g.gene_name in key_genes and (
                any(e.avg_depth > max_cov for e in g.get_exons()) or
                any(a.avg_depth > max_cov for a in g.get_amplicons()))]
        for region_type in ['exons', 'amplicons']:
            selected_regions = []
            for gene in selected_genes:
                if coverage_type == 'low':
                    cur_regions = [a for a in (gene.get_amplicons() if region_type == 'amplicons' else gene.get_exons())
                                   if a.rates_within_threshs[depth_cutoff] < MIN_DEPTH_PERCENT_AT_THRESH and 'Multi' not in a.feature]
                else:
                    cur_regions = [a for a in (gene.get_amplicons() if region_type == 'amplicons' else gene.get_exons())
                                   if a.avg_depth > max_cov and 'Multi' not in a.feature]
                selected_regions.extend(cur_regions)

            if selected_regions:
                selected_regions_bed_fpath = join(sample.flagged_regions_dirpath, coverage_type + '_cov_' + region_type + '.bed')
                save_regions_to_bed(cnf, selected_regions, selected_regions_bed_fpath)

                # Report cov for Hotspots
                for db in vcf_dbs:
                    _report_normalize_coverage_for_variant_sites(cnf, sample, ave_depth, db, selected_regions_bed_fpath, selected_regions,
                                                                 depth_cutoff, region_type, coverage_type)

            report = make_flat_region_report(sample, selected_regions, depth_threshs)
            flagged_txt_fpath = add_suffix(add_suffix(sample.flagged_txt, region_type), coverage_type)
            flagged_tsv_fpath = add_suffix(add_suffix(sample.flagged_tsv, region_type), coverage_type)
            report.save_txt(flagged_txt_fpath)
            report.save_tsv(flagged_tsv_fpath)

            info('')
            info(coverage_type + ' covered ' + region_type + '(total ' + str(len(selected_regions)) + ') for sample ' + sample.name + ' saved into:')
            info('  ' + flagged_txt_fpath + ', ' + flagged_tsv_fpath)

    return report

def min_and_max_based_on_outliers(regions):
    rs_by_depth = sorted(regions, key=lambda r: r.avg_depth)
    l = len(rs_by_depth)
    q1 = rs_by_depth[int((l - 1) / 4)].avg_depth
    q3 = rs_by_depth[int((l - 1) * 3 / 4)].avg_depth
    d = q3 - q1
    _min_cov = max(0, q1 - 3 * d)
    _max_cov = q1 + 3 * d

    info('Using outliers mechanism. l = {}, q1 = {}, q3 = {}, d = {},'
         'min_cov = {}, max_cov = {}'.format(l, q1, q3, d, _min_cov, _max_cov))
    return _min_cov, _max_cov

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

    sambamba = get_system_path(cnf, 'sambamba')
    bedtools = get_system_path(cnf, 'bedtools')

    info()
    info('Depth of coverage for regions in BED ' + bed_fpath)
    cov_bg = join(cnf.work_dir, 'coverage.bg')
    cmdline = '{sambamba} view -f bam -t {cnf.threads} -L {bed_fpath} {bam_fpath} | {bedtools} genomecov -ibam stdin -bg'.format(**locals())
    call(cnf, cmdline, output_fpath=cov_bg, exit_on_error=False)

    info()
    info('Intersecting depth regions with VCF ' + clipped_gz_vcf_fpath)
    vcf_depth_numbers_fpath = join(cnf.work_dir, 'vcf_bg.intersect')
    if not cnf.reuse_intermediate or not verify_file(vcf_depth_numbers_fpath, silent=True, is_critical=False):
        cmdline = '{bedtools} intersect -a {clipped_gz_vcf_fpath} -b {cov_bg} -wao'.format(**locals())
        res = call(cnf, cmdline, output_fpath=vcf_depth_numbers_fpath, exit_on_error=False)
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

    # Getting average depth of coverage of each variant (exactly for those parts that were in BED)
    depth_by_var = {var: (sum(depths) / len(depths)) if len(depths) != 0 else None
                    for var, depths in depths_per_var.iteritems()}

    return depth_by_var


def _read_vcf_records_per_bed_region_and_clip_vcf(cnf, vcf_fpath, bed_fpath, region_type, sample):
    info()
    info('Intersecting VCF ' + vcf_fpath + ' using BED ' + bed_fpath)

    vcf_columns_num = count_bed_cols(vcf_fpath)
    bed_columns_num = count_bed_cols(bed_fpath)

    vcf_bed_intersect = join(cnf.work_dir, splitext(basename(vcf_fpath))[0] + '_' + region_type + '_vcf_bed.intersect')
    bedtools = get_system_path(cnf, 'bedtools')
    if not cnf.reuse_intermediate or not verify_file(vcf_bed_intersect, silent=True, is_critical=False):
        cmdline = '{bedtools} intersect -header -a {vcf_fpath} -b {bed_fpath} -wo'.format(**locals())
        res = call(cnf, cmdline, output_fpath=vcf_bed_intersect)

    regions_in_order = []
    regions_set = set()
    vars_by_region = defaultdict(dict)
    var_by_site = dict()

    clipped_vcf_fpath = intermediate_fname(cnf, vcf_fpath, '_' + region_type + '_clip')

    with open(vcf_bed_intersect) as f, open(clipped_vcf_fpath, 'w') as clip_vcf:
        for l in f:
            l = l.strip()
            if not l or l.startswith('#'):
                clip_vcf.write(l + '\n')
                continue
            fs = l.split('\t')
            chrom, pos, id_, ref, alt, qual, filt, info_fields = fs[:8]
            chrom_b, start_b, end_b, symbol, strand, feature, biotype = None, None, None, None, None, None, None
            if bed_columns_num >= 8:
                chrom_b, start_b, end_b, symbol, _, strand, feature, biotype, _ = fs[-(bed_columns_num+1):][:9]
            elif bed_columns_num >= 4:
                chrom_b, start_b, end_b, symbol, _ = fs[-(bed_columns_num+1):][:5]
            assert chrom == chrom_b, l
            r = chrom, id_, start_b, end_b, symbol, strand, feature, biotype
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


class LowCoveredRegion(object):
    def __init__(self, gene, chrom, ave_depth, exon_amplicone=None, start=None, end=None, mut=None, sample=None):
        self.gene = gene
        self.chrom = chrom
        self.ave_depth = ave_depth
        self.exon_amplicone = exon_amplicone
        self.start = start
        self.end = end
        self.mut = mut
        self.sample = sample


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
    Metric('Gene'),
    Metric('Chr'),
    Metric('Start'),
    Metric('End'),
    Metric('Strand'),
    Metric('Feature'),
    Metric('Biotype'),
    Metric('ID'),
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


def _intersect_with_tricky_regions(cnf, selected_bed_fpath, sample):
    info()
    info('Detecting problematic regions for ' + sample)
    tricky_regions = {'low_gc.bed': 'Low GC', 'high_gc.bed': 'High GC', 'low_complexity.bed': 'Low complexity',
                      'bad_promoter.bed': 'Bad promoter'}
    bed_filenames = tricky_regions.keys()

    bed_fpaths = [join(cnf.tricky_regions, bed_filename) for bed_filename in bed_filenames]

    info('Intersecting BED ' + selected_bed_fpath + ' using BED files with tricky regions')

    vcf_bed_intersect = join(cnf.work_dir, splitext(basename(selected_bed_fpath))[0] + '_tricky_vcf_bed.intersect')
    if not cnf.reuse_intermediate or not verify_file(vcf_bed_intersect, silent=True, is_critical=False):
        bedtools = get_system_path(cnf, 'bedtools')
        cmdline = bedtools + ' intersect -header -a ' + selected_bed_fpath + ' -b ' + ' '.join(bed_fpaths) + ' -wa -wb -filenames'
        call(cnf, cmdline, output_fpath=vcf_bed_intersect, exit_on_error=False)

    regions_by_reasons = {}

    with open(vcf_bed_intersect) as f:
        for l in f:
            l = l.strip()
            if not l or l.startswith('#'):
                continue
            fs = l.split('\t')
            pos = fs[1]
            filename = fs[-4]
            regions_by_reasons.setdefault(pos, set()).add(tricky_regions[os.path.basename(filename)])

    return regions_by_reasons

def _report_normalize_coverage_for_variant_sites(cnf, sample, ave_sample_depth, vcf_key, bed_fpath, selected_regions,
                                                 depth_cutoff, region_type, coverage_type):
    info()
    info('Normalized coverage for ' + vcf_key + ' hotspots (' + region_type + ')')
    vcf_fpath = cnf.genome.get(vcf_key)
    if not vcf_fpath:
        err('Error: no ' + vcf_key + ' for ' + cnf.genome.name + ' VCF fpath specified in ' + cnf.sys_cnf)
        return None

    clipped_gz_vcf_fpath, regions_in_order, vars_by_region, var_by_site = \
        _read_vcf_records_per_bed_region_and_clip_vcf(cnf, vcf_fpath, bed_fpath, region_type, sample)

    depth_by_var = _get_depth_for_each_variant(cnf, var_by_site, clipped_gz_vcf_fpath, bed_fpath, sample.bam)

    info()
    info('Saving report for ' + coverage_type + ' covered ' + region_type + ' for sample ' + sample.name)
    report = PerRegionSampleReport(sample=sample, metric_storage=single_report_metric_storage)
    report.add_record('Sample', sample.name)
    report.add_record('Average sample depth', ave_sample_depth)

    total_variants = 0
    total_variants_below_cutoff = 0
    total_regions = 0
    for r in regions_in_order:
        total_regions += 1

        (chrom, id_, start, end, gene, strand, feature, biotype) = r
        for region in selected_regions:
            if region.start == int(start) and region.end == int(end):
                strand = region.strand
                feature = region.feature
                biotype = region.biotype
                break

        rep_region = report.add_row()
        rep_region.add_record('Gene', gene)
        rep_region.add_record('Chr', chrom)
        rep_region.add_record('Start', start)
        rep_region.add_record('End', end)
        rep_region.add_record('Strand', strand)
        rep_region.add_record('Feature', feature)
        rep_region.add_record('Biotype', biotype)
        rep_region.add_record('ID', id_)
        variants, depths = [], []
        for var in sorted(vars_by_region[r].values(), key=lambda v: v.pos):
            total_variants += 1

            depth = depth_by_var.get((var.get_site()))
            if depth >= depth_cutoff:
                continue

            total_variants_below_cutoff += 1
            variants.append(var)

            norm_depth = depth / ave_sample_depth if depth and ave_sample_depth > 0 else None
            depths.append((depth, norm_depth))

            total_variants += 1
            if total_variants % 10000 == 0:
                info('Processed {0:,} variants, {0:,} regions.'.format(total_variants, total_regions))

        rep_region.add_record('Hotspots num', len(variants))
        rep_region.add_record('Hotspots list', variants)
        rep_region.add_record('Hotspots depths/norm depths', depths)

    best_report_basename = coverage_type + '_cov_' + region_type + '.' + vcf_key
    report.save_txt(join(sample.flagged_regions_dirpath, best_report_basename + '.txt'))
    report.save_tsv(join(sample.flagged_regions_dirpath, best_report_basename + '.tsv'))
    info('')
    info(vcf_key + ' variants coverage report (total: {0:,} variants, {1:,} variants below cut-off {2}, {3:,} regions) '
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

