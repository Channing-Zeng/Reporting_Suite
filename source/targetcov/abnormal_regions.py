# coding=utf-8
from collections import OrderedDict
import os
from os.path import join
import shutil
from ext_modules import vcf_parser
from source.bcbio_structure import BCBioStructure
from source.calling_process import call
from source.file_utils import add_suffix, verify_file

from source.logger import info, err
from source.reporting import Metric, MetricStorage, ReportSection, Record, SquareSampleReport
from source.targetcov.Region import Region, _proc_regions, _save_regions_to_bed
from source.tools_from_cnf import get_tool_cmdline
from source.utils import median


def make_abnormal_regions_reports(cnf, sample, filtered_vcf_by_callername=None):
    if not filtered_vcf_by_callername:
        info('No variants, skipping flagging regions reports...')
        return []

    detail_gene_rep_fpath = join(
        cnf.output_dir,
        sample.name + '.' + BCBioStructure.targetseq_dir +
        BCBioStructure.detail_gene_report_baseending + '.tsv')

    info('Reading regions from ' + detail_gene_rep_fpath)
    regions = _read_regions(detail_gene_rep_fpath)

    abnormal_regions_reports = _generate_abnormal_regions_reports(
        cnf, sample, regions, filtered_vcf_by_callername or dict())

    return abnormal_regions_reports


def _read_regions(gene_report_fpath):
    regions = []
    with open(gene_report_fpath) as f:
        for line in f:
            tokens = line.strip().split('\t')
            if line.startswith('Sample') or line.startswith('#'):
                depth_threshs = [int(v[:-1]) for v in tokens[12:]]
                continue

            sample_name, chrom, start, end, gene, exon_num, strand, feature, size, avg_depth, \
                std_dev, percent_within_normal = tokens[:12]
            region = Region(
                sample_name=sample_name,
                chrom=chrom,
                start=int(''.join([c for c in start if c.isdigit()])),
                end=int(''.join([c for c in end if c.isdigit()])),
                gene_name=gene,
                exon_num=int(exon_num) if exon_num else None,
                strand=strand or None,
                feature=feature,
                size=int(''.join([c for c in size if c.isdigit()])),
                avg_depth=float(avg_depth),
                std_dev=float(std_dev),
                percent_within_normal=float(percent_within_normal[:-1])
            )

            percents = []
            for token in tokens[12:]:
                try:
                    v = float(token[:-1])
                except:
                    v = '-'
                percents.append(v)

            for thresh, percent in zip(depth_threshs, percents):
                region.percent_within_threshs[thresh] = percent
            regions.append(region)
    return regions


def _generate_abnormal_regions_reports(cnf, sample, regions, filtered_vcf_by_callername):
    info('Calculation median coverage...')
    median_cov = median((r.avg_depth for r in regions))
    sample.median_cov = median_cov
    if median_cov == 0:
        err('Median coverage is 0')
        return None, None

    info('Median: ' + str(median_cov))
    info()

    info('Extracting abnormally covered regions.')
    # minimal_cov = min(cnf.coverage_reports.min_cov_factor * median_cov, cnf.coverage_reports.min_cov)
    # maximal_cov = cnf.coverage_reports.max_cov_factor * median_cov
    #
    # info('Assuming abnormal if below min(median*' +
    #      str(cnf.coverage_reports.min_cov_factor) + ', ' +
    #      str(cnf.coverage_reports.min_cov) + ') = ' + str(minimal_cov) +
    #      ', or above median*' + str(cnf.coverage_reports.max_cov_factor) +
    #      ' = ' + str(maximal_cov))

    # Outliers method
    rs_by_depth = sorted(regions, key=lambda r: r.avg_depth)
    l = len(rs_by_depth)
    q1 = rs_by_depth[int((l - 1) / 4)].avg_depth
    q3 = rs_by_depth[int((l - 1) * 3 / 4)].avg_depth
    d = q3 - q1
    minimal_cov = q1 - 3 * d
    maximal_cov = q1 + 3 * d

    info('Classifying regions...')
    low_regions, high_regions = [], []
    _proc_regions(regions, _classify_region, low_regions, high_regions,
                  median_cov, minimal_cov, maximal_cov)

    abnormal_regions_report_fpaths = []
    for caller_name, vcf_fpath in filtered_vcf_by_callername:
        vcf_dbs = [
            # VCFDataBase('dbsnp', 'DBSNP', cnf.genome),
            VCFDataBase('cosmic', 'Cosmic', cnf.genome),
            VCFDataBase('oncomine', 'Oncomine', cnf.genome),
        ]
        for kind, regions, f_basename in zip(
                ['low', 'high'],
                [low_regions, high_regions],
                [BCBioStructure.detail_lowcov_gene_report_baseending,
                 BCBioStructure.detail_highcov_gene_report_baseending]):
            if len(regions) == 0:
                err('No flagged regions with too ' + kind + ' coverage, skipping counting missed variants.')
                continue

            report_base_name = cnf.name + '_' + caller_name + '.' + BCBioStructure.targetseq_name + f_basename

            abnormal_regions_bed_fpath = _save_regions_to_bed(cnf, regions, report_base_name)
            for vcf_db in vcf_dbs:
                vcf_db.missed_vcf_fpath, vcf_db.annotated_bed = _find_missed_variants(cnf, vcf_db,
                    caller_name, vcf_fpath, report_base_name, abnormal_regions_bed_fpath)

            report = _make_flagged_region_report(cnf, sample, regions, vcf_fpath, caller_name, vcf_dbs)

            regions_html_rep_fpath = report.save_html(cnf.output_dir, report_base_name,
                caption='Regions with ' + kind + ' coverage for ' + caller_name)

            regions_txt_rep_fpath = report.save_txt(cnf.output_dir, report_base_name,
                [report.metric_storage.sections_by_name[kind + '_cov']])

            abnormal_regions_report_fpaths.append(regions_txt_rep_fpath)
            info('Too ' + kind + ' covered regions (total ' + str(len(regions)) + ') saved into:')
            info('  HTML: ' + str(regions_html_rep_fpath))
            info('  TXT:  ' + regions_txt_rep_fpath)

    return abnormal_regions_report_fpaths


def _classify_region(region, low_regions, high_regions, median_cov, minimal_cov, maximal_cov):
    if region.avg_depth < minimal_cov:
        region.cov_factor = region.avg_depth / median_cov
        low_regions.append(region)

    if region.avg_depth > maximal_cov:
        region.cov_factor = region.avg_depth / median_cov
        high_regions.append(region)


class VCFDataBase():
    def __init__(self, name, descriptive_name, genome_cnf):
        self.name = name
        self.descriptive_name = descriptive_name
        self.vcf_fpath = genome_cnf.get(name)
        self.annotated_bed = None
        self.missed_vars_by_region_info = OrderedDict()


def _make_flagged_region_report(cnf, sample, regions, filtered_vcf_fpath, caller_name, vcf_dbs):
    regions_metrics = [
        Metric('Sample'),
        Metric('Chr'),
        Metric('Start'),
        Metric('End'),
        Metric('Gene'),
        Metric('Feature'),
        Metric('Size'),
        Metric('AvgDepth'),
        Metric('StdDev', description='Coverage depth standard deviation'),
        Metric('Wn 20% of Mean', unit='%', description='Persentage of the region that lies within 20% of an avarage depth.'),
        Metric('Depth/Median', description='Average depth of coveraege of the region devided by median coverage of the sample'),
        Metric('Total variants'),
    ]
    for vcf_db in vcf_dbs:
        regions_metrics.append(Metric(vcf_db.descriptive_name + ' missed'))

    for thres in cnf.coverage_reports.depth_thresholds:
        regions_metrics.append(Metric('{}x'.format(thres), unit='%',
            description='Bases covered by at least {} reads'.format(thres)))

    general_metrics = [
        Metric('Median depth'),
        Metric('Variants'),
    ]
    for vcf_db in vcf_dbs:
        general_metrics.append(
            Metric(vcf_db.descriptive_name + ' missed variants'))

    region_metric_storage = MetricStorage(
        general_section=ReportSection(metrics=general_metrics),
        sections_by_name=OrderedDict(
            low_cov=ReportSection('Low covered regions', 'Low covered regions', regions_metrics[:]),
            high_cov=ReportSection('Regions', 'Amplcons and exons coverage depth statistics', regions_metrics[:],
        )))

    report = SquareSampleReport(sample, metric_storage=region_metric_storage)

    for vcf_db in vcf_dbs:
        if not vcf_db.annotated_bed:
            report.add_record(vcf_db.descriptive_name + ' missed variants', None)
            continue

        with open(vcf_db.annotated_bed) as f:
            for line in f:
                if not line.startswith('#'):
                    tokens = line.strip().split('\t')
                    chrom, start, end, gene, feature, missed_vars = tokens[:9]
                    vcf_db.missed_vars_by_region_info[(chrom, int(start), int(end), feature)] = int(missed_vars)

        info('Missed variants from ' + vcf_db.descriptive_name + ': ' +
             str(sum(vcf_db.missed_vars_by_region_info.values())))
        report.add_record(vcf_db.descriptive_name + ' missed variants', vcf_db.missed_vcf_fpath)

        info()

    # cosmic_missed_vcf_fpath = add_suffix(filtered_vcf_fpath, 'cosmic')
    # shutil.copy(filtered_vcf_fpath, cosmic_missed_vcf_fpath)

    report.add_record('Median depth', sample.median_cov)
    report.add_record('Variants', filtered_vcf_fpath)

    for r in regions:
        for vcf_db in vcf_dbs:
            if vcf_db.descriptive_name not in r.missed_by_db:
                r.missed_by_db[vcf_db.descriptive_name] = 0
            r.missed_by_db[vcf_db.descriptive_name] += \
                vcf_db.missed_vars_by_region_info.get((r.chrom, r.start, r.end, r.feature), 0)

    rec_by_name = {metric.name: Record(metric, []) for metric in regions_metrics}
    for rec in rec_by_name.values():
        report.records.append(rec)

    _proc_regions(regions, _fill_in_record_info, rec_by_name, sample.median_cov)

    return report


def _add_vcf_header(source_vcf_fpath, dest_vcf_fpath):
    with open(source_vcf_fpath) as inp, open(dest_vcf_fpath, 'w') as out:
        for line in inp:
            if line.startswith('#'):
                out.write(line)
    with open(dest_vcf_fpath + '_tmp') as inp, open(dest_vcf_fpath, 'a') as out:
        for line in inp:
            out.write(line)
    os.remove(dest_vcf_fpath + '_tmp')
    info(dest_vcf_fpath)


def _find_missed_variants(cnf, vcf_db, caller_name, input_vcf_fpath, report_base_name, bed_fpath):
    if not vcf_db.vcf_fpath and not verify_file(vcf_db.vcf_fpath, vcf_db.descriptive_name):
        info('Skipping counting variants from ' + vcf_db.descriptive_name)
        return None, None

    info('Counting missed variants from ' + vcf_db.descriptive_name)
    db_in_roi_vcf_fpath = join(cnf.work_dir, report_base_name + '_' + vcf_db.name + '_roi.vcf')
    missed_vcf_fpath = join(cnf.output_dir, report_base_name + '_' + caller_name + '_' + vcf_db.name + '_missed.vcf')
    annotated_regions_bed_fpath = join(cnf.work_dir, report_base_name + '_' + caller_name + '_' + vcf_db.name + '.bed')

    bedtools = get_tool_cmdline(cnf, 'bedtools')

    db_vcf_fpath = vcf_db.vcf_fpath
    # Take variants from DB VCF that lie in BED
    # (-wa means to take whole regions from a, not the parts that overlap)
    cmdline = '{bedtools} intersect -wa -a {db_vcf_fpath} -b {bed_fpath}'.format(**locals())
    call(cnf, cmdline, output_fpath=db_in_roi_vcf_fpath + '_tmp')
    _add_vcf_header(db_vcf_fpath, db_in_roi_vcf_fpath)
    info()

    # Take variants from VCF missed in DB (but only in the regions from BED)
    cmdline = '{bedtools} intersect -v -a {input_vcf_fpath} -b {db_in_roi_vcf_fpath}'.format(**locals())
    call(cnf, cmdline, output_fpath=missed_vcf_fpath + '_tmp')
    _add_vcf_header(input_vcf_fpath, missed_vcf_fpath)
    info()

    # Take variants from VCF missed in DB (but only in the regions from BED)
    cmdline = '{bedtools} intersect -c -a {bed_fpath} -b {missed_vcf_fpath}'.format(**locals())
    call(cnf, cmdline, output_fpath=annotated_regions_bed_fpath)
    info()

    return missed_vcf_fpath, annotated_regions_bed_fpath


# def _count_variants(region, vcf_fpath):
#     num = 0
#     with open(vcf_fpath) as inp:
#         reader = vcf_parser.Reader(inp)
#         for rec in reader:
#             if region.start <= rec.POS <= region.end:
#                 num += 1
#     return num


def _fill_in_record_info(region, rec_by_name, median_depth):
    rec_by_name['Sample'].value.append(         region.sample_name)
    rec_by_name['Chr'].value.append(            region.chrom)
    rec_by_name['Start'].value.append(          region.start)
    rec_by_name['End'].value.append(            region.end)
    rec_by_name['Gene'].value.append(           region.gene_name)
    rec_by_name['Feature'].value.append(        region.feature)
    rec_by_name['Size'].value.append(           region.get_size())
    rec_by_name['AvgDepth'].value.append(      region.avg_depth)
    rec_by_name['StdDev'].value.append(         region.std_dev)
    rec_by_name['Wn 20% of Mean'].value.append( region.percent_within_normal)
    rec_by_name['Depth/Median'].value.append(   region.avg_depth / median_depth if median_depth != 0 else None)
    rec_by_name['Total variants'].value.append( region.var_num)

    for db_descriptive_name, num in region.missed_by_db.items():
        rec_by_name[db_descriptive_name + ' missed'].value.append(num)

    for depth_thres, percent in region.percent_within_threshs.items():
        rec_by_name['{}x'.format(depth_thres)].value.append(percent)

