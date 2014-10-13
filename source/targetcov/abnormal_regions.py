# coding=utf-8
from collections import OrderedDict
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


def make_abnormal_regions_reports(cnf, bcbio_structure, sample, filtered_vcf_by_callername=None):
    gene_rep_fpath = bcbio_structure.get_gene_report_fpaths_by_sample().get(sample.name)

    regions = _read_regions(gene_rep_fpath)

    abnormal_regions_reports = _generate_abnormal_regions_reports(cnf, sample, regions,
        filtered_vcf_by_callername or dict())

    return abnormal_regions_reports


def _read_regions(gene_report_fpath):
    regions = []
    with open(gene_report_fpath) as f:
        for line in f:
            chrom, start, end, gene, feature, size, avg_depth, std_dev, percent_within_normal = line.strip().split()
            region = Region(
                chrom,
                int(start),
                int(end),
                gene,
                feature,
                int(size),
                float(avg_depth),
                float(std_dev),
                float(percent_within_normal)
            )
            regions.append(region)
    return regions


def _generate_abnormal_regions_reports(cnf, sample, regions, filtered_vcf_by_callername):
    info('Calculation median coverage...')
    median_cov = median((r.avg_depth for r in regions))
    if median_cov == 0:
        err('Median coverage is 0')
        return None, None

    info('Median: ' + str(median_cov))
    info()

    info('Extracting abnormally covered regions.')
    minimal_cov = min(cnf.coverage_reports.min_cov_factor * median_cov, cnf.coverage_reports.min_cov)
    maximal_cov = cnf.coverage_reports.max_cov_factor * median_cov

    info('Assuming abnormal if below min(median*' +
         str(cnf.coverage_reports.min_cov_factor) + ', ' +
         str(cnf.coverage_reports.min_cov) + ') = ' + str(minimal_cov) +
         ', or above median*' + str(cnf.coverage_reports.max_cov_factor) +
         ' = ' + str(maximal_cov))

    info('Classifying regions...')
    low_regions, high_regions = [], []
    _proc_regions(regions, _classify_region, low_regions, high_regions,
                  median_cov, minimal_cov, maximal_cov)

    abnormal_regions_reports = []
    if filtered_vcf_by_callername:  # No VCFs passed, so no need to find missed variants.
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

                report_base_name = cnf.name + '_' + caller_name + '.' + BCBioStructure.targetseq_name + f_basename

                abnormal_regions_fpath = _save_regions_to_bed(cnf, regions, report_base_name)
                for vcf_db in vcf_dbs:
                    vcf_db.missed_vcf_fpath, \
                    vcf_db.annotated_bed = _find_missed_variants(cnf, vcf_db,
                        caller_name, vcf_fpath, report_base_name, abnormal_regions_fpath)

                report = _make_flagged_region_report(cnf, sample, vcf_fpath, caller_name, vcf_dbs)

                regions_html_rep_fpath = report.save_html(cnf.output_dir, report_base_name,
                    caption='Regions with ' + kind + ' coverage for ' + caller_name)

                regions_txt_rep_fpath = report.save_txt(cnf.output_dir, report_base_name,
                    [report.metric_storage.sections_by_name[kind + '_cov']])

                abnormal_regions_reports.append(regions_txt_rep_fpath)
                info('Too ' + kind + ' covered regions (total ' + str(len(regions)) + ') saved into:')
                info('  ' + str(regions_html_rep_fpath))
                info('  ' + regions_txt_rep_fpath)

    return abnormal_regions_reports


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


def _make_flagged_region_report(cnf, sample, filtered_vcf_fpath, caller_name, vcf_dbs):
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
        regions_metrics.append(Metric('{}x'.format(thres), unit='%', description='Bases covered by at least {} reads'.format(thres)))

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
        vcf_db.missed_vars_by_region_info = OrderedDict()  # indexed by name (chr, start, end, feature)

        with open(vcf_db.annotated_bed) as f:
            for line in f:
                print line
                tokens = line.strip().split()
                chrom, start, end, gene, feature, missed_vars = tokens[:5]
                vcf_db.missed_vars_by_region_info[(chrom, start, end, feature)] = missed_vars

        if not vcf_db.missed_vcf_fpath and not verify_file(vcf_db.missed_vcf_fpath):
            info('Skipping counting variants from ' + vcf_db.descriptive_name)
            report.add_record(vcf_db.descriptive_name + ' missed variants', None)
        else:
            info('Counting missed variants from ' + vcf_db.descriptive_name)
            report.add_record(vcf_db.descriptive_name + ' missed variants',
                              vcf_db.missed_vcf_fpath)

        info()

    # cosmic_missed_vcf_fpath = add_suffix(filtered_vcf_fpath, 'cosmic')
    # shutil.copy(filtered_vcf_fpath, cosmic_missed_vcf_fpath)
    #
    # report.add_record('Median depth', median_depth)
    # report.add_record('Variants', filtered_vcf_fpath)
    #
    # rec_by_name = {metric.name: Record(metric, []) for metric in regions_metrics}
    # for rec in rec_by_name.values():
    #     report.records.append(rec)
    #
    # _proc_regions(regions, _fill_in_record_info, rec_by_name, median_depth)

    return report


def _add_vcf_header(source_vcf_fpath, dest_vcf_fpath):
    with open(source_vcf_fpath) as inp, open(dest_vcf_fpath, 'w') as out:
        for line in inp:
            if line.startswith('#'):
                out.write(line)
    with open(dest_vcf_fpath + '_tmp') as inp, open(dest_vcf_fpath, 'a') as out:
        for line in inp:
            out.write(line)


def _find_missed_variants(cnf, vcf_db, caller_name, input_vcf_fpath, report_base_name, bed_fpath):
    if not vcf_db.vcf_fpath and not verify_file(vcf_db.vcf_fpath):
        info('Skipping counting variants from ' + vcf_db.descriptive_name)
        return None

    info('Counting missed variants from ' + vcf_db.descriptive_name)
    db_in_roi_vcf_fpath = join(cnf.work_dir, report_base_name + '_' + vcf_db.name + '_roi.vcf')
    missed_vcf_fpath = join(cnf.output_dir, report_base_name + '_' + caller_name + '_' + vcf_db.name + '_missed.vcf')
    annotated_regions_fpath = join(cnf.work_dir, report_base_name + '_' + caller_name + '_' + vcf_db.name + '.bed')

    bedtools = get_tool_cmdline(cnf, 'bedtools')

    db_vcf_fpath = vcf_db.vcf_fpath
    # Take variants from DB VCF that lie in BED
    # (-wa means to take whole regions from a, not the parts that overlap)
    cmdline = '{bedtools} intersect -wa -a {db_vcf_fpath} -b {bed_fpath}'.format(**locals())
    call(cnf, cmdline, output_fpath=db_in_roi_vcf_fpath + '_tmp')
    _add_vcf_header(db_vcf_fpath, db_in_roi_vcf_fpath)

    # Take variants from VCF missed in DB (but only in the regions from BED)
    cmdline = '{bedtools} intersect -v -a {input_vcf_fpath} -b {db_in_roi_vcf_fpath}'.format(**locals())
    call(cnf, cmdline, output_fpath=missed_vcf_fpath + '_tmp')
    _add_vcf_header(input_vcf_fpath, missed_vcf_fpath)

    # Take variants from VCF missed in DB (but only in the regions from BED)
    cmdline = '{bedtools} intersect -wao -a {bed_fpath} -b {missed_vcf_fpath}'.format(**locals())
    call(cnf, cmdline, output_fpath=missed_vcf_fpath + '_tmp')
    _add_vcf_header(input_vcf_fpath, annotated_regions_fpath)

    return missed_vcf_fpath, annotated_regions_fpath


# def _count_variants(region, vcf_fpath):
#     num = 0
#     with open(vcf_fpath) as inp:
#         reader = vcf_parser.Reader(inp)
#         for rec in reader:
#             if region.start <= rec.POS <= region.end:
#                 num += 1
#     return num


def _fill_in_record_info(region, rec_by_name, median_depth):
    rec_by_name['Sample'].value.append(         region.sample)
    rec_by_name['Chr'].value.append(            region.chrom)
    rec_by_name['Start'].value.append(          region.start)
    rec_by_name['End'].value.append(            region.end)
    rec_by_name['Gene'].value.append(           region.gene_name)
    rec_by_name['Feature'].value.append(        region.feature)
    rec_by_name['Size'].value.append(           region.get_size())
    rec_by_name['MeanDepth'].value.append(      region.avg_depth)
    rec_by_name['StdDev'].value.append(         region.std_dev)
    rec_by_name['Wn 20% of Mean'].value.append( region.percent_within_normal)
    rec_by_name['Depth/Median'].value.append(   region.avg_depth / median_depth if median_depth != 0 else None)
    rec_by_name['Total variants'].value.append( region.var_num)

    for db_descriptive_name, num in region.missed_by_db.items():
        rec_by_name[db_descriptive_name + ' missed'].value.append(num)

    for depth_thres, bases in region.bases_within_threshs.items():
        if int(region.get_size()) == 0:
            percent = None
        else:
            percent = 100.0 * bases / region.get_size()
        rec_by_name['{}x'.format(depth_thres)].value.append(percent)



