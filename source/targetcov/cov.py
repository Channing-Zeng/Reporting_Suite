# coding=utf-8
from collections import OrderedDict
import copy
from itertools import izip, chain, repeat
import os
from os.path import join, basename, isfile
import shutil
import sys
from ext_modules import vcf_parser
from source.bcbio_structure import BCBioStructure, Sample
from source.calling_process import call, call_check_output, call_pipe
from source.file_utils import intermediate_fname, splitext_plus, add_suffix, verify_file

from source.logger import step_greetings, critical, info, err
from source.ngscat.bed_file import verify_bed
from source.reporting import Metric, SampleReport, FullReport, MetricStorage, ReportSection, write_txt_rows, \
    write_tsv_rows, Record, SquareSampleReport
from source.targetcov.Region import Region
from source.tools_from_cnf import get_tool_cmdline
from source.utils import get_chr_len_fpath, median


def generate_targetcov_reports(cnf, sample, filtered_vcf_by_callername=None):
    if not (sample.bam or sample.bed):
        if not sample.bam: err(sample.name + ': BAM file is required.')
        if not sample.bed: err(sample.name + ': BED file is required.')
        sys.exit(1)

    summary_report = None
    summary_report_txt_path = None
    summary_report_json_fpath = None
    summary_report_html_fpath = None
    gene_report_fpath = None
    abnormal_regions_reports = []

    info()
    info('Calculation of coverage statistics for the regions in the input BED file...')
    amplicons, combined_region, max_depth, total_bed_size = bedcoverage_hist_stats(cnf, sample.bam, sample.bed)

    chr_len_fpath = get_chr_len_fpath(cnf)
    exons_bed = cnf['genome']['exons']

    if 'summary' in cnf['coverage_reports']['report_types']:
        step_greetings('Target coverage summary report')

        summary_report = generate_summary_report(
            cnf, sample, chr_len_fpath,
            cnf.coverage_reports.depth_thresholds, cnf.padding,
            combined_region, max_depth, total_bed_size)

        summary_report_json_fpath = join(cnf.output_dir, cnf.name + '.' + BCBioStructure.targetseq_name + '.json')
        summary_report.dump(summary_report_json_fpath)

        summary_report_txt_fpath  = summary_report.save_txt (cnf.output_dir, cnf.name + '.' + BCBioStructure.targetseq_name)
        summary_report_html_fpath = summary_report.save_html(cnf.output_dir, cnf.name + '.' + BCBioStructure.targetseq_name)
        info()
        info('Saved to ')
        info('  ' + summary_report_txt_fpath)


    if 'genes' in cnf['coverage_reports']['report_types']:
        step_greetings('Analysing regions.')

        if not exons_bed:
            if cnf['reports'] == 'genes':
                critical('Error: no genes or exons specified for the genome in system config, '
                         'cannot run per-exon report.')
            else:
                err('Warning: no genes or exons specified for the genome in system config, '
                    'cannot run per-exon report.')
        else:
            # log('Annotating amplicons.')
            # annotate_amplicons(amplicons, genes_bed)

            info('Sorting exons BED file.')
            exons_bed = sort_bed(cnf, exons_bed)
            info('Getting the exons that overlap amplicons.')
            overlapped_exons_bed = intersect_bed(cnf, exons_bed, sample.bed)
            info('Adding other exons for the genes of overlapped exons.')
            target_exons_bed = _add_other_exons(cnf, exons_bed, overlapped_exons_bed)

            info('Calculation of coverage statistics for exons of the genes ovelapping with the input regions...')
            exons, _, _, _ = bedcoverage_hist_stats(cnf, sample.bam, target_exons_bed)
            for exon in exons:
                exon.gene_name = exon.extra_fields[0]

            info('Saving region coverage report...')
            gene_report_fpath, abnormal_regions_reports = _generate_region_cov_report(
                cnf, filtered_vcf_by_callername or dict(), sample, cnf.output_dir,
                cnf.name, cnf.coverage_reports.depth_thresholds, amplicons, exons)

    if summary_report:
        report = summary_report
        summary_report_html_fpath = report.save_html(cnf.output_dir,
            cnf.name + '.' + BCBioStructure.targetseq_name,
            caption='Target coverage statistics for ' + cnf.name)

        # info('\t' + summary_report_html_fpath)

    return summary_report_txt_path, summary_report_json_fpath, summary_report_html_fpath, \
           gene_report_fpath, abnormal_regions_reports


def _add_other_exons(cnf, exons_bed, overlapped_exons_bed):
    gene_names = set()
    with open(overlapped_exons_bed) as overl_f:
        for line in overl_f:
            if not line or not line.strip() or line.startswith('#') or len(line.split()) < 4:
                continue
            gene_names.add(line.split()[3])

    new_overlp_exons_bed = intermediate_fname(cnf, overlapped_exons_bed, 'by_genes')
    with open(exons_bed) as exons_f, \
         open(new_overlp_exons_bed, 'w') as new_overl_f:
        for line in exons_f:
            if not line or not line.strip() or line.startswith('#') or len(line.split()) < 4:
                new_overl_f.write(line)
            elif line.split()[3] in gene_names:
                new_overl_f.write(line)

    return sort_bed(cnf, new_overlp_exons_bed)


def get_records_by_metrics(records, metrics):
    _records = []
    for rec in records:
        if rec.metric.name in metrics:
            rec.metric = metrics[rec.metric.name]
            _records.append(rec)
    return _records


header_metric_storage = MetricStorage(
    general_section=ReportSection('general_section', '', [
        Metric('Bases in target', short_name='Target bp', common=True)
    ]),
    sections_by_name=OrderedDict(
        basic_metrics=ReportSection('basic_metrics', '', [
            Metric('Reads'),
            Metric('Mapped reads', short_name='Mapped'),
            Metric('Percentage of mapped reads', short_name='%', unit='%'),

            Metric('Unmapped reads', short_name='Unmapped', quality='Less is better'),
            Metric('Percentage of unmapped reads', short_name='%', unit='%', quality='Less is better'),

            Metric('Covered bases in target', short_name='Covered'),
            Metric('Percentage of target covered by at least 1 read', short_name='%', unit='%'),

            Metric('Reads mapped on target', short_name='Reads on trg'),
            Metric('Percentage of reads mapped on target', short_name='%', unit='%'),

            Metric('Reads mapped on padded target', 'On padded trg'),
            Metric('Percentage of reads mapped on padded target', short_name='%', unit='%'),

            Metric('Read bases mapped on target', short_name='Read bp on trg')
        ]),
        depth_metrics=ReportSection('depth_metrics', 'Target coverage depth', [
            Metric('Average target coverage depth', short_name='Avg'),
            Metric('Std. dev. of target coverage depth', short_name='Std dev', quality='Less is better'),
            Metric('Maximum target coverage depth', short_name='Max'),
            Metric('Percentage of target within 20% of mean depth', short_name='&#177;20% avg', unit='%')
        ])))


def generate_summary_report(cnf, sample, chr_len_fpath,
                            depth_thresholds, padding,
                            combined_region, max_depth, total_bed_size):

    for depth in depth_thresholds:
        name = 'Part of target covered at least by ' + str(depth) + 'x'
        header_metric_storage.add_metric(
            Metric(name, short_name=str(depth) + 'x', description=name, unit='%'),
            'depth_metrics')

    report = SampleReport(sample, metric_storage=header_metric_storage)

    info('* General coverage statistics *')
    info('Getting number of reads...')
    v_number_of_reads = number_of_reads(cnf, sample.bam)
    report.add_record('Reads', v_number_of_reads)

    info('Getting number of mapped reads...')
    v_mapped_reads = number_of_mapped_reads(cnf, sample.bam)
    v_percent_mapped = 100.0 * v_mapped_reads / v_number_of_reads if v_number_of_reads else None
    v_percent_unmapped = 100.0 * (v_number_of_reads - v_mapped_reads) / v_number_of_reads if v_number_of_reads else None
    report.add_record('Mapped reads', v_mapped_reads)
    report.add_record('Percentage of mapped reads', v_percent_mapped)
    report.add_record('Unmapped reads', v_number_of_reads - v_mapped_reads)
    report.add_record('Percentage of unmapped reads', v_percent_unmapped)
    info('')

    info('* Target coverage statistics *')
    report.add_record('Bases in target', total_bed_size)

    combined_region.sum_up(depth_thresholds)

    v_covered_bases_in_targ = combined_region.bases_within_threshs.items()[0][1]
    report.add_record('Covered bases in target', v_covered_bases_in_targ)

    v_percent_covered_bases_in_targ = 100.0 * v_covered_bases_in_targ / total_bed_size if total_bed_size else None
    report.add_record('Percentage of target covered by at least 1 read', v_percent_covered_bases_in_targ)
    info('Getting number of mapped reads on target...')

    v_mapped_reads_on_target = number_mapped_reads_on_target(cnf, sample.bed, sample.bam)
    report.add_record('Reads mapped on target', v_mapped_reads_on_target)

    v_percent_mapped_on_target = 100.0 * v_mapped_reads_on_target / v_mapped_reads if v_mapped_reads else None
    report.add_record('Percentage of reads mapped on target ', v_percent_mapped_on_target)

    info('Making bed file for padded regions...')
    padded_bed = get_padded_bed_file(cnf, sample.bed, chr_len_fpath, padding)

    info('Getting number of mapped reads on padded target...')
    v_reads_on_padded_targ = number_mapped_reads_on_target(cnf, padded_bed, sample.bam)
    report.add_record('Reads mapped on padded target', v_reads_on_padded_targ)

    v_percent_mapped_on_padded_target = 100.0 * v_reads_on_padded_targ / v_mapped_reads if v_mapped_reads else None
    report.add_record('Percentage of reads mapped on padded target', v_percent_mapped_on_padded_target)

    v_read_bases_on_targ = combined_region.avg_depth * total_bed_size  # sum of all coverages
    report.add_record('Read bases mapped on target', v_read_bases_on_targ)

    info('')
    report.add_record('Average target coverage depth', combined_region.avg_depth)
    report.add_record('Std. dev. of target coverage depth', combined_region.std_dev)
    report.add_record('Maximum target coverage depth', max_depth)
    report.add_record('Percentage of target within 20% of mean depth', combined_region.percent_within_normal)

    for depth, bases in combined_region.bases_within_threshs.items():
        percent_val = 100.0 * bases / total_bed_size if total_bed_size else 0
        report.add_record('Part of target covered at least by ' + str(depth) + 'x', percent_val)

    return report


def _generate_region_cov_report(cnf, filtered_vcf_by_callername, sample, output_dir,
                                sample_name, depth_threshs, amplicons, exons):
    for ampl in amplicons:
        ampl.feature = 'Amplicon'
        ampl.sample = sample_name

    for exon in exons:
        exon.feature = 'Exon'
        exon.sample = sample_name

    info('Groupping exons per gene...')
    exon_genes = _get_exons_merged_by_genes(cnf, exons)
    info()

    info('Groupping amplicons per gene...')
    amplicon_genes_by_name = _get_amplicons_merged_by_genes(amplicons, exon_genes)
    info()

    regions = []
    info('Combining...')
    i = 0
    for exon_gene in exon_genes:
        if i and i % 10000 == 0:
            info('Processed {0:,} genes, current gene {1}'.format(i, exon_gene.gene_name))
        i += 1

        for exon in exon_gene.subregions:
            regions.append(exon)
        regions.append(exon_gene)

        amplicon_gene = amplicon_genes_by_name.get(exon_gene.gene_name)
        if amplicon_gene:
            for amplicon in amplicon_gene.subregions:
                regions.append(amplicon)
            regions.append(amplicon_gene)
    info('Processed {0:,} genes.'.format(i))

    info()
    info('Summing up region stats...')
    i = 0
    for region in regions:
        i += 1
        if i % 10000 == 0:
            info('Processed {0:,} regions.'.format(i))
        region.sum_up(depth_threshs)
        # if region.percent_within_normal < 100:
        #     bad_regions.append(region)

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
         str(cnf.coverage_reports.min_cov_factor) +
         str(minimal_cov) + ', ' + str(cnf.coverage_reports.min_cov) +
         ') or above median*' + str(cnf.coverage_reports.max_cov_factor) +
         ' = ' + str(maximal_cov))

    low_regions, high_regions = [], []
    _proc_regions(regions, _classify_region, low_regions, high_regions,
                  median_cov, minimal_cov, maximal_cov)

    rows = _make_flat_region_report(regions, depth_threshs)

    # Saving gene details report
    gene_report_basename = cnf.name + '.' + \
                           BCBioStructure.targetseq_name + \
                           BCBioStructure.detail_gene_report_baseending
    txt_rep_fpath = write_txt_rows(rows, output_dir, gene_report_basename)
    tsv_rep_fpath = write_tsv_rows(rows, output_dir, gene_report_basename)
    info('')
    info('Regions (total ' + str(len(regions)) + ') saved into:')
    info('  ' + txt_rep_fpath)

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

                bed_fpath = _save_regions_to_bed(cnf, regions, report_base_name)
                for vcf_db in vcf_dbs:
                    vcf_db.missed_vcf_fpath = _prepare_missed_vcfs(cnf, vcf_db, caller_name, vcf_fpath, report_base_name, bed_fpath)

                report = _make_flagged_region_report(
                    cnf, sample, vcf_fpath, caller_name, regions, bed_fpath, vcf_dbs, depth_threshs)

                regions_html_rep_fpath = report.save_html(output_dir, report_base_name,
                    caption='Regions with ' + kind + ' coverage for ' + caller_name)

                regions_txt_rep_fpath = report.save_txt(output_dir, report_base_name,
                    [report.metric_storage.sections_by_name[kind + '_cov']])

                abnormal_regions_reports.append(regions_txt_rep_fpath)
                info('Too ' + kind + ' covered regions (total ' + str(len(regions)) + ') saved into:')
                info('  ' + str(regions_html_rep_fpath))
                info('  ' + regions_txt_rep_fpath)

    return tsv_rep_fpath, abnormal_regions_reports


def _classify_region(region, low_regions, high_regions, median_cov, minimal_cov, maximal_cov):
    if region.avg_depth < minimal_cov:
        region.cov_factor = region.avg_depth / median_cov
        low_regions.append(region)

    if region.avg_depth > maximal_cov:
        region.cov_factor = region.avg_depth / median_cov
        high_regions.append(region)


def _save_regions_to_bed(cnf, regions, f_basename):
    bed_fpath = join(cnf.work_dir, f_basename + '.bed')
    info()
    info('Saving regions to ' + bed_fpath)

    if isfile(bed_fpath):
        if cnf.reuse_intermediate:
            if not verify_bed(bed_fpath):
                sys.exit(1)
            return bed_fpath
        else:
            os.remove(bed_fpath)

    with open(bed_fpath, 'w') as f:
        for region in regions:
            f.write('\t'.join([region.chrom, str(region.start), str(region.end), str(region.avg_depth)]) + '\n')

    info('Saved to ' + bed_fpath)
    return bed_fpath


def _get_amplicons_merged_by_genes(amplicons, exon_genes):
    amplicon_genes_by_name = dict()

    i = 0
    for exon_gene in exon_genes:
        if i and i % 1000 == 0:
            info('Processed {0:,} regions, current gene {1}'.format(i, exon_gene.gene_name))
        i += 1

        for amplicon in amplicons:
            if exon_gene.intersect(amplicon):
                amplicon_gene = amplicon_genes_by_name.get(exon_gene.gene_name)
                if amplicon_gene is None:
                    amplicon_gene = Region(
                        sample=amplicon.sample, chrom=amplicon.chrom,
                        start=amplicon.start, end=amplicon.end, size=0,
                        feature='Gene-' + amplicon.feature,
                        gene_name=exon_gene.gene_name)
                    amplicon_genes_by_name[amplicon_gene.gene_name] = amplicon_gene
                amplicon_copy = copy.copy(amplicon)
                amplicon_gene.add_subregion(amplicon_copy)
                amplicon_copy.gene_name = amplicon_gene.gene_name

    info('Processed {0:,} regions.'.format(i))
    return amplicon_genes_by_name


def _get_exons_merged_by_genes(cnf, subregions):
    genes_by_name = dict()

    i = 0
    for exon in subregions:
        if not exon.gene_name:
            info()
            err('No gene name info in the record: ' + str(exon) + '. Skipping.')
            continue

        if i and i % 10000 == 0:
            info('Processed {0:,} exons, current gene {1}'.format(i, exon.gene_name))
        i += 1

        gene = genes_by_name.get(exon.gene_name)
        if gene is None:
            gene = Region(
                sample=exon.sample, chrom=exon.chrom,
                start=exon.start, end=exon.end, size=0,
                feature='Gene-' + exon.feature,
                gene_name=exon.gene_name)
            genes_by_name[exon.gene_name] = gene
        gene.add_subregion(exon)

    sorted_genes = sorted(genes_by_name.values(), key=lambda g: (g.chrom, g.start, g.end))

    info('Processed {0:,} exons.'.format(i))
    return sorted_genes


# def _make_region_report(sample, regions, depth_threshs):
#     metrics = [
#         Metric('Sample'),
#         Metric('Chr'),
#         Metric('Start'),
#         Metric('End'),
#         Metric('Gene'),
#         Metric('Feature'),
#         Metric('Size'),
#         Metric('Mean Depth'),
#         Metric('Std Dev', description='Coverage depth standard deviation'),
#         Metric('Wn 20% of Mean Depth', description='Persentage of the region that lies within 20% of an avarage depth.'),
#         Metric('Depth/Median', description='Average depth of coveraege of the region devided by median coverage of the sample'),
#     ]
#     for thres in depth_threshs:
#         metrics.append(Metric('{}x'.format(thres), description='Bases covered by at least {} reads'.format(thres)))
#
#     region_metric_storage = MetricStorage(sections=[ReportSection(metrics=metrics)])
#
#     report = SquareSampleReport(sample, metric_storage=region_metric_storage,
#         report_name='Amplicons and exons coverage statistics.')
#
#     i = 0
#     for region in regions:
#         i += 1
#         if i % 10000 == 0:
#             info('Processed {0:,} regions.'.format(i))
#
#         region.sum_up(depth_threshs)
#         report.add_row('Sample', region.sample)
#
#         # row = map(str, [region.sample,
#         # row = map(str, [region.sample,
#
#         # report
#         row = map(str, [region.sample,
#                         region.chrom,
#                         region.start,
#                         region.end,
#                         region.gene_name])
#         row += [region.feature]
#         row += [str(region.get_size())]
#         row += ['{0:.2f}'.format(region.avg_depth)]
#         row += ['{0:.2f}'.format(region.std_dev)]
#         row += ['{0:.2f}%'.format(region.percent_within_normal) if region.percent_within_normal is not None else '-']
#         if print_cov_factor:
#             row += ['{0:.2f}%'.format(region.cov_factor)]
#
#         for depth_thres, bases in region.bases_within_threshs.items():
#             if int(region.get_size()) == 0:
#                 percent_str = '-'
#             else:
#                 percent = 100.0 * bases / region.get_size()
#                 percent_str = '{0:.2f}%'.format(percent)
#             row.append(percent_str)
#
#         all_rows.append(row)
#         # max_lengths = map(max, izip(max_lengths, chain(map(len, line_fields), repeat(0))))
#     info('Processed {0:,} regions.'.format(i))
#     return all_rows


class VCFDataBase():
    def __init__(self, name, descriptive_name, genome_cnf):
        self.name = name
        self.descriptive_name = descriptive_name
        self.vcf_fpath = genome_cnf.get(name)


def _make_flagged_region_report(cnf, sample, filtered_vcf_fpath, caller_name, regions, bed_fpath,
                                vcf_dbs, depth_threshs):
    regions_metrics = [
        Metric('Sample'),
        Metric('Chr'),
        Metric('Start'),
        Metric('End'),
        Metric('Gene'),
        Metric('Feature'),
        Metric('Size'),
        Metric('MeanDepth'),
        Metric('StdDev', description='Coverage depth standard deviation'),
        Metric('Wn 20% of Mean', unit='%', description='Persentage of the region that lies within 20% of an avarage depth.'),
        Metric('Depth/Median', description='Average depth of coveraege of the region devided by median coverage of the sample'),
        Metric('Total variants')
    ]
    for vcf_db in vcf_dbs:
        regions_metrics.append(Metric(vcf_db.descriptive_name + ' missed'))

    for thres in depth_threshs:
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
        _proc_regions(regions, lambda r: r.sum_up(depth_threshs))

        if not vcf_db.missed_vcf_fpath and not verify_file(vcf_db.missed_vcf_fpath):
            info('Skipping counting variants from ' + vcf_db.descriptive_name)
            report.add_record(vcf_db.descriptive_name + ' missed variants', None)
        else:
            info('Counting missed variants from ' + vcf_db.descriptive_name)
            report.add_record(vcf_db.descriptive_name + ' missed variants', vcf_db.missed_vcf_fpath)

            def _(r): r.missed_by_db[vcf_db.descriptive_name] = _count_variants(r, vcf_db.vcf_fpath)
            _proc_regions(regions, _)
        info()

    def _(r): r.var_num = _count_variants(r, filtered_vcf_fpath)
    _proc_regions(regions, _)

    median_depth = median(r.avg_depth for r in regions)
    cosmic_missed_vcf_fpath = add_suffix(filtered_vcf_fpath, 'cosmic')

    shutil.copy(filtered_vcf_fpath, cosmic_missed_vcf_fpath)

    report.add_record('Median depth', median_depth)
    report.add_record('Variants', filtered_vcf_fpath)

    rec_by_name = {metric.name: Record(metric, []) for metric in regions_metrics}
    for rec in rec_by_name.values():
        report.records.append(rec)

    _proc_regions(regions, _fill_in_record_info, rec_by_name, median_depth)

    return report


def _add_vcf_header(source_vcf_fpath, dest_vcf_fpath):
    with open(source_vcf_fpath) as inp, open(dest_vcf_fpath, 'w') as out:
        for line in inp:
            if line.startswith('#'):
                out.write(line)
    with open(dest_vcf_fpath + '_tmp') as inp, open(dest_vcf_fpath, 'a') as out:
        for line in inp:
            out.write(line)


def _prepare_missed_vcfs(cnf, vcf_db, caller_name, input_vcf_fpath, report_base_name, bed_fpath):
    if not vcf_db.vcf_fpath and not verify_file(vcf_db.vcf_fpath):
        info('Skipping counting variants from ' + vcf_db.descriptive_name)
    else:
        info('Counting missed variants from ' + vcf_db.descriptive_name)
        db_in_roi_vcf_fpath = join(cnf.output_dir, report_base_name + '_' + vcf_db.name + '_roi.vcf')
        missed_vcf_fpath = join(cnf.output_dir, report_base_name + '_' + caller_name + '_' + vcf_db.name + '_missed.vcf')

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

        return missed_vcf_fpath


def _proc_regions(regions, fn, *args, **kwargs):
    i = 0
    for region in regions:
        i += 1

        fn(region, *args, **kwargs)

        if i % 10000 == 0:
            info('Processed {0:,} regions.'.format(i))

    if not i % 10000:
        info('Processed {0:,} regions.'.format(i))


def _count_variants(region, vcf_fpath):
    num = 0
    with open(vcf_fpath) as inp:
        reader = vcf_parser.Reader(inp)
        for rec in reader:
            if region.start <= rec.POS <= region.end:
                num += 1
    return num


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


def _make_flat_region_report(regions, depth_threshs):
    header_fields = ['Sample', 'Chr', 'Start', 'End', 'Gene', 'Feature', 'Size',
                     'Mean Depth', 'Standard Dev.', 'Within 20% of Mean']
    for thres in depth_threshs:
        header_fields.append('{}x'.format(thres))

    all_rows = [header_fields]

    i = 0
    for region in regions:
        i += 1
        if i % 10000 == 0:
            info('Processed {0:,} regions.'.format(i))

        region.sum_up(depth_threshs)

        row = map(str, [region.sample, region.chrom, region.start, region.end, region.gene_name])
        row += [region.feature]
        row += [str(region.get_size())]
        row += ['{0:.2f}'.format(region.avg_depth)]
        row += ['{0:.2f}'.format(region.std_dev)]
        row += ['{0:.2f}%'.format(region.percent_within_normal) if region.percent_within_normal is not None else '-']

        for depth_thres, bases in region.bases_within_threshs.items():
            if int(region.get_size()) == 0:
                percent_str = '-'
            else:
                percent = 100.0 * bases / region.get_size()
                percent_str = '{0:.2f}%'.format(percent)
            row.append(percent_str)

        all_rows.append(row)
        # max_lengths = map(max, izip(max_lengths, chain(map(len, line_fields), repeat(0))))
    info('Processed {0:,} regions.'.format(i))
    return all_rows


def bedcoverage_hist_stats(cnf, bam, bed):
    if not bam or not bed:
        err()
        if not bam: err('BAM file is required.')
        if not bed: err('BED file is required.')
        sys.exit(1)

    regions, max_depth, total_bed_size = [], 0, 0

    bedtools = get_tool_cmdline(cnf, 'bedtools')
    cmdline = '{bedtools} coverage -abam {bam} -b {bed} -hist'.format(**locals())
    bedcov_output = join(cnf.work_dir,
                         splitext_plus(basename(bed))[0] + '_' +
                         splitext_plus(basename(bam))[0] + '_bedcov_output.txt')
    call(cnf, cmdline, bedcov_output)

    _total_regions_count = 0

    with open(bedcov_output) as f:
        for next_line in f:
            if not next_line.strip() or next_line.startswith('#'):
                continue

            line_tokens = next_line.strip().split()
            chrom = line_tokens[0]
            start, end = None, None
            try:
                depth, bases, region_size = map(int, line_tokens[-4:-1])
            except:
                critical('Undexpected error: incorrect line in the coverageBed output:\n' + next_line)

            if next_line.startswith('all'):
                max_depth = max(max_depth, depth)
                total_bed_size += bases
                extra_fields = []
            else:
                start, end = map(int, line_tokens[1:3])
                extra_fields = line_tokens[3:-4]

            line_region_key_tokens = (None, chrom, start, end)

            if regions == [] or hash(line_region_key_tokens) != regions[-1].key():
                region = Region(sample=None, chrom=chrom,
                                start=start, end=end, size=region_size,
                                extra_fields=extra_fields)
                regions.append(region)

                _total_regions_count += 1

                if _total_regions_count > 0 and _total_regions_count % 100000 == 0:
                    info('  Processed {0:,} regions'.format(_total_regions_count))

            regions[-1].add_bases_for_depth(depth, bases)

    if _total_regions_count % 100000 != 0:
        info('Processed {0:,} regions'.format(_total_regions_count))

    return regions[:-1], regions[-1], max_depth, total_bed_size


def sort_bed(cnf, bed_fpath):
    bedtools = get_tool_cmdline(cnf, 'bedtools')
    cmdline = '{bedtools} sort -i {bed_fpath}'.format(**locals())
    output_fpath = intermediate_fname(cnf, bed_fpath, 'sorted')
    call(cnf, cmdline, output_fpath)
    return output_fpath


def intersect_bed(cnf, bed1, bed2):
    bed1_fname, _ = splitext_plus(basename(bed1))
    bed2_fname, _ = splitext_plus(basename(bed2))
    output_fpath = join(cnf['work_dir'], bed1_fname + '__' + bed2_fname + '.bed')
    bedtools = get_tool_cmdline(cnf, 'bedtools')
    cmdline = '{bedtools} intersect -a {bed1} -b {bed2} -u'.format(**locals())
    call(cnf, cmdline, output_fpath)
    return output_fpath


#TODO: works slow.
def number_of_mapped_reads(cnf, bam):
    samtools = get_tool_cmdline(cnf, 'samtools')
    cmdline = '{samtools} view -c -F 4 {bam}'.format(**locals())
    res = call_check_output(cnf, cmdline)
    return int(res)

#TODO: works slow.
def number_of_unmapped_reads(cnf, bam):
    samtools = get_tool_cmdline(cnf, 'samtools')
    cmdline = '{samtools} view -c -f 4 {bam}'.format(**locals())
    res = call_check_output(cnf, cmdline)
    return int(res)

#TODO: works slow.
def number_of_reads(cnf, bam):
    samtools = get_tool_cmdline(cnf, 'samtools')
    cmdline = '{samtools} view -c {bam}'.format(**locals())
    res = call_check_output(cnf, cmdline)
    return int(res)

#TODO: works slow.
def number_mapped_reads_on_target(cnf, bed, bam):
    samtools = get_tool_cmdline(cnf, 'samtools')
    cmdline = '{samtools} view -c -F 4 -L {bed} {bam}'.format(**locals())
    res = call_check_output(cnf, cmdline)
    return int(res)

#TODO: works slow.
def number_bases_in_aligned_reads(cnf, bam):
    samtools = get_tool_cmdline(cnf, 'samtools')
    cmdline = '{samtools} depth {bam}'.format(**locals())
    proc = call_pipe(cnf, cmdline)
    count = 0
    while True:
        coverage_line = proc.stdout.readline()
        if coverage_line:
            values = coverage_line.strip().split('\t')
            count += int(values[2])
    return count


def get_padded_bed_file(cnf, bed, genome, padding):
    bedtools = get_tool_cmdline(cnf, 'bedtools')
    cmdline = '{bedtools} slop -i {bed} -g {genome} -b {padding}'.format(**locals())
    output_fpath = intermediate_fname(cnf, bed, 'padded')
    call(cnf, cmdline, output_fpath)
    return output_fpath


def _bases_by_depth(depth_vals, depth_thresholds):
    bases_by_min_depth = {depth: 0 for depth in depth_thresholds}

    for depth_value in depth_vals:
        for threshold in depth_thresholds:
            if depth_value >= threshold:
                bases_by_min_depth[threshold] += 1

        return [100.0 * bases_by_min_depth[thres] / len(depth_vals) if depth_vals else 0
                for thres in depth_thresholds]


