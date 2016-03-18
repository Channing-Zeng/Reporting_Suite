from collections import OrderedDict
from genericpath import isfile
from os.path import join, basename, splitext
from random import random
from time import sleep

from source.calling_process import call
from source.file_utils import splitext_plus, add_suffix, safe_mkdir, file_transaction, verify_file, intermediate_fname
from source.logger import critical, info, warn, send_email, err
from source.qsub_utils import submit_job
from source.targetcov.bam_and_bed_utils import verify_bam, verify_bed, bam_to_bed
from source.targetcov.Region import Region
from source.targetcov.bam_and_bed_utils import count_bed_cols, bedtools_version
from source.tools_from_cnf import get_system_path
from source.utils import get_chr_len_fpath, get_chr_lengths


class BedCov:
    def __init__(self, chrom, chrom_len, bedcov_output_fpath):
        self.chrom = chrom
        self.chrom_len = chrom_len
        self.bedcov_output_fpath = bedcov_output_fpath


def split_bed_by_chrom(cnf, bed_fpath):
    info('Splitting the BED file ' + bed_fpath + ' by chromosome: ', ending='')
    bed_fpath_by_chrom = dict()
    cur_chr_f = None
    cur_chr = None

    with open(bed_fpath) as f:
        for l in f:
            fs = l.strip().split('\t')
            if fs:
                if fs[0] != cur_chr:
                    if cur_chr:
                        info(str(cur_chr), ending=', ', print_date=False)
                    cur_chr = fs[0]
                    cur_chr_fpath = intermediate_fname(cnf, bed_fpath, cur_chr)
                    cur_chr_f = open(cur_chr_fpath, 'w')
                    bed_fpath_by_chrom[cur_chr] = cur_chr_fpath
                cur_chr_f.write(l)
    info('Done.', print_date=False)
    return bed_fpath_by_chrom


# def bedcoverage_hist_stats(cnf, sample_name, bed, bam):
#     if not bam or not bed:
#         info()
#         msgs = []
#         if not bam: msgs.append('BAM file is required.')
#         if not bed: msgs.append('BED file is required.')
#         if msgs:
#             critical(msgs)
#
#     # bedtools = get_system_path(cnf, 'bedtools')
#     # gzip = get_system_path(cnf, 'gzip')
#     # chr_lengths = get_chr_len_fpath(cnf)
#
#     # regions = []
#     # bamtools = get_system_path(cnf, 'bamtools')
#     # if not bamtools:
#     info('Running bedcoverage -hist...')
#     bedcov_output_fpath = launch_bedcoverage_hist(cnf, bed, bam)
#     info()
#     info('Analysing bedcoverage -hist output...')
#     regions = summarize_bedcoverage_hist_stats(bedcov_output_fpath, sample_name, count_bed_cols(bed))
#
#     # else:
#     #     chroms = get_chr_lengths(cnf).keys()
#     #
#     #     bed_fpath_by_chrom = split_bed_by_chrom(cnf, bed_fpath)
#     #
#     #     stub = join(cnf.work_dir, basename(splitext_plus(bam_fpath)[0]))
#     #     if cnf.reuse_intermediate and all(verify_bam(stub + '.REF_' + chrom + '.bam', silent=True) for chrom in chroms):
#     #         info('BAM ' + bam_fpath + ' is split, reusing...')
#     #     else:
#     #         info('Splitting the BAM file, writing as ' + stub + '.REF_#.bam')  # TODO do spltting once in a run (not for exome then for target)
#     #         cmdline = '{bamtools} split -in {bam_fpath} -stub {stub} -reference'.format(**locals())
#     #         call(cnf, cmdline)
#     #
#     #     for chrom in chroms:
#     #         chrom_bed_fpath = verify_bed(bed_fpath_by_chrom.get(chrom), silent=True)
#     #         chrom_bam_fpath = verify_bam(stub + '.REF_' + chrom + '.bam', silent=True)
#     #
#     #         if not chrom_bed_fpath:
#     #             info('No regions for ' + chrom)
#     #         if not chrom_bam_fpath:
#     #             info('No coverage for ' + chrom)
#     #         if chrom_bed_fpath and chrom_bam_fpath:
#     #             bedcov_output_fpath = launch_bedcoverage_hist(cnf, chrom_bed_fpath, chrom_bam_fpath)
#     #             if not verify_file(bedcov_output_fpath):
#     #                 info('No coverage for ' + chrom)
#     #             else:
#     #                 info('Anylising bedcoverage output for ' + str(chrom) + '...')
#     #                 rs = summarize_bedcoverage_hist_stats(bedcov_output_fpath, sample_name, bed_col_num)
#     #                 regions.extend(rs)
#
#         # with open(bedcov_output_fpath, 'w') as f:
#         #     for chrom, bedov_output in bedcov_by_chrom.items():
#
#     # for r in regions:
#     return regions


# TODO:
# split bam and bed by chromosome
# run bedtools hist
# merge. for general stats, merge more sofisticated.


# def grep(cnf, fpath, pattern):
#     res =


# def launch_bedcoverage_hist(cnf, bed, bam, bedcov_output_fpath=None, qsub=False, **kwargs):
#     # import pybedtools
#     # bed = pybedtools.BedTool(bed_fpath)
#     # bam = pybedtools.BedTool(bam_fpath)
#     # return bed.coverage(bam)
#
#     if not bedcov_output_fpath:
#         bedcov_output_fpath = join(cnf.work_dir,
#             splitext_plus(basename(bed))[0] + '__' +
#             splitext_plus(basename(bam))[0] + '_bedcov_output.txt')
#
#     if cnf.reuse_intermediate and verify_file(bedcov_output_fpath, silent=True):
#         info(bedcov_output_fpath + ' exists, reusing.')
#         if qsub:
#             return None
#         else:
#             return bedcov_output_fpath
#
#     bedtools = get_system_path(cnf, 'bedtools')
#     chr_lengths = get_chr_len_fpath(cnf)
#
#     if bam.endswith('bam'):
#         bam = bam_to_bed(cnf, bam)
#
#     v = bedtools_version(bedtools)
#     if v and v >= 24:
#         cmdline = '{bedtools} coverage -sorted -g {chr_lengths} -a {bed} -b {bam} -hist'.format(**locals())
#     else:
#         cmdline = '{bedtools} coverage -a {bam} -b {bed} -hist'.format(**locals())
#
#     res = None
#     tries = 0
#     MAX_TRIES = 2
#     WAIT_MINUTES = int(random() * 60) + 30
#     err_fpath = join(cnf.work_dir, 'bedtools_cov_' + splitext(basename(bedcov_output_fpath))[0] + '.err')
#
#     if qsub:
#         job_name = splitext_plus(basename(bed))[0] + '__' +\
#                    splitext_plus(basename(bam))[0] + '_bedcov'
#         return submit_job(cnf, cmdline, job_name, output_fpath=bedcov_output_fpath, **kwargs)
#     else:
#         return call(cnf, cmdline, bedcov_output_fpath, exit_on_error=False)





def summarize_bedcoverage_hist_stats(sambamba_depth_output_fpath, sample_name, bed_col_num):
    """
    :param sambamba_depth_fpath: file path
    :param sample_name: sting
    :param bed_col_num: int
    :return: regions, total_region, max_depth
        regions: [Region(
            sample_name: string
            chrom: string
            start: int
            end: int
            size: int
            gene_name: string
            bases_by_depth: dict(depth:int->bases:int)
            extra_fields: [string])]
        total_region: Region(--|--)
        max_depth: int
    """
    read_count_col = None
    mean_cov_col = None
    min_depth_col = None
    std_dev_col = None
    total_regions_count = 0

    cur_unannotated_gene = None
    info('Reading coverage statistics...')

    with open(sambamba_depth_output_fpath) as sambabma_depth_file:
        for line in sambabma_depth_file:
            if line.startswith('#'):
                read_count_col = line.split('\t').index('readCount')
                mean_cov_col = line.split('\t').index('meanCoverage')
                min_depth_col = line.split('\t').index('minDepth')
                std_dev_col = line.split('\t').index('stdDev')
                continue
            line_tokens = line.replace('\n', '').split()
            chrom = line_tokens[0]
            start, end = map(int, line_tokens[1:3])
            region_size = end - start
            gene_name = line_tokens[3] if read_count_col != 3 else None
            ave_depth = float(line_tokens[mean_cov_col])
            min_depth = int(line_tokens[min_depth_col])
            std_dev = float(line_tokens[std_dev_col])
            rates_within_threshs = line_tokens[std_dev_col + 1:-1]

            extra_fields = tuple(line_tokens[4:read_count_col]) if read_count_col > 4 else ()

            region = Region(
                sample_name=sample_name, chrom=chrom,
                start=start, end=end, size=region_size,
                avg_depth=ave_depth,
                gene_name=gene_name, extra_fields=extra_fields)

            region.rates_within_threshs = OrderedDict((depth, float(rate) / 100.0) for (depth, rate) in zip(depth_thresholds, rates_within_threshs))
            region.min_depth = min_depth
            region.std_dev = std_dev

            region.feature = 'Amplicon'
            if gene_name != '.':
                cur_unannotated_gene = None
                gene = gene_by_name_and_chrom[(gene_name, chrom)]
                if (gene.gene_name, gene.chrom) not in ready_to_report_set:
                    ready_to_report_genes.append(gene)
                    ready_to_report_set.add((gene.gene_name, gene.chrom))
                gene.add_amplicon(region)
            else:
                if cur_unannotated_gene is None:
                    cur_unannotated_gene = GeneInfo(sample_name=sample_name, gene_name=gene_name, chrom=chrom, feature='NotAnnotatedSummary')
                    ready_to_report_genes.append(cur_unannotated_gene)
                cur_unannotated_gene.add_amplicon(region)

            row = [region.chrom, region.start, region.end, region.get_size(), region.gene_name, region.strand,
                   region.feature, region.biotype, region.transcript_id, region.min_depth, region.avg_depth, region.std_dev,
                   region.rate_within_normal]
            row = [Metric.format_value(val, human_readable=True) for val in row]
            rates = [Metric.format_value(val, unit='%', human_readable=True) for val in region.rates_within_threshs.values()]
            row.extend(rates)
            col_widths = [max(len(v), w) for v, w in izip(row, col_widths)]

            total_regions_count += 1
            if total_regions_count > 0 and total_regions_count % 10000 == 0:
                 info('  Processed {0:,} regions'.format(total_regions_count))
