#!/usr/bin/env python

import math
from collections import defaultdict

from source.utils import mean, median


def _get_factors_by_sample(mapped_reads_by_sample):
    factor_by_sample = dict()
    mean_reads = mean(mapped_reads_by_sample.values())
    for k, v in mapped_reads_by_sample.items():
        factor_by_sample[k] = mean_reads / v if v else 0
    return factor_by_sample


def _get_factors_by_gene(gene_records, med_depth):
    min_depth_by_genes = defaultdict(list)
    [min_depth_by_genes[gene_record.name].append(gene_record.min_depth) for gene_record in gene_records]
    factors_by_gene = dict()
    for k, v in min_depth_by_genes.items():
        med = median(v)
        factors_by_gene[k] = med_depth / med if med else 0
    return factors_by_gene


def get_norm_depths_by_seq_distr(mapped_reads_by_sample, record_by_sample):
    norm_depths = defaultdict(dict)  # gene -> { sample -> [] }

    factor_by_sample = _get_factors_by_sample(mapped_reads_by_sample)

    for sample, factor in factor_by_sample.items():
        for gene_info in [gene_infos for gene_infos in record_by_sample if gene_infos.sample_name == sample]:
            norm_depths[gene_info.name][sample] = gene_info.min_depth * factor

    return norm_depths


def get_report_data(list_genes_info, norm2, norm3, norm_depths_by_gene, norm_depths_by_seq_distr):
    report_data = []
    for gene_info in list_genes_info:
        gene_name = gene_info.name
        sample = gene_info.sample_name
        report_data.append(map(str,[sample, gene_name, gene_info.chrom, gene_info.start_position, gene_info.end_position, gene_info.size,
                gene_info.min_depth, norm_depths_by_seq_distr[gene_name][sample],
                norm_depths_by_gene[gene_name][sample],
                norm2[gene_name][sample], norm3[gene_name][sample]]))
    return report_data


def run_copy_number(sample_mapped_reads, gene_depth):

    list_genes_info = report_row_to_objects(gene_depth)

    med_depth = median([rec.min_depth for rec in list_genes_info])
    norm_depths_by_seq_distr = get_norm_depths_by_seq_distr(sample_mapped_reads, list_genes_info)
    factors_by_gene = _get_factors_by_gene(list_genes_info, med_depth)
    norm_depths_by_gene = defaultdict(dict)
    norm2 = defaultdict(dict)
    norm3 = defaultdict(dict)
    median_depth_by_sample = dict()

    for sample_name in sample_mapped_reads:
        median_depth_by_sample[sample_name] = median([gene_info.min_depth
                                                      for gene_info in list_genes_info
                                                      if gene_info.sample_name == sample_name])

    for gene, norm_depth_by_sample in norm_depths_by_seq_distr.items():

        for gene_info in list_genes_info:
            sample = gene_info.sample_name
            if sample not in norm_depth_by_sample:
                continue
            gene_norm_depth = norm_depth_by_sample[sample] * factors_by_gene[gene] + 0.1

            norm_depths_by_gene[gene][sample] = gene_norm_depth

            norm2[gene][sample] = math.log(gene_norm_depth / med_depth, 2) if med_depth else 0

            norm3[gene][sample] = math.log(gene_norm_depth / median_depth_by_sample[sample], 2)

    return get_report_data(list_genes_info, norm2, norm3, norm_depths_by_gene, norm_depths_by_seq_distr)


def report_row_to_objects(gene_depth):
    gene_details = []
    for read in gene_depth:
        gene_details.append(GeneDetail(*read))
    return gene_details


class GeneDetail():
    def __init__(self, sample_name=None, chrom=None, start_position=None, end_position=None, name=None,
                 type="Gene-Amplicon", size=None, min_depth=None):
        self.sample_name = sample_name
        self.chrom = chrom
        self.start_position = int(start_position)
        self.end_position = int(end_position)
        self.name = name
        self.type = type
        self.size = int(size)
        self.min_depth = float(min_depth)

    def __str__(self):
        values = [self.sample_name, self.chrom, self.start_position, self.end_position, self.name, self.type, self.size,
                  self.min_depth]
        return '"' + '\t'.join(map(str, values)) + '"'




