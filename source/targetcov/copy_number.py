#!/usr/bin/env python

import math
from collections import defaultdict

#from source.utils import mean, median
from operator import itemgetter, attrgetter
from numpy import mean, median

# Normalize the coverage from targeted sequencing to CNV log2 ratio. The algorithm assumes the medium
# is diploid, thus not suitable for homogeneous samples (e.g. parent-child).
def norm_depths_by_gene(sample_mapped_reads, gene_depth):
    list_genes_info = _report_row_to_objects(gene_depth)
    norm_depths_by_seq_distr, med_depth = _get_norm_depths_by_sample(sample_mapped_reads, list_genes_info)
    factors_by_gene = _get_factors_by_gene(list_genes_info, med_depth)
    norm_depths_by_gene = defaultdict(dict)
    norm2 = defaultdict(dict)
    norm3 = defaultdict(dict)
    median_depth_by_sample = dict()

    for sample_name in sample_mapped_reads:
        depths = [norm_depth_by_sample[sample_name] for gene, norm_depth_by_sample in norm_depths_by_seq_distr.items()]
        median_depth_by_sample[sample_name] = median(depths) if depths else 0

    for gene, norm_depth_by_sample in norm_depths_by_seq_distr.items():
        for sample in sample_mapped_reads:

            if sample not in norm_depth_by_sample:
                continue
            #norm1b
            norm_depths_by_gene[gene][sample] = norm_depth_by_sample[sample] * factors_by_gene[gene] + 0.1
            norm2[gene][sample] = math.log(norm_depths_by_gene[gene][sample] / med_depth, 2) if med_depth else 0

            norm3[gene][sample] = math.log(norm_depths_by_gene[gene][sample] / median_depth_by_sample[sample], 2) if \
                median_depth_by_sample[sample] else 0

    return _get_report_data(list_genes_info, norm2, norm3, norm_depths_by_gene, norm_depths_by_seq_distr)


# mean_reads =  mean for all the samples mapped reads
# return factor_by_sample = sample mapped read/mean for all the samples mapped reads
def _get_factors_by_sample(mapped_reads_by_sample):
    factor_by_sample = dict()
    mapped_reads_no_undeter = _removekey(mapped_reads_by_sample, "Undetermined")
    mean_reads = mean(mapped_reads_no_undeter.values())
    for k, v in mapped_reads_by_sample.items():
        factor_by_sample[k] = mean_reads / v if v else 0
    return factor_by_sample


#med_depth = median for all gene
#return factor_by_gene = median gene for all samples/ med_depth
def _get_factors_by_gene(gene_records, med_depth):
    min_depth_by_genes = defaultdict(list)

    [min_depth_by_genes[gene_record.name].append(gene_record.min_depth) for gene_record in gene_records]
    factors_by_gene = dict()
    for k, v in min_depth_by_genes.items():
        med = median(v)
        factors_by_gene[k] = med_depth / med if med else 0
    return factors_by_gene


#MeanDepth_Norm1
# gene -> { sample -> [] }
def _get_norm_depths_by_sample(mapped_reads_by_sample, record_by_sample):
    norm_depths = defaultdict(dict)
    depths = []
    factor_by_sample = _get_factors_by_sample(mapped_reads_by_sample)
    for sample, factor in factor_by_sample.items():
        for gene_info in [gene_infos for gene_infos in record_by_sample if gene_infos.sample_name == sample]:
            depth_by_factor = gene_info.min_depth * factor
            norm_depths[gene_info.name][sample] = depth_by_factor
            depths.append(depth_by_factor)
    return norm_depths, median(depths)


def _get_report_data(list_genes_info, norm2, norm3, norm_depths_by_gene, norm_depths_by_seq_distr):
    header = ["Sample", "Gene", "Chr", "Start", "Stop", "Length", "MeanDepth", "MeanDepth_Norm1",
              "MeanDepth_Norm2", "log2Ratio_norm1", "log2Ratio_norm2"]
    report_data = []

    for gene_info in list_genes_info:
        gene_name = gene_info.name
        sample = gene_info.sample_name
        if sample != "Undetermined":
            report_data.append(map(str,
                                   (sample, gene_name, gene_info.chrom, gene_info.start_position,
                                    gene_info.end_position, gene_info.size,
                                    '{0:.3f}'.format(gene_info.min_depth),
                                    '{0:.3f}'.format(norm_depths_by_seq_distr[gene_name][sample]),
                                    '{0:.3f}'.format(norm_depths_by_gene[gene_name][sample]),
                                    '{0:.3f}'.format(norm2[gene_name][sample]),
                                    '{0:.3f}'.format(norm3[gene_name][sample]))))

    report_data = sorted(report_data, key=itemgetter(2, 1))
    report_data = [header] + report_data

    return report_data


def _report_row_to_objects(gene_depth):
    gene_details = []
    for read in gene_depth:
        gene_details.append(GeneDetail(*read))
    return gene_details


def _removekey(d, key):
    r = dict(d)
    del r[key]
    return r


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

    def __repr__(self):
        return repr(
            (self.name, self.sample_name))

