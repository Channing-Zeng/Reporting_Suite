#!/usr/bin/env python

import math
from collections import defaultdict, OrderedDict
from source.logger import info, err

from source.utils import OrderedDefaultDict
# from numpy import median, mean
from source.utils import median, mean

# Normalize the coverage from targeted sequencing to CNV log2 ratio. The algorithm assumes the medium
# is diploid, thus not suitable for homogeneous samples (e.g. parent-child).


def run_copy_number(mapped_reads_by_sample, gene_depth):
    mapped_reads_by_sample = {k: v for k, v in mapped_reads_by_sample.items() if 'Undetermined' not in k}

    info('Parsing rows...')
    records = _report_row_to_objects(gene_depth)

    info('Calculating depths normalized by samples...')
    norm_depths_by_sample = _get_norm_depths_by_sample(mapped_reads_by_sample, records)

    med_depth = median([depth for gn, vs in norm_depths_by_sample.items() for sn, depth in vs.items()])

    info('Getting factors by genes...')
    factors_by_gene = _get_factors_by_gene(norm_depths_by_sample, med_depth)

    info('Getting median norm...')
    median_depth_by_sample = _get_med_norm_depths(mapped_reads_by_sample, norm_depths_by_sample)

    norm_depths_by_gene = OrderedDefaultDict(dict)
    norm2 = OrderedDefaultDict(dict)
    norm3 = OrderedDefaultDict(dict)

    info()
    info('Norm depth by gene, norm2, norm3...')
    for gene, norm_depth_by_sample in norm_depths_by_sample.items():
        info(str(gene), ending=' ', print_date=False)
        for gene_info in records:
            sample = gene_info.sample_name
            if sample not in norm_depth_by_sample:
                continue

            gene_norm_depth = norm_depth_by_sample[sample] * factors_by_gene[gene] + 0.1  # norm1b

            norm_depths_by_gene[gene][sample] = gene_norm_depth

            norm2[gene][sample] = math.log(gene_norm_depth / med_depth, 2) if med_depth else 0

            norm3[gene][sample] = math.log(gene_norm_depth / median_depth_by_sample[sample], 2) if \
                median_depth_by_sample[sample] else 0
    info(print_date=False)

    return _get_report_data(records, norm2, norm3, norm_depths_by_gene, norm_depths_by_sample)


def _get_med_norm_depths(mapped_reads_by_sample, norm_depths_by_sample):
    mean_depths_by_sample = defaultdict(list)

    for sample_name in mapped_reads_by_sample:
        for gn, v_by_sample in norm_depths_by_sample.items():
            if sample_name not in v_by_sample:
                err(sample_name + ' not in ' + str(v_by_sample))
                continue
            mean_depths_by_sample[sample_name].append(v_by_sample[sample_name])

    median_depth_by_sample = OrderedDict()
    for sample_name in mapped_reads_by_sample:
        median_depth_by_sample[sample_name] = median(mean_depths_by_sample[sample_name])

    return median_depth_by_sample


# mean_reads =  mean for all the samples mapped reads
# return factor_by_sample = sample mapped read/mean for all the samples mapped reads
def _get_factors_by_sample(mapped_reads_by_sample):
    factor_by_sample = dict()

    mean_reads = mean(mapped_reads_by_sample.values())
    for k, v in mapped_reads_by_sample.items():
        factor_by_sample[k] = mean_reads / v if v else 0

    return factor_by_sample


#med_depth = median for all gene
#return factor_by_gene = median gene for all samples/ med_depth
def _get_factors_by_gene(norm_depths_1, med_depth):
    # mean_depth_by_genes = defaultdict(list)
    # for gene_name, values in norm_depths_by_gene:
    #     mean_depth_by_genes[gene_name].append(rec.mean_depth)

    factors_by_gene = dict()
    for gene_name, values in norm_depths_1.items():
        med = median(values.values())
        factors_by_gene[gene_name] = med_depth / med if med else 0

    return factors_by_gene


#MeanDepth_Norm1
# gene -> { sample -> [] }
def _get_norm_depths_by_sample(mapped_reads_by_sample, records):
    norm_depths = defaultdict(dict)

    factor_by_sample = _get_factors_by_sample(mapped_reads_by_sample)

    for sample, factor in factor_by_sample.items():
        info('  ' + sample)
        for rec in [rec for rec in records if rec.sample_name == sample]:
            norm_depths[rec.gene_name][sample] = rec.mean_depth * factor

    return norm_depths


def _get_report_data(records, norm2, norm3, norm_depths_by_gene, norm_depths_by_seq_distr):
    header = ["Sample", "Gene", "Chr", "Start", "Stop", "Length", "MeanDepth", "MeanDepth_Norm1",
              "MeanDepth_Norm2", "log2Ratio_norm1", "log2Ratio_norm2"]
    report_data = []

    for rec in records:
        gene_name = rec.gene_name
        sample = rec.sample_name
        report_data.append(map(str,
           (sample, gene_name, rec.chrom, rec.start_position,
            rec.end_position, rec.size,
            '{0:.3f}'.format(rec.mean_depth, norm_depths_by_seq_distr[gene_name][sample]),
            '{0:.3f}'.format(norm_depths_by_gene[gene_name][sample]),
            '{0:.3f}'.format(norm2[gene_name][sample]),
            '{0:.3f}'.format(norm3[gene_name][sample]))))

    # report_data = sorted(report_data, key=itemgetter(2, 1))
    report_data = [header] + report_data

    return report_data


# sample_order = '''PC9_IRLR PC9_9291-L0B_1 PC9_IR-GM PC9_9291_6 PC9_9291-5 PC9_IR-4
# PC9_IR-6 PC9_9291-2 PC9 PC9_9291-3 gDNA'''.split()
#
# gene_order = 'PIK3CA KRAS NRAS HRAS BRAF EGFR'.split()

def _report_row_to_objects(gene_depth):
    inputs = OrderedDefaultDict(OrderedDict)
    # for gene_name in gene_order:
    #     for sample_name in sample_order:
    #         inputs[gene_name][sample_name] = None

    for read in gene_depth:
        rec = Record(*read)
        if 'Undetermined' not in rec.sample_name:
            inputs[rec.gene_name][rec.sample_name] = rec

    details_list = []
    for gene_dict in inputs.values():
        for rec in gene_dict.values():
            details_list.append(rec)

    return details_list


class Record():
    def __init__(self, sample_name=None, chrom=None, start_position=None, end_position=None, gene_name=None,
                 type="Gene-Amplicon", size=None, mean_depth=None):
        self.sample_name = sample_name
        self.chrom = chrom
        self.start_position = int(start_position)
        self.end_position = int(end_position)
        self.gene_name = gene_name
        self.type = type
        self.size = int(size)
        self.mean_depth = float(mean_depth)

    def __str__(self):
        values = [self.sample_name, self.chrom, self.start_position, self.end_position, self.gene_name, self.type, self.size,
                  self.mean_depth]
        return '"' + '\t'.join(map(str, values)) + '"'

    def __repr__(self):
        return repr(
            (self.gene_name, self.sample_name))

