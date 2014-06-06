#!/usr/bin/env python
from collections import defaultdict

from itertools import groupby, chain
import math
import numpy

_SAMPLE_MAPPED_READS = {'Sample1': 100, "Sample2": 200}  # cnt

_GENE_RECORDS = [
    ["Sample1", "chrM", 0, 1000, "DQ582201", "Gene-Amplicon", 1000, 1052.66],
    ["Sample1", "chrM", 0, 1000, "JB137816", "Gene-Amplicon", 1000, 1052.66],
    ["Sample1", "chrM", 2000, 4500, "TVAS5", "Gene-Amplicon", 1500, 984.02],
    ["Sample1", "chrM", 6000, 8000, "BC018860", "Gene-Amplicon", 2000, 0.00],
    ["Sample1", "chrM", 6000, 8000, "OK/SW-cl.16", "Gene-Amplicon", 2000, 0.00],
    ["Sample1", "chrM", 6000, 8000, "JA760602", "Gene-Amplicon", 2000, 0.00],
    ["Sample2", "chrM", 0, 1000, "DQ582201", "Gene-Amplicon", 1000, 1052.66],
    ["Sample2", "chrM", 0, 1000, "JB137816", "Gene-Amplicon", 1000, 1052.66]]


#sample2	chrM	2000	4500	TVAS5	Gene-Amplicon	1500	984.02
#sample2	chrM	6000	8000	BC018860	Gene-Amplicon	2000	0.00
#sample2	chrM	6000	8000	OK/SW-cl.16	Gene-Amplicon	2000	0.00
#sample2	chrM	6000	8000	JA760602	Gene-Amplicon	2000	0.00) }


def flatten(dict_of_dicts):
    try:
        for dic in dict_of_dicts.itervalues():
            for nested_v in flatten(dic):
                yield nested_v
    except (AttributeError, TypeError):
        for list_v in dict_of_dicts:
            yield list_v


# group by element
def group_by(lists, element_num):
    groups = dict()
    for l, v in groupby(sorted(lists, key=lambda x: x[element_num]), lambda x: x[element_num]):
        groups[l] = list(v)
    #res = [list(v) for l, v in groupby(sorted(lists, key=lambda x: x[element]), lambda x: x[element])]
    return groups


def mean_reads_for_all_samples(mapped_reads_by_sample):
    return numpy.mean(mapped_reads_by_sample.values())


#TODO chanf array index to get method ?
def median_by_column(lists, col):
    col = [l[col] for l in lists]
    return numpy.median(col)


def get_factors_by_sample(mapped_reads_by_sample):
    factor_by_sample = dict()
    meanreads = mean_reads_for_all_samples(mapped_reads_by_sample)
    for k, v in mapped_reads_by_sample.items():
        factor_by_sample[k] = meanreads/v
    return factor_by_sample


#TODO make an object from _GENE_MINDEPTH values
def get_norm_depths_by_seq_distr(mapped_reads_by_sample, records_by_sample):
    norm_depths = defaultdict(dict)  # gene -> { sample -> [] }
    depths_by_sample = defaultdict(list)
    factor_by_sample = get_factors_by_sample(mapped_reads_by_sample)

    for sample, factor in factor_by_sample.items():
        gene_infos = records_by_sample[sample]
        for gene_info in gene_infos:
            gene_name = gene_info[4]
            depth = gene_info[-1]
            norm_depths[gene_name][sample] = depth * factor
            depths_by_sample[sample].append(depth)

    return norm_depths


def get_factors_by_gene(gene_records, med_depth):
    values_by_gene = group_by(gene_records, 4)

    factors_by_gene = dict()

    for gene, values in values_by_gene.iteritems():
        median = median_by_column(values, 7)
        factors_by_gene[gene] = med_depth / median if median else 0

    return factors_by_gene


def main():
    records_by_sample = group_by(_GENE_RECORDS, 0)

    norm_depths_by_seq_distr = get_norm_depths_by_seq_distr(_SAMPLE_MAPPED_READS, records_by_sample)

    med_depth = numpy.median([rec[7] for rec in _GENE_RECORDS])

    factors_by_gene = get_factors_by_gene(_GENE_RECORDS, med_depth)

    norm_depths_by_gene = defaultdict(dict)

    norm2 = defaultdict(dict)

    norm3 = defaultdict(dict)

    for gene, norm_depth_by_sample in norm_depths_by_seq_distr.items():
        for sample in records_by_sample.keys():
            if sample not in norm_depth_by_sample:
                continue

            gene_norm_depth = norm_depth_by_sample[sample] * factors_by_gene[gene] + 0.1

            norm_depths_by_gene[gene][sample] = gene_norm_depth

            norm2[gene][sample] = math.log(gene_norm_depth / med_depth, 2)

            median_depth = numpy.median([rec[7] for rec in records_by_sample[sample]])

            norm3[gene][sample] = math.log(gene_norm_depth / median_depth, 2)




if __name__ == '__main__':
    main()




















