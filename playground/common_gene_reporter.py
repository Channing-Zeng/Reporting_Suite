#!/usr/bin/env python

import sys
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are using %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from site import addsitedir
import os
from os.path import dirname, join, splitext, realpath, isdir
from optparse import OptionParser

this_py_fpath = splitext(__file__)[0] + '.py'
this_py_real_fpath = realpath(this_py_fpath)
project_dir = dirname(dirname(this_py_real_fpath))

addsitedir(join(project_dir))
addsitedir(join(project_dir, 'ext_modules'))
from source.file_utils import file_exists
from source.logger import info, err, warn, critical


def proc_args(argv):
    group1_name = 'Resistant'
    group2_name = 'Sensitive'

    description = 'This script find genes with mutations presented in (almost) all samples in one groups' \
                  'and (almost) not presented in another group' \
                  ' (default group names: Resistant vs Sensitive). Input is PASS.txt files from bcbio-postproc.'
    parser = OptionParser(description=description)
    parser.add_option('-n', '--num-samples-limit', dest='ns', default=1, type=int,
                      help='For each reported gene: max number of samples WITHOUT the gene in group1, '
                           'max number of samples WITH the gene in group2')

    (opts, args) = parser.parse_args(argv)
    
    if len(args) == 0:
        critical('No PASS.txt files provided to input.')
        
    variants_fpaths = [fpath for fpath in args if file_exists(fpath)]
    return opts, [group1_name, group2_name], variants_fpaths


def parse_variants(fpath):
    sample_column_name = 'Sample'
    gene_column_name = 'Gene'

    genes_per_sample = dict()
    with open(fpath) as f:
        header = f.readline().split('\t')
        if sample_column_name not in header:
            warn('"' + sample_column_name + '" is not found in ' + fpath + ' header, skipping this file!')
            return genes_per_sample
        else:
            sample_column_id = header.index(sample_column_name)
        if gene_column_name not in header:
            warn('"' + gene_column_name + '" is not found in ' + fpath + ' header, skipping this file!')
            return genes_per_sample
        else:
            gene_column_id = header.index(gene_column_name)
        for line in f:
            line = line.split('\t')
            sample = line[sample_column_id]
            gene = line[gene_column_id]
            if sample not in genes_per_sample:
                genes_per_sample[sample] = set()
            genes_per_sample[sample].add(gene)
    info('Found info for %d samples:' % len(genes_per_sample))
    for k, v in genes_per_sample.items():
        info('\t%s (%d unique genes)' % (k, len(v)))
    return genes_per_sample


def inverse_genes_per_sample_dict(genes_per_sample, min_ns=2):
    samples_per_gene = dict()
    all_genes = set()
    for gene_set in genes_per_sample.values():
        all_genes |= gene_set
    info('Total number of unique genes: %d' % len(all_genes))
    for gene in all_genes:
        samples_per_gene[gene] = []
        for sample, gene_set in genes_per_sample.items():
            if gene in gene_set:
                samples_per_gene[gene].append(sample)
        if len(samples_per_gene[gene]) < min_ns:
            del samples_per_gene[gene]
    info('Total number of unique genes presented in at least %d samples: %d' % (min_ns, len(samples_per_gene)))
    return samples_per_gene


def print_results(genes_per_group, limit):
    for cur_group in genes_per_group:
        sys.stderr.write('In %s group (limit = %d):\n' % (cur_group, limit))
        # unique_genes = set(genes_per_group[cur_group])
        # for group, genes in genes_per_group.items():
        #     if cur_group != group:
        #         unique_genes -= genes
        unique_genes = genes_per_group[cur_group]
        sys.stderr.write(', '.join(unique_genes) + '\n')


def __is_sample_in_group(sample, group):
    return group in sample


def find_suitable_genes_for_group(samples_per_gene, group, samples_per_group, limit):
    ns_this_group = len(samples_per_group[group])
    ns_other_group = sum(map(len, samples_per_group.values())) - ns_this_group
    info('Processing group ' + group)
    suitable_genes = []
    info('  %d total genes' % len(samples_per_gene))
    for gene, samples in samples_per_gene.items():
        if len([s for s in samples if __is_sample_in_group(s, group)]) >= ns_this_group - limit:
            suitable_genes.append(gene)
    info('  %d genes with mutations in >=%d samples of this group (out of %d)' %
         (len(suitable_genes), ns_this_group - limit, ns_this_group))
    #info(str(suitable_genes))
    for gene in list(suitable_genes):
        samples = samples_per_gene[gene]
        if len([s for s in samples if not __is_sample_in_group(s, group)]) > limit:
            suitable_genes.remove(gene)
    info('  %d remained after removing genes with mutations in >%d samples of another group (out of %d)' %
         (len(suitable_genes), limit, ns_other_group))
    #info(str(suitable_genes))
    info('  detailed:')
    for gene in suitable_genes:
        samples = samples_per_gene[gene]
        in_this_group = len([s for s in samples if __is_sample_in_group(s, group)])
        in_other_group = len([s for s in samples if not __is_sample_in_group(s, group)])
        info('     %s present in %d samples of this group and in %d samples of another group' %
             (gene, in_this_group, in_other_group))
    return suitable_genes


def main():
    opt, groups, variants_fpaths = proc_args(sys.argv[1:])
    if opt.ns < 0:
        opt.ns = 0

    genes_per_sample = dict()
    for fpath in variants_fpaths:
        genes_per_sample.update(parse_variants(fpath))
    info('')

    samples = genes_per_sample.keys()
    samples_per_group = {group: [sample for sample in samples if __is_sample_in_group(sample, group)]
                          for group in groups}
    info('Samples per groups: ')
    for group, samples in samples_per_group.items():
        info(group + ': ' + str(samples))

    samples_per_gene = inverse_genes_per_sample_dict(genes_per_sample,
                                                     min_ns=min(map(len, samples_per_group.values())) - opt.ns)

    genes_per_group = dict()
    for group in groups:
        genes_per_group[group] = find_suitable_genes_for_group(samples_per_gene, group, samples_per_group, opt.ns)

    # genes_per_group = dict()
    # info('')
    # for group in groups:
    #     for sample in [s for s in genes_per_sample.keys() if group in s]:
    #         if group not in genes_per_group:
    #             genes_per_group[group] = genes_per_sample[sample]
    #         else:
    #             genes_per_group[group] &= genes_per_sample[sample]
    #         info('%s sample added to %s group, total common genes: %d'
    #                     % (sample, group, len(genes_per_group[group])))
    # info('Intersected into %d groups: %s' % (len(groups), ', '.join(groups)))
    # for k, v in genes_per_group.items():
    #     info('\t%s (%d unique genes)' % (k, len(v)))

    print_results(genes_per_group, opt.ns)


if __name__ == '__main__':
    main()