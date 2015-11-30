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
from source import logger


def proc_args(argv):
    group1_name = 'Resistant'
    group2_name = 'Sensitive'

    description = 'This script generates gene-level report of differences between two groups (Resistant vs Sensitive).'
    parser = OptionParser(description=description)
    #parser.add_option('-o', dest='output_dir', metavar='DIR', default=join(os.getcwd(), 'targetqc'))

    (opts, args) = parser.parse_args()
    variants_fpaths = [fpath for fpath in args if file_exists(fpath)]
    return opts, [group1_name, group2_name], variants_fpaths


def parse_variants(fpath):
    sample_column_name = 'Sample'
    gene_column_name = 'Gene'

    genes_per_sample = dict()
    with open(fpath) as f:
        header = f.readline().split('\t')
        if sample_column_name not in header:
            logger.warn('"' + sample_column_name + '" is not found in ' + fpath + ' header, skipping this file!')
            return genes_per_sample
        else:
            sample_column_id = header.index(sample_column_name)
        if gene_column_name not in header:
            logger.warn('"' + gene_column_name + '" is not found in ' + fpath + ' header, skipping this file!')
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
    logger.info('Found info for %d samples:' % len(genes_per_sample))
    for k, v in genes_per_sample.items():
        logger.info('\t%s (%d unique genes)' % (k, len(v)))
    return genes_per_sample


def print_results(genes_per_group):
    for cur_group in genes_per_group:
        sys.stderr.write('In %s group only:\n' % cur_group)
        unique_genes = set(genes_per_group[cur_group])
        for group, genes in genes_per_group.items():
            if cur_group != group:
                unique_genes -= genes
        sys.stderr.write(', '.join(unique_genes) + '\n')


def main():
    opt, groups, variants_fpaths = proc_args(sys.argv)
    genes_per_sample = dict()
    for fpath in variants_fpaths:
        genes_per_sample.update(parse_variants(fpath))
    genes_per_group = dict()
    logger.info('')
    for group in groups:
        for sample in [s for s in genes_per_sample.keys() if group in s]:
            if group not in genes_per_group:
                genes_per_group[group] = genes_per_sample[sample]
            else:
                genes_per_group[group] &= genes_per_sample[sample]
            logger.info('%s sample added to %s group, total common genes: %d'
                        % (sample, group, len(genes_per_group[group])))
    logger.info('Intersected into %d groups: %s' % (len(groups), ', '.join(groups)))
    for k, v in genes_per_group.items():
        logger.info('\t%s (%d unique genes)' % (k, len(v)))
    print_results(genes_per_group)


if __name__ == '__main__':
    main()