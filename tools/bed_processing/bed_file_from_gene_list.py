#!/usr/bin/env python

from os.path import abspath, dirname, realpath, join

from site import addsitedir

project_dir = abspath(dirname(dirname(dirname(realpath(__file__)))))
addsitedir(join(project_dir))
addsitedir(join(project_dir, 'ext_modules'))
# import scripts.__check_python_version  # do not remove it: checking for python version and adding site dirs inside

import sys
from tools.bed_processing.make_exons import read_approved_genes, get_approved_gene_symbol


def main():
    if len(sys.argv) < 2:
        sys.stderr.write('Usage: ' + __file__ + ' genes.txt [regions.bed] [HGNC_gene_synonyms.txt] [--not-check] > regions_for_genes.txt')
        sys.stderr.write('    the script filters regions to have gene names in genes.txt')
        sys.exit(1)

    genes_list_fpath = sys.argv[1]
    exons_fpath = '/ngs/reference_data/genomes/Hsapiens/hg19/bed/Exons/Exons.bed'
    synonyms_fpath = '/ngs/reference_data/genomes/Hsapiens/common/HGNC_gene_synonyms.txt'
    if len(sys.argv) > 2: exons_fpath = sys.argv[2]
    if len(sys.argv) > 3: synonyms_fpath = sys.argv[3]
    not_check = len(sys.argv) > 4

    approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym = \
        read_approved_genes(synonyms_fpath)

    gnames = set()
    with open(genes_list_fpath) as f:
        for l in f:
            gname = l.strip()

            if not_check:
                gnames.add(gname)
            else:
                approved_gname, status = get_approved_gene_symbol(
                    approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym, gname)

                if not approved_gname:
                    gname2 = gname.split('.')[0]

                    if gname2 != gname:
                        log(gname + ' was not found with dot, trying without dot, as ' + gname2)

                        approved_gname, status2 = get_approved_gene_symbol(
                            approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym, gname2)

                if approved_gname:
                    if gname != approved_gname:
                        log(gname + ' found as ' + approved_gname)
                    else:
                        log(gname)
                    gnames.add(approved_gname)

                else:
                    log(gname + ' not found, skipping...')

    log('Found ' + str(len(gnames)) + ' uniq gene names in ' + genes_list_fpath)
    log()

    total_lines_written = 0
    genes_found_in_ref = set()

    with open(exons_fpath) as f:
        for l in f:
            if not l.startswith('#'):
                fs = l.strip().split('\t')
                if len(fs) < 4:
                    pass
                else:
                    g = fs[3]
                    if g in gnames:
                        total_lines_written += 1
                        genes_found_in_ref.add(g)
                        print '\t'.join(fs)

    log('Written ' + str(total_lines_written) + ' regions')
    log('Found in reference ' + str(len(genes_found_in_ref)) + ' genes out of ' + str(len(gnames)))
    log('Not found in reference: ' + ', '.join(gnames - genes_found_in_ref))


def log(msg=''):
    sys.stderr.write(msg + '\n')


if __name__ == '__main__':
    main()