#!/usr/bin/env python

from os.path import abspath, dirname, realpath, join, exists
from site import addsitedir
project_dir = abspath(dirname(dirname(realpath(__file__))))
addsitedir(join(project_dir))
addsitedir(join(project_dir, 'ext_modules'))
import sub_scripts.__check_python_version  # do not remove it: checking for python version and adding site dirs inside

import sys
from tools.make_exons import read_approved_genes, get_approved_gene_symbol


def main():
    if len(sys.argv) < 2:
        sys.stderr.write('Usage: ' + __file__ + ' genes.txt [HGNC_gene_synonyms.txt] [regions.bed] > regions_for_genes.txt')
        sys.stderr.write('    the script filters regions to have gene names in genes.txt')
        sys.exit(1)

    genes_list_fpath = sys.argv[1]
    exons_fpath = '/ngs/reference_data/genomes/Hsapiens/hg19/bed/Exons/Exons.bed'
    synonyms_fpath = '/ngs/reference_data/genomes/Hsapiens/hg19/bed/Exons/HGNC_gene_synonyms.txt'
    if len(sys.argv) > 2: synonyms_fpath = sys.argv[2]
    if len(sys.argv) > 3: exons_fpath = sys.argv[3]

    approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym = \
        read_approved_genes(synonyms_fpath)

    genes = set()
    with open(genes_list_fpath) as f:
        for l in f:
            gname = l.strip()

            approved_gname, status = get_approved_gene_symbol(
                approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym, gname)

            if not approved_gname:
                gname2 = gname.split('.')[0]

                if gname2 != gname:
                    log(gname + ' was not found with dot, trying without dot, as ' + gname2)

                    approved_gname, status2 = get_approved_gene_symbol(
                        approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym, gname2)

            if approved_gname:
                log(gname + ' found' + ((' as ' + approved_gname) if approved_gname else ''))
                genes.add(approved_gname)

            else:
                log(gname + ' not found, skipping...')

    log('Found ' + str(len(genes)) + ' uniq gene names in ' + genes_list_fpath)
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
                    if g in genes:
                        total_lines_written += 1
                        genes_found_in_ref.add(g)
                        print '\t'.join(fs)

    log('Written ' + str(total_lines_written) + ' regions')
    log('Found in reference ' + str(len(genes_found_in_ref)) + ' genes out of ' + str(len(genes)))


def log(msg=''):
    sys.stderr.write(msg + '\n')


if __name__ == '__main__':
    main()