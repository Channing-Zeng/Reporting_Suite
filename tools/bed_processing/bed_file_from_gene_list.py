#!/usr/bin/env python
import bcbio_postproc  # do not remove it: checking for python version and adding site dirs inside

import sys
from source.logger import err
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

    not_found_in_HGNC_genes = list()
    all_genes = set()
    gnames = set()
    duplicated_genes = list()
    doublicated_after_approving = list()
    with open(genes_list_fpath) as f:
        for i, l in enumerate(f):
            gname = l.strip()
            sys.stderr.write(str(i + 1) + ' ' + gname + ' ')

            if gname in all_genes:
                duplicated_genes.append(gname)
                sys.stderr.write(' duplicated gene\n')
                continue
            else:
                all_genes.add(gname)

            gnames.add(gname)
            sys.stderr.write('\n')

    log('Found ' + str(len(gnames)) + ' uniq gene names in ' + genes_list_fpath)
    log()

    genes_found_in_ref, lines = exons_for_gene_list(exons_fpath, gnames)
    genes_not_found_in_ref = gnames - genes_found_in_ref
    log('Written ' + str(len(lines)) + ' regions')
    log('Total genes: ' + str(len(all_genes)) + ((', ' + str(len(duplicated_genes)) +
        ' duplicated: ' + ', '.join(duplicated_genes)) if duplicated_genes else ''))
    # log(str(len(doublicated_after_approving)) + ' duplicated after approval: ' + ', '.join(doublicated_after_approving))
    # log('Not found in HGNC ' + str(len(not_found_in_HGNC_genes)) + ' genes' +
    #     (': ' + ', '.join(not_found_in_HGNC_genes) if not_found_in_HGNC_genes else ''))
    log('Found in reference ' + str(len(genes_found_in_ref)) + ' genes out of ' + str(len(gnames)))
    log('Not found in reference: ' + ', '.join(genes_not_found_in_ref))

    corrected_genes_not_found_in_ref = set()
    if genes_not_found_in_ref:
        log()
        log('Correcting genes not found in Exons DB')
        for gname in genes_not_found_in_ref:
            sys.stderr.write('\n')
            sys.stderr.write(gname + '\n')

            approved_gname, status = get_approved_gene_symbol(
                approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym, gname,
                indent='    ')

            if not approved_gname:
                if '.' in gname:
                    gname2 = gname.split('.')[0]
                    sys.stderr.write('    Not found with dot, trying without dot, as ' + gname2 + '...\n')
                    approved_gname, status2 = get_approved_gene_symbol(
                        approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym, gname2,
                        indent='    ')

            if approved_gname:
                if gname != approved_gname:
                    sys.stderr.write('    Found as ' + approved_gname)
                else:
                    sys.stderr.write(gname )
                if approved_gname in gnames:
                    sys.stderr.write('    Approved version is already met: ' + approved_gname)
                    doublicated_after_approving.append(gname + '->' + approved_gname)

                corrected_genes_not_found_in_ref.add(approved_gname)

            else:
                sys.stderr.write('    Not found, skipping...')
                not_found_in_HGNC_genes.append(gname)
            sys.stderr.write('\n')

    genes_found_in_ref, lines_2 = exons_for_gene_list(exons_fpath, corrected_genes_not_found_in_ref)

    log('After correction, written ' + str(len(lines_2)) + ' regions')
    # log(str(len(doublicated_after_approving)) + ' duplicated after approval: ' + ', '.join(doublicated_after_approving))
    # log('Not found in HGNC ' + str(len(not_found_in_HGNC_genes)) + ' genes' +
    #     (': ' + ', '.join(not_found_in_HGNC_genes) if not_found_in_HGNC_genes else ''))
    log('Found in reference ' + str(len(genes_found_in_ref)) + ' genes out of ' + str(len(corrected_genes_not_found_in_ref)))
    log('Not found in reference: ' + ', '.join(corrected_genes_not_found_in_ref - genes_found_in_ref))

    lines.extend(lines_2)
    for l in lines:
        print l


def exons_for_gene_list(exons_fpath, gnames):
    genes_found_in_ref = set()
    lines = []

    with open(exons_fpath) as f:
        for l in f:
            if not l.startswith('#'):
                fs = l.strip().split('\t')
                if len(fs) < 4:
                    pass
                else:
                    gname = fs[3]
                    if gname not in gnames:
                        pass
                        # genes_not_found_in_ref.add(gname)
                    else:
                        genes_found_in_ref.add(gname)
                        lines.append('\t'.join(fs))

    return genes_found_in_ref, lines


def log(msg=''):
    sys.stderr.write(msg + '\n')


if __name__ == '__main__':
    main()