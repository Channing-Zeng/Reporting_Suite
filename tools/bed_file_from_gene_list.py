#!/usr/bin/env python

import sys
from tools.make_exons import read_approved_genes, get_approved_gene_symbol
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

genes = []
with open(genes_list_fpath) as f:
    for l in f:
        gname = l.strip()

        approved_gname, status = get_approved_gene_symbol(
            approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym, gname)

        if not approved_gname:
            gname2 = gname.split('.')[0]
            if gname2 != gname:
                approved_gname, status2 = get_approved_gene_symbol(
                    approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym, gname2)

        if approved_gname:
            genes.append(approved_gname)


with open(exons_fpath) as f:
    for l in f:
        if not l.startswith('#'):
            fs = l.strip().split('\t')
            if len(fs) < 4:
                pass
            else:
                g = fs[3]
                if g in genes:
                    print '\t'.join(fs)
