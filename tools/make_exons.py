#!/usr/bin/env python

from collections import defaultdict
import sys


class ApprovedGene:
    def __init__(self, gname, prev_names, synonyms, chrom, ucsc_id=None):
        self.gname = gname
        self.prev_names = prev_names
        self.synonyms = synonyms
        self.chrom = chrom
        self.ucsc_id = ucsc_id


def _read_approved_genes(synonyms_fpath):
    approved_gene_by_name = dict()
    approved_gnames_by_prev_gname = defaultdict(list)
    approved_gnames_by_synonym = defaultdict(list)

    with open(synonyms_fpath) as syn:
        i = 1
        for l in syn:
            if l and not l.startswith('#'):
                approved_gn, prev_names, synonyms, chrom, ucsc_id = l[:-1].split('\t')
                if chrom:
                    chrom = 'chr' + chrom.split('p')[0]
                approved_gene = ApprovedGene(approved_gn, prev_names, synonyms, chrom, ucsc_id)
                approved_gene_by_name[approved_gn] = approved_gene

                for gn in prev_names.split(', '):
                    if gn:
                        approved_gnames_by_prev_gname[gn].append(approved_gene)

                for gn in synonyms.split(', '):
                    if gn:
                        approved_gnames_by_synonym[gn].append(approved_gene)

            if i % 10000 == 0:
                sys.stderr.write('Processed ' + str(i) + ' lines from ' + synonyms_fpath + '\n')
            i += 1
        sys.stderr.write('Processed ' + str(i) + ' lines from ' + synonyms_fpath + '\n')
        sys.stderr.write('\n')

    return approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym


def _check_gene_symbol(approved_gene, gene_symbol, ucsc_g_id, chrom):
    if ucsc_g_id != approved_gene.ucsc_id:
        sys.stderr.write(gene_symbol + ': HGNC ucsc_id = ' + approved_gene.ucsc_id + ', UCSC id = ' + ucsc_g_id + '\n')
    if chrom != approved_gene.chrom:
        sys.stderr.write(gene_symbol + ': HGNC chrom = ' + approved_gene.chrom + ', UCSC chrom = ' + chrom + '\n')


def _get_approved_gene_symbol(approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym,
                              gene_symbol, ucsc_g_id, chrom):
    if gene_symbol in approved_gene_by_name:
        _check_gene_symbol(approved_gene_by_name[gene_symbol], gene_symbol, ucsc_g_id, chrom)
        return approved_gene_by_name[gene_symbol]

    def _get_approved_by_prev_gname(approved_genes, kind):
        if approved_genes:
            if len(approved_genes) == 1:
                _check_gene_symbol(approved_genes[0], gene_symbol, ucsc_g_id, chrom)
                return approved_genes[0]

            if len(approved_genes) > 1:
                approved_genes_same_ucsc = [g for g in approved_genes if g.ucsc_id == ucsc_g_id]
                if len(approved_genes_same_ucsc) == 1:
                    _check_gene_symbol(approved_genes_same_ucsc[0], gene_symbol, ucsc_g_id, chrom)
                    return approved_genes_same_ucsc[0]
                if len(approved_genes_same_ucsc) > 1:
                    sys.stderr.write('Err: multiple approved gene names for ' + gene_symbol + ' (as ' + kind + ') with ucsc_id ' + ucsc_g_id)
                    return 'AMBIGOUS'

                if len(approved_genes_same_ucsc) == 0:
                    approved_genes_same_chrom = [g for g in approved_genes if g.chrom == chrom]
                    if len(approved_genes_same_chrom) == 1:
                        g = approved_genes_same_chrom[0]
                        sys.stderr.write('Only ' + g.gname + ' for ' + gene_symbol + ' (as ' + kind + ') has the same chrom ' + chrom + ', picking it')
                        _check_gene_symbol(g, gene_symbol, ucsc_g_id, chrom)
                        return g
                    if len(approved_genes_same_chrom) > 1:
                        sys.stderr.write('Err: multiple approved gene names for ' + gene_symbol + ' (as ' + kind + ') with chrom ' + chrom)
                        return 'AMBIGOUS'
                    if len(approved_genes_same_chrom) == 0:
                        sys.stderr.write('Err: no approved gene names for ' + gene_symbol + ' (as ' + kind + ') with same chrom ' + chrom)
                        return 'NOT FOUND'

    res = _get_approved_by_prev_gname(approved_gnames_by_prev_gname.get(gene_symbol), 'prev')
    if res == 'AMBIGOUS':
        return None
    elif res == 'NOT FOUND':
        res = _get_approved_by_prev_gname(approved_gnames_by_synonym.get(gene_symbol), 'synonym')
        if res in ['AMBIGOUS', 'NOT FOUND']:
            return None
        else:
            return res
    else:
        return res


def main():
    if len(sys.argv) < 1:
        sys.stderr.write('The script writes all exons for all known UCSC genes, with associated gene symbols.\n')
        sys.stderr.write('When the gene name is found in HGNC, it get replaced with an approved name.\n')
        sys.stderr.write('If the gene is not charactirized (like LOC729737), this symbol is just kept as is.\n')
        sys.stderr.write('\n')
        sys.stderr.write('Usage:\n')
        sys.stderr.write('    ' + __file__ + ' HGNC_gene_synonyms.txt [file_to_write_not_approved_genes.txt] < UCSC_knownGene.txt > UCSC_exons.bed\n')
        sys.stderr.write('\n')
        sys.stderr.write('    where HGNC_gene_synonyms.txt (from http://www.genenames.org/cgi-bin/download) is:\n')
        sys.stderr.write('      #Approved Symbol  Previous Symbols                    Synonyms                          Chromosome   UCSC ID(supplied by UCSC)\n')
        sys.stderr.write('      OR7E26P           OR7E67P, OR7E69P, OR7E70P, OR7E68P  OR1-51, OR1-72, OR1-73, OR912-95  19q13.43	 uc002qsg.3\n')
        sys.stderr.write('      ...\n')
        sys.stderr.write('\n')
        sys.stderr.write('    and UCSC_knownGene.txt (from http://genome.ucsc.edu/cgi-bin/hgTables) is:\n')
        sys.stderr.write('      #hg19.knownGene.name  hg19.knownGene.chrom  hg19.knownGene.strand  hg19.knownGene.txStart  hg19.knownGene.txEnd  hg19.knownGene.exonCount  hg19.knownGene.exonStarts  hg19.knownGene.exonEnds  hg19.kgXref.geneSymbol\n')
        sys.stderr.write('      uc001aaa.3	          chr1	                +	                   11873                   14409                 3                         11873,12612,13220,	      12227,12721,14409,	   DDX11L1\n')
        sys.stderr.write('      ...\n')
        sys.stderr.write('\n')
        sys.stderr.write('See more info in http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Making+the+full+list+of+UCSC+exons+with+approved+HUGO+gene+symbols\n')
        sys.exit(1)

    synonyms_fpath = sys.argv[1]
    approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym = _read_approved_genes(synonyms_fpath)

    not_approved_fpath = sys.argv[2] if len(sys.argv) > 2 else None
    not_approved_gene_names = set()

    i = 1
    with sys.stdin as inp, sys.stdout as out:
        for l in inp:
            if l and not l.startswith('#'):
                g_id, chrom, strand, txStart, txEnd, exonCount, exonStarts, exonEnds, geneSymbol = l[:-1].split('\t')

                approved_gene_symbol = _get_approved_gene_symbol(approved_gene_by_name, approved_gnames_by_prev_gname,
                                                                 approved_gnames_by_synonym, geneSymbol, g_id, chrom)
                if approved_gene_symbol:
                    for j, s, e in zip(range(int(exonCount)), [e for e in exonStarts.split(',') if e], [e for e in exonEnds.split(',') if e]):
                        out.write('\t'.join([chrom, s, e, approved_gene_symbol.gname, str(j), strand]) + '\n')
                else:
                    not_approved_gene_names.add(geneSymbol)

            if i % 10000 == 0:
                sys.stderr.write('Processed ' + str(i) + ' lines from stdin\n')
            i += 1
        sys.stderr.write('Processed ' + str(i) + ' lines from stdin.\n')


    sys.stderr.write('Not approved ' + str(len(not_approved_gene_names)) + ' genes.\n')
    if not_approved_fpath:
        with open(not_approved_fpath, 'w') as f:
            for not_approved_gn in not_approved_gene_names:
                f.write(not_approved_gn + '\n')
        sys.stderr.write('Saved not approved to ' + not_approved_fpath + '\n')


if __name__ == '__main__':
    main()