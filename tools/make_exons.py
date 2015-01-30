#!/usr/bin/env python

import sys

if len(sys.argv) < 1:
    sys.stderr.write('The script writes all exons for all known UCSC genes, with associated gene symbols.\n')
    sys.stderr.write('When the gene name is found in HGNC, it get replaced with an approved name.\n')
    sys.stderr.write('If the gene is not charactirized (like LOC729737), this symbol is just kept as is.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('    ' + __file__ + ' HGNC_gene_synonyms.txt [file_to_write_not_approved_genes.txt] < UCSC_knownGene.txt > UCSC_exons.bed\n')
    sys.stderr.write('\n')
    sys.stderr.write('    where HGNC_gene_synonyms.txt (from http://www.genenames.org/cgi-bin/download) is:\n')
    sys.stderr.write('      #Approved Symbol  Previous Symbols                    Synonyms                          UCSC ID(supplied by UCSC)\n')
    sys.stderr.write('      OR7E26P           OR7E67P, OR7E69P, OR7E70P, OR7E68P  OR1-51, OR1-72, OR1-73, OR912-95\n')
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
not_approved_fpath = sys.argv[2] if len(sys.argv) > 2 else None

approved_gname_by_gname = dict()
with open(synonyms_fpath) as syn:
    i = 1
    for l in syn:
        if l and not l.startswith('#'):
            ts = l[:-1].split('\t')
            approved_gn, prev_names, synonyms, ucsc_ids = l.split('\t')
            approved_gname_by_gname[approved_gn] = approved_gn
            for gn in prev_names.split(', ') + synonyms.split(', '):
                if gn:
                    approved_gname_by_gname[gn] = approved_gn
        if i % 10000 == 0:
            sys.stderr.write('Processed ' + str(i) + ' lines from ' + synonyms_fpath + '\n')
        i += 1
    sys.stderr.write('Processed ' + str(i) + ' lines from ' + synonyms_fpath + '\n')
    sys.stderr.write('\n')

# inp = sys.stdin
# out = sys.stdout
i = 1

not_approved_f = open(not_approved_fpath, 'w') if not_approved_fpath else None

with sys.stdin as inp, sys.stdout as out:
    for l in inp:
        if l and not l.startswith('#'):
            g_id, chrom, strand, txStart, txEnd, exonCount, exonStarts, exonEnds, geneSymbol = l[:-1].split('\t')
            approved_gene_symbol = approved_gname_by_gname.get(geneSymbol)
            if approved_gene_symbol:
                for j, s, e in zip(range(int(exonCount)), [e for e in exonStarts.split(',') if e], [e for e in exonEnds.split(',') if e]):
                    out.write('\t'.join([chrom, s, e, approved_gene_symbol, str(j), strand]) + '\n')
            else:
                if not_approved_f:
                    not_approved_f.write(geneSymbol + '\n')
        if i % 10000 == 0:
            sys.stderr.write('Processed ' + str(i) + ' lines from stdin\n')
        i += 1
    sys.stderr.write('Processed ' + str(i) + ' lines from stdin.\n')

if not_approved_f:
    not_approved_f.close()
# inp.close()
# out.close()