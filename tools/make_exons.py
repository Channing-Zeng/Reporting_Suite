#!/usr/bin/env python

from os.path import abspath, dirname, realpath, join
from site import addsitedir
project_dir = abspath(dirname(dirname(realpath(__file__))))
addsitedir(join(project_dir))
addsitedir(join(project_dir, 'ext_modules'))
import sub_scripts.__check_python_version  # do not remove it: checking for python version and adding site dirs inside

from collections import defaultdict
import sys
from source.logger import err, is_local


class ApprovedGene:
    def __init__(self, gname, prev_names, synonyms, chrom, ucsc_id=None, ensembl_id=None):
        self.gname = gname
        self.prev_names = prev_names
        self.synonyms = synonyms
        self.chrom = chrom

        self.db_id = ensembl_id


def parse_hgnc_chrom(chrom):
    if chrom in ['reserved', 'c10_B']:
        return None

    CHROMS = ['Y', 'X', 'mitochondria']
    for i in range(22, 0, -1):
        CHROMS.append(str(i))

    for c in CHROMS:
        if chrom.startswith(c):
            if c == 'mitochondria':
                return 'chrM'
            return 'chr' + c

    sys.stderr.write('Cannot parse chromosome ' + chrom + '\n')
    return None


def parse_ensembl_chrom(chrom):
    CHROMS = ['Y', 'X', 'MT']
    for i in range(22, 0, -1):
        CHROMS.append(str(i))

    for c in CHROMS:
        if chrom.startswith(c):
            if c == 'MT':
                return 'chrM'
            return 'chr' + c

    return None


def _read_approved_genes(synonyms_fpath):
    approved_gene_by_name = dict()
    approved_gnames_by_prev_gname = defaultdict(list)
    approved_gnames_by_synonym = defaultdict(list)

    with open(synonyms_fpath) as f:
        i = 0
        for l in f:
            if l and not l.startswith('#'):
                approved_gn, prev_names, synonyms, hgnc_chrom, ensembl_id, ucsc_id = l[:-1].split('\t')
                if hgnc_chrom:
                    hgnc_chrom = parse_hgnc_chrom(hgnc_chrom)

                approved_gene = ApprovedGene(approved_gn, prev_names, synonyms, hgnc_chrom, ucsc_id, ensembl_id)
                approved_gene_by_name[approved_gn] = approved_gene

                for gn in prev_names.split(', '):
                    if gn:
                        approved_gnames_by_prev_gname[gn].append(approved_gene)

                for gn in synonyms.split(', '):
                    if gn:
                        approved_gnames_by_synonym[gn].append(approved_gene)
            i += 1
        sys.stderr.write('Processed ' + str(i) + ' lines from ' + synonyms_fpath + '\n')
        sys.stderr.write('\n')

    return approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym


def _check_gene_symbol(approved_gene, gene_symbol, db_id, chrom):
    if db_id != approved_gene.db_id:
        sys.stderr.write('Discordant db ids for ' + gene_symbol + ': db id = ' + str(db_id) + ', in HGNC it is ' + str(approved_gene.db_id) + '\n')
    # else:
        # sys.stderr.write('Accordant db ids for ' + gene_symbol + ': db id = ' + str(db_id) + ', in HGNC it is ' + str(approved_gene.db_id) + '\n')

    if chrom != approved_gene.chrom:
        sys.stderr.write('Discordant chroms for ' + gene_symbol + ': chrom = ' + chrom + ', in HGNC chrom is ' + approved_gene.chrom + '\n')
        return None

    return approved_gene


def _get_approved_gene_symbol(approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym,
                              gene_symbol, db_id, db_chrom):
    if gene_symbol in approved_gene_by_name:
        if _check_gene_symbol(approved_gene_by_name[gene_symbol], gene_symbol, db_id, db_chrom):
            return approved_gene_by_name[gene_symbol].gname, None

    sys.stderr.write('Gene name ' + gene_symbol + ' is not approved, searching for an approved verion.\n')

    def _get_approved_genes_by_kind(approved_genes, kind):
        if not approved_genes:
            return 'NOT FOUND'

        if len(approved_genes) > 1:
            approved_genes_same_ucsc = [g for g in approved_genes if g.db_id == db_id]

            if len(approved_genes_same_ucsc) > 1:
                sys.stderr.write('  Error: multiple approved gene names for ' + gene_symbol + ' (as ' + kind + ') with ucsc_id ' + db_id + ': ' + ', '.join(g.gname for g in approved_genes_same_ucsc) + '\n')
                return 'AMBIGOUS'

            if len(approved_genes_same_ucsc) == 1:
                if _check_gene_symbol(approved_genes_same_ucsc[0], gene_symbol, db_id, db_chrom):
                    sys.stderr.write('  Found approved gene for ' + gene_symbol + ' (as ' + kind + ') with ucsc_id ' + db_id + '\n')
                    return approved_genes_same_ucsc[0].gname

            # Ok, no genes with same ucsc id, or not the same chromosome for them.

            approved_genes_same_chrom = [g for g in approved_genes if g.chrom == db_chrom]

            if len(approved_genes_same_chrom) > 1:
                sys.stderr.write('  Error: multiple approved gene names for ' + gene_symbol + ' (as ' + kind + ') with chrom ' + db_chrom + ', '.join(g.gname for g in approved_genes_same_ucsc) + '\n')
                return 'AMBIGOUS'

            if len(approved_genes_same_chrom) == 1:
                g = approved_genes_same_chrom[0]
                sys.stderr.write('  Only ' + g.gname + ' for ' + gene_symbol + ' (as ' + kind + ') has the same chrom ' + db_chrom + ', picking it\n')
                if _check_gene_symbol(g, gene_symbol, db_id, db_chrom):
                    return g.gname
                else:
                    return 'NOT FOUND'

            if len(approved_genes_same_chrom) == 0:
                sys.stderr.write('  Error: no approved gene names for ' + gene_symbol + ' (as ' + kind + ') with same chrom ' + db_chrom + '\n')
                return 'NOT FOUND'

        if len(approved_genes) == 1:
            if _check_gene_symbol(approved_genes[0], gene_symbol, db_id, db_chrom):
                sys.stderr.write('  Found approved gene for ' + gene_symbol + ' (as ' + kind + ')\n')
                return approved_genes[0].gname

        return 'NOT FOUND'

    res = _get_approved_genes_by_kind(approved_gnames_by_prev_gname.get(gene_symbol), 'prev')
    if res == 'AMBIGOUS':
        return None, 'AMBIGOUS\tAS PREV'
    elif res == 'NOT FOUND':
        res = _get_approved_genes_by_kind(approved_gnames_by_synonym.get(gene_symbol), 'synonym')
        if res == 'AMBIGOUS':
            return None, res + '\tAS SYNONYM'
        if res == 'NOT FOUND':
            return None, res
        else:
            sys.stderr.write('  Finally found approved gene for ' + gene_symbol + ' (as synonym):' + res + '\n')
            return res, None
    else:
        sys.stderr.write('  Finally found approved gene for ' + gene_symbol + ' (as prev):' + res + '\n')
        return res, None


def _proc_ucsc(inp, out, approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym):
    not_approved_gene_names = list()

    for l in inp:
        if l and not l.startswith('#'):
            ucsc_id, ucsc_chrom, strand, txStart, txEnd, exonCount, exonStarts, exonEnds, geneSymbol = l[:-1].split('\t')

            approved_gene_symbol, status = _get_approved_gene_symbol(
                approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym,
                geneSymbol, ucsc_id, ucsc_chrom)

            if approved_gene_symbol:
                for j, s, e in zip(range(int(exonCount)),
                   [e for e in exonStarts.split(',') if e], [
                    e for e in exonEnds.split(',') if e]):
                    out.write('\t'.join([ucsc_chrom, s, e, approved_gene_symbol, '.', strand, 'gene']) + '\n')
            else:
                not_approved_gene_names.append(geneSymbol + '\t' + status)

    return not_approved_gene_names


class Gene:
    def __init__(self, name, chrom, start, end, strand, biotype, db_id):
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.biotype = biotype
        self.db_id = db_id

        self.exons = []


class Exon:
    def __init__(self, gene, start, end, biotype):
        self.gene = gene
        self.start = start
        self.end = end
        self.biotype = biotype


def _rm_quotes(l):
    return l[1:-1]


def _proc_ensembl(inp, out, approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym):
    gene_by_name = dict()
    genes = []

    zero_size_gene_lines = []
    zero_size_cds_lines = []

    sys.stderr.write('Parsing Ensembl input...\n')
    i = 0
    for l in inp:
        if l and not l.startswith('#'):
            chrom, biotype, feature, start, end, _, strand, _, props_line = l[:-1].split('\t')

            if is_local:
                if chrom != '21':
                    continue

            if feature not in ['gene', 'CDS']:
                continue

            start, end = str(int(start) - 1), end

            if int(end) <= int(start):
                sys.stderr.write('Error: start > end: ' + l + '\n')
                continue

            chrom = parse_ensembl_chrom(chrom)
            if not chrom:
                continue

            _prop_dict = dict(t.strip().split(' ') for t in props_line.split(';') if t.strip())
            gene_symbol = _rm_quotes(_prop_dict['gene_name'])
            db_id = _rm_quotes(_prop_dict['gene_id'])
            gene_biotype = _rm_quotes(_prop_dict['gene_biotype'])

            if feature == 'gene':
                assert gene_biotype == biotype, 'Gene: gene_biotype "' + gene_biotype + '" do not match biotype "' + biotype + '" for ' + gene_symbol
                assert gene_symbol not in genes, 'Error: duplicated gene ' + gene_symbol
                gene = Gene(gene_symbol, chrom, start, end, strand, biotype, db_id)
                gene_by_name[gene_symbol] = gene
                genes.append(gene)

            elif feature == 'CDS':
                assert gene_symbol in gene_by_name, 'Error: CDS record before gene record ' + gene_symbol
                gene = gene_by_name[gene_symbol]
                assert gene_biotype == gene.biotype, 'CDS: gene_biotype "' + gene_biotype + '" do not match biotype "' + gene.biotype + '" for ' + gene_symbol
                exon = Exon(gene, start, end, biotype)
                gene.exons.append(exon)

        i += 1
        if i % 1000 == 0:
            sys.stderr.write('processed ' + str(i) + ' lines, ' + str(len(genes)) + ' genes found\n')
            sys.stderr.flush()

    sys.stderr.write('\n')
    sys.stderr.write('Processed ' + str(i) + ' lines, ' + str(len(genes)) + ' genes found\n')
    sys.stderr.write('Found ' + str(len(zero_size_cds_lines)) + ' zero-size CDS, ' + str(len(zero_size_gene_lines)) + ' zero-size genes.\n')
    sys.stderr.write('\n')

    not_approved_gene_names = dict()

    exons_num = 0
    j = 0
    for g in genes:
        if len(g.exons) == 0:
            continue

        approved_gname, status = _get_approved_gene_symbol(
            approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym,
            g.name, g.db_id, g.chrom)

        if not approved_gname:
            gname2 = g.name.split('.')[0]
            if gname2 != g.name:
                approved_gname, status2 = _get_approved_gene_symbol(
                    approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym,
                    gname2, g.db_id, g.chrom)

            sys.stderr.write('  Not found approved gene for ' + g.name + ' (' + status + ')')
            if gname2 != g.name:
                if not approved_gname:
                    sys.stderr.write(' and ' + gname2 + ' (' + status2 + ')')
                else:
                    sys.stderr.write(', but approved for ' + gname2 + ' as ' + approved_gname)
            sys.stderr.write('\n')

        if approved_gname:
            out.write('\t'.join([g.chrom, g.start, g.end, approved_gname, '.', g.strand, 'gene', g.biotype]) + '\n')
            for e in g.exons:
                exons_num += 1
                out.write('\t'.join([g.chrom, e.start, e.end, approved_gname, '.', g.strand, 'CDS', e.biotype]) + '\n')
        else:
            if g.name not in not_approved_gene_names:
                not_approved_gene_names[g.name] = status

        j += 1
        if j % 1000 == 0:
            sys.stderr.write('processed ' + str(j / 1000) + 'k genes...\n')

    sys.stderr.write('\n')
    sys.stderr.write('Processed ' + str(j) + ' genes, ' + str(exons_num) + ' CDSs\n')
    sys.stderr.write('Found ' + str(len(zero_size_cds_lines)) + ' zero-size CDS, ' + str(len(zero_size_gene_lines)) + ' zero-size genes.\n')

    return not_approved_gene_names


def main():
    if len(sys.argv) < 1:
        sys.stderr.write('The script writes all CDS for all known Ensembl genes, with associated gene symbols.\n')
        sys.stderr.write('When the gene name is found in HGNC, it get replaced with an approved name.\n')
        sys.stderr.write('If the gene is not charactirized (like LOC729737), this symbol is just kept as is.\n')
        sys.stderr.write('\n')
        sys.stderr.write('Usage:\n')
        sys.stderr.write('    ' + __file__ + ' HGNC_gene_synonyms.txt [file_to_write_not_approved_genes.txt] < UCSC_knownGene.txt > UCSC_HGNC_exons.bed\n')
        sys.stderr.write('\n')
        sys.stderr.write('    where HGNC_gene_synonyms.txt (from http://www.genenames.org/cgi-bin/download) is:\n')
        sys.stderr.write('      #Approved Symbol  Previous Symbols                    Synonyms                          Chromosome   Ensembl Gene ID   UCSC ID(supplied by UCSC)\n')
        sys.stderr.write('      OR7E26P           OR7E67P, OR7E69P, OR7E70P, OR7E68P  OR1-51, OR1-72, OR1-73, OR912-95  19q13.43	 ENSG00000121410   uc002qsg.3\n')
        sys.stderr.write('      ...\n')
        sys.stderr.write('\n')
        sys.stderr.write('    and UCSC_knownGene.txt (from http://genome.ucsc.edu/cgi-bin/hgTables) is:\n')
        sys.stderr.write('      #hg19.knownGene.name  hg19.knownGene.chrom  hg19.knownGene.strand  hg19.knownGene.txStart  hg19.knownGene.txEnd  hg19.knownGene.exonCount  hg19.knownGene.exonStarts  hg19.knownGene.exonEnds  hg19.kgXref.geneSymbol\n')
        sys.stderr.write('      uc001aaa.3	          chr1	                +	                   11873                   14409                 3                         11873,12612,13220,	      12227,12721,14409,	   DDX11L1\n')
        sys.stderr.write('      ...\n')
        sys.stderr.write('    or Ensembl.gtf (ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz)')
        sys.stderr.write('      1       pseudogene      gene    11869   14412   .       +       .       gene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene";')
        sys.stderr.write('      1       processed_transcript    transcript      11869   14409   .       +       .       gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene"; transcript_name "DDX11L1-002"; transcript_source "havana";')
        sys.stderr.write('\n')
        sys.stderr.write('See more info in http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Making+the+full+list+of+UCSC+exons+with+approved+HUGO+gene+symbols\n')
        sys.exit(1)

    if is_local:
        sys.stderr.write('Local: will run only for chr21\n')
        sys.stderr.write('\n')

    synonyms_fpath = sys.argv[1]
    approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym = _read_approved_genes(synonyms_fpath)

    not_approved_fpath = sys.argv[2] if len(sys.argv) > 2 else None

    with sys.stdin as inp, sys.stdout as out:
        l = inp.readline()
        if l.startswith('#!genome-build'):
            not_approved_gene_names = _proc_ensembl(inp, out, approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym)
        else:
            not_approved_gene_names = _proc_ucsc(inp, out, approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym)

    sys.stderr.write('\n')
    sys.stderr.write('Not approved ' + str(len(not_approved_gene_names.keys())) + ' genes.\n')
    if not_approved_fpath:
        with open(not_approved_fpath, 'w') as f:
            f.write('#Searched as\tStatus\n')
            f.writelines((gn + '\t' + st + '\n' for gn, st in not_approved_gene_names.items()))
        sys.stderr.write('Saved not approved to ' + not_approved_fpath + '\n')


if __name__ == '__main__':
    main()