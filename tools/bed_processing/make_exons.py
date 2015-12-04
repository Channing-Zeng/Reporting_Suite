#!/usr/bin/env python
import bcbio_postproc  # do not remove it: checking for python version and adding site dirs inside

from collections import defaultdict, OrderedDict
import sys
from traceback import format_exc
from source.file_utils import adjust_path, verify_file
from source.logger import err

us_syn_path = '/ngs/reference_data/genomes/Hsapiens/common/HGNC_gene_synonyms.txt'

ALL_EXONS = True


def main():
    if len(sys.argv) < 2:
        err('The script writes all CDS, stop codon, and ncRNA exon regions for all known Ensembl genes, with '
            'associated gene symbols.')
        err('When the gene name is found in HGNC, it get replaced with an approved name.')
        err('If the gene is not charactirized (like LOC729737), this symbol is just kept as is.')
        err('')
        err('Usage:')
        err('    ' + __file__ + ' Ensembl.gtf [HGNC_gene_synonyms.txt=' + us_syn_path + '] [additional_feature_list]'
                                                                                        ' > Exons.bed')
        err('')
        err('   where HGNC_gene_synonyms.txt (from http://www.genenames.org/cgi-bin/download) is:')
        err('     #Approved Symbol  Previous Symbols                    Synonyms                          '
            'Chromosome   Ensembl Gene ID   UCSC ID(supplied by UCSC)')
        err('     OR7E26P           OR7E67P, OR7E69P, OR7E70P, OR7E68P  OR1-51, OR1-72, OR1-73, OR912-95  '
            '19q13.43	    ENSG00000121410   uc002qsg.3')
        err('     ...')
        err('')
        err('   feature_list is by default empty, but could be transcript')
        err('')
        err('   and UCSC_knownGene.txt (from http://genome.ucsc.edu/cgi-bin/hgTables) is:')
        err('     #hg19.knownGene.name  hg19.knownGene.chrom  hg19.knownGene.strand  hg19.knownGene.txStart  '
            'hg19.knownGene.txEnd  hg19.knownGene.exonCount  hg19.knownGene.exonStarts  hg19.knownGene.exonEnds'
            '  hg19.kgXref.geneSymbol')
        err('     uc001aaa.3	          chr1	                +	                   11873                   '
            '14409                 3                         11873,12612,13220,	      12227,12721,14409,	   DDX11L1')
        err('     ...')
        err('   or Ensembl.gtf (ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz)')
        err('     1  pseudogene            gene        11869  14412  .  +  .  gene_id "ENSG00000223972"; '
            'gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene";')
        err('     1  processed_transcript  transcript  11869  14409  .  +  .  gene_id "ENSG00000223972"; '
            'transcript_id "ENST00000456328"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype '
            '"pseudogene"; transcript_name "DDX11L1-002"; transcript_source "havana";')
        err('     ...')
        err('')
        err('   Writes to Exons.bed')
        err('')
        err('See more info in http://wiki.rd.astrazeneca.net/display/NG/SOP+-+Making+the+full+list+of+UCSC+exons'
            '+with+approved+HUGO+gene+symbols')
        sys.exit(1)

    input_fpath = verify_file(sys.argv[1])

    synonyms_fpath = None
    if len(sys.argv) > 2:
        synonyms_fpath = verify_file(sys.argv[2])
        err('Synonyms file provided ' + synonyms_fpath + '')
    else:
        err('No synonyms file provided, skipping approving')

    not_approved_fpath = None
    if len(sys.argv) > 3:
        not_approved_fpath = adjust_path(sys.argv[3])

    out = sys.stdout
    with open(input_fpath) as inp:
        l = inp.readline()
        if l.startswith('#!genome-build'):
            gene_by_name = _proc_ensembl(inp, out)
        else:
            gene_by_name = _proc_ucsc(inp, out)

    if synonyms_fpath and synonyms_fpath != "''":
        gene_by_name, not_approved_gene_names = _approve(gene_by_name, synonyms_fpath)

        err('')
        err('Not approved by HGNC - ' + str(len(not_approved_gene_names)) + ' genes.')
        if not_approved_fpath:
            with open(not_approved_fpath, 'w') as f:
                f.write('#Searched as\tStatus\n')
                f.writelines((l + '\n' for l in not_approved_gene_names))
            err('Saved not approved to ' + not_approved_fpath)

        with open('serialized_genes.txt', 'w') as f:
            for g in gene_by_name.values():
                f.write(str(g) + '\t' + str(g.db_id) + '\n')
                for e in g.exons:
                    f.write('\t' + str(e) + '\n')

    for g in gene_by_name.values():
        out.write(g)
        for e in g.exons:
            out.write(e)


def _approve(gene_by_name, synonyms_fpath):
    approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym = \
        read_approved_genes(synonyms_fpath)

    not_approved_gene_names = list()
    gene_after_approving_by_name = OrderedDict()
    total_approved = 0
    total_not_approved = 0
    j = 0
    for g in gene_by_name.values():
        if len(g.exons) == 0:
            continue

        gene_after_approving_by_name[g.name] = g
        if is_approved_symbol(g.name, approved_gene_by_name):
            gene_after_approving_by_name[g.name] = g
            total_approved += 1
        else:
            not_approved_gene_names.append(g.name)
            total_not_approved += 1

        j += 1
        if j % 1000 == 0:
            err('processed ' + str(j / 1000) + 'k genes...')

    err('-----')
    err('Total: ' + str(j))
    if approved_gene_by_name:
        err('Total approved: ' + str(total_approved))
        err('Total not approved: ' + str(total_not_approved))
    err()
    err('Saving genes...')

    gene_features = 0
    features_counter = defaultdict(int)
    biotypes_counter = defaultdict(int)
    no_exon_gene_num = 0

    filtered_gene_after_approving_by_name = OrderedDict()
    for g in gene_after_approving_by_name.values():
        if len(g.exons) == 0:
            no_exon_gene_num += 1
        else:
            filtered_gene_after_approving_by_name[g.name] = g

            gene_features += 1
            features_counter[g.feature] += 1
            biotypes_counter[g.biotype] += 1

            for e in g.exons:
                features_counter[e.feature] += 1

                if e.feature == 'exon': e.feature = 'Exon'
                elif e.feature == 'stop_codon': e.feature = 'CDS'
                else: e.feature = e.feature[0].upper() + e.feature[1:]

    err('Skipped {} genes with no sub-features.'.format(no_exon_gene_num))
    err('Approved {} genes, including:'.format(gene_features))
    err('    Gene: {}'.format(features_counter['Gene']))
    err('    Multi_Gene: {}'.format(features_counter['Multi_Gene']))
    err('')

    err('Out of total: {} protein coding genes, {} ncRNA genes, including:'.format(
        biotypes_counter['protein_coding'], sum(biotypes_counter.values()) - biotypes_counter['protein_coding']))
    for bt, cnt in biotypes_counter.items():
        if bt != 'protein_coding':
            err('    ' + bt + ': ' + str(cnt))

    err()
    if ALL_EXONS:
        err('Found {} exons.'.format(features_counter['exon']))
    else:
        err('Also found {} CDS, {} stop codons, and {} ncRNA exons.'.format(
            features_counter['CDS'], features_counter['stop_codon'], features_counter['exon']))

    return filtered_gene_after_approving_by_name, not_approved_gene_names


class ApprovedGene:
    def __init__(self, name, prev_names, synonyms, chrom, ucsc_id=None, ensembl_id=None):
        self.name = name
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

    err('  Notice: cannot parse chromosome ' + chrom)
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


def read_approved_genes(synonyms_fpath):
    approved_gene_by_name = dict()
    approved_gnames_by_prev_gname = defaultdict(list)
    approved_gnames_by_synonym = defaultdict(list)

    err('Parsing HGNC database ' + synonyms_fpath + '...')
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
        err('  Processed ' + str(i) + ' lines from ' + synonyms_fpath)
        err()

    return approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym


def _check_gene_symbol(approved_gene, gene_symbol, db_id, chrom):
    if db_id and db_id != approved_gene.db_id:
        # sys.stderr.write('Discordant db ids for ' + gene_symbol + ': db id = ' +
        # str(db_id) + ', in HGNC it is ' + str(approved_gene.db_id) + '\n')
        pass
    # else:
        # sys.stderr.write('Accordant db ids for ' + gene_symbol + ': db id = ' +
        # str(db_id) + ', in HGNC it is ' + str(approved_gene.db_id) + '\n')

    if chrom and chrom != approved_gene.chrom:
        # sys.stderr.write('Discordant chroms for ' + gene_symbol + ': chrom = ' +
        # chrom + ', in HGNC chrom is ' + approved_gene.chrom + '\n')
        return None

    return approved_gene


def get_approved_gene_symbol(approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym,
                             gene_symbol, db_id='', db_chrom='', indent=''):
    if gene_symbol in approved_gene_by_name:
        if _check_gene_symbol(approved_gene_by_name[gene_symbol], gene_symbol, db_id, db_chrom):
            return approved_gene_by_name[gene_symbol].name, None

    err(indent + 'Gene name ' + gene_symbol + ' is not approved, searching for an approved version... ',
        ending='', print_date=False)

    def _get_approved_genes_by_kind(approved_genes, kind):
        if not approved_genes:
            return 'NOT FOUND'

        if len(approved_genes) > 1:
            approved_genes_same_ucsc = [g for g in approved_genes if g.db_id == db_id]

            if len(approved_genes_same_ucsc) > 1:
                err(' ERROR: multiple approved gene names for ' + gene_symbol + ' (as ' + kind + ') with ucsc_id ' +
                    db_id + ': ' + ', '.join(g.name for g in approved_genes_same_ucsc) + '', print_date=False)
                return 'AMBIGUOUS'

            if len(approved_genes_same_ucsc) == 1:
                if _check_gene_symbol(approved_genes_same_ucsc[0], gene_symbol, db_id, db_chrom):
                    err(' found approved gene for ' + gene_symbol + ' (as ' + kind + ') with ucsc_id ' + db_id,
                        print_date=False)
                    return approved_genes_same_ucsc[0].name

            # Ok, no genes with same ucsc id, or not the same chromosome for them.

            approved_genes_same_chrom = [g for g in approved_genes if g.chrom == db_chrom]

            if len(approved_genes_same_chrom) > 1:
                err(' ERROR: multiple approved gene names for ' + gene_symbol + ' (as ' + kind + ') with chrom ' +
                    db_chrom + ', '.join(g.name for g in approved_genes_same_ucsc) + '', print_date=False)
                return 'AMBIGUOUS'

            if len(approved_genes_same_chrom) == 1:
                g = approved_genes_same_chrom[0]
                err(' only ' + g.name + ' for ' + gene_symbol + ' (as ' + kind + ') has the same chrom '
                    + db_chrom + ', picking it', print_date=False)
                if _check_gene_symbol(g, gene_symbol, db_id, db_chrom):
                    return g.name
                else:
                    return 'NOT FOUND'

            if len(approved_genes_same_chrom) == 0:
                err(' ERROR: no approved gene names for ' + gene_symbol + ' (as ' + kind + ') with same chrom '
                    + db_chrom + '', print_date=False)
                return 'NOT FOUND'

        if len(approved_genes) == 1:
            if _check_gene_symbol(approved_genes[0], gene_symbol, db_id, db_chrom):
                err(' found approved gene symbol for ' + gene_symbol + ': ' + approved_genes[0].name + ' (as '
                    + kind + ')', print_date=False)
                return approved_genes[0].name

        return 'NOT FOUND'

    res = _get_approved_genes_by_kind(approved_gnames_by_prev_gname.get(gene_symbol), 'prev')
    if res == 'AMBIGUOUS':
        return None, 'AMBIGUOUS\tAS PREV'
    elif res == 'NOT FOUND':
        res = _get_approved_genes_by_kind(approved_gnames_by_synonym.get(gene_symbol), 'synonym')
        if res == 'AMBIGUOUS':
            return None, res + '\tAS SYNONYM'
        if res == 'NOT FOUND':
            err(' not found.', print_date=False)
            return None, res
        else:
            err(indent + 'Finally found approved gene for ' + gene_symbol + ' (as synonym): ' + res, print_date=False)
            return res, None
    else:
        err(indent + 'Finally found approved gene for ' + gene_symbol + ' (as prev): ' + res, print_date=False)
        return res, None


def _proc_ucsc(inp, out):  #, approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym):
    gene_by_name = dict()

    for l in inp:
        if l and not l.startswith('#'):
            ucsc_id, ucsc_chrom, strand, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, gene_symbol =\
                l[:-1].split('\t')
            cdsStart = int(cdsStart)
            cdsEnd = int(cdsEnd)
            exonCount = int(exonCount)
            exonStarts = [int(v) + 1 for v in exonStarts.split(',') if v]
            exonEnds = map(int, filter(None, exonEnds.split(',')))

            # approved_gene_symbol, status = get_approved_gene_symbol(
            #     approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym,
            #     gene_symbol, ucsc_id, ucsc_chrom)
            #
            # if not approved_gene_symbol:
            #     not_approved_gene_names.append(gene_symbol + '\t' + status)
            #     if DO_APPROVE:
            #         continue
            #     else:
            #         approved_gene_symbol = gene_symbol

            txStart = exonStarts[0] - 1
            txEnd = exonEnds[exonCount - 1]

            # out.write('\t'.join([ucsc_chrom, str(min(txStart, cdsStart)), str(max(txEnd, cdsEnd)),
            #                      gene_symbol, '.', strand, 'Gene', '.']) + '\n')
            if gene_symbol not in gene_by_name:
                gene = Gene(gene_symbol, ucsc_chrom, min(txStart, cdsStart), str(max(txEnd, cdsEnd)), strand)
                gene_by_name[gene_symbol] = gene
            gene = gene_by_name[gene_symbol]

            for j, eStart, eEnd in zip(
                   range(exonCount),
                   [s for s in exonStarts if s],
                   [e for e in exonEnds if e]):
                eStart -= 1

                exon = None
                if eEnd <= cdsStart or eStart > cdsEnd:  # usually it means cdsStart = 0,
                                                         # no CDS for this gene, thus reporting exons
                    exon = Exon(gene, eStart, eEnd, '', 'Exon')
                    # out.write('\t'.join([ucsc_chrom, str(eStart), str(eEnd), gene_symbol, '.',
                    #                      strand, 'Exon', '.']) + '\n')
                else:
                    if cdsStart <= eStart:
                        exon = Exon(gene, eStart, eEnd, '', 'CDS')
                        # out.write('\t'.join([ucsc_chrom, str(eStart), str(eEnd), gene_symbol, '.',
                        #                      strand, 'CDS', '.']) + '\n')
                    elif eEnd > cdsStart:
                        exon = Exon(gene, cdsStart, eEnd, '', 'CDS')
                        # out.write('\t'.join([ucsc_chrom, str(cdsStart), str(eEnd), gene_symbol, '.',
                        #                      strand, 'CDS', '.']) + '\n')
                    else:
                        err('Warn: exon ' + str(eStart) + ':' + str(eEnd) +
                            ' does not contain CDS, CDS start = ' + str(cdsStart))
                if exon:
                    gene.exons.append(exon)

    return gene_by_name


class Gene:
    def __init__(self, name, chrom, start, end, strand, biotype='', db_id='', source=''):
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.biotype = biotype
        self.db_id = db_id
        self.feature = 'Gene'
        self.source = source

        self.approved_gname = None

        self.exons = []

    def __str__(self):
        fs = [self.chrom,
              '{}'.format(self.start) if self.start else '.',
              '{}'.format(self.end) if self.end else '.',
              self.name or '.', '.', self.strand or '.',
              self.feature or '.', self.biotype or '.']
        return '\t'.join(fs) + '\n'

    def __repr__(self):
        return '{self.name} {self.chrom}:{self.start}-{self.end} {self.biotype} ' \
               '{self.db_id} {self.source}'.format(self=self)


class Exon:
    def __init__(self, gene, start, end, biotype=None, feature=None):
        self.gene = gene
        self.start = start
        self.end = end
        self.biotype = biotype
        self.feature = feature

    def __str__(self):
        fs = [self.gene.chrom,
              '{}'.format(self.start) if self.start else '.',
              '{}'.format(self.end) if self.end else '.',
              self.gene.name or '.', '.', self.gene.strand or '.',
              self.feature or '.', self.biotype or '.']
        return '\t'.join(fs) + '\n'


def _rm_quotes(l):
    return l[1:-1]


def is_approved_symbol(gname, approved_gene_by_name):
    if gname not in approved_gene_by_name:
        # gname2 = gname.split('.')[0]
        # if gname != gname2:
        #     if gname2 not in approved_gene_by_name:
        return False
    return True


def _proc_ensembl(inp, out, additional_feature_list=None):
    if additional_feature_list is None:
        additional_feature_list = []

    err('additional_feature_list = ' + str(additional_feature_list))

    gene_by_name = OrderedDict()
    gene_by_id = OrderedDict()

    err('Parsing Ensembl input...')
    total_lines = 0
    total_non_coding_genes = 0

    for l in inp:
        if l and not l.startswith('#'):
            chrom, _, feature, start, end, _, strand, _, props_line = l[:-1].split('\t')

            # if is_local():
            #     if chrom != '21':
            #         continue

            total_lines += 1
            if total_lines % 1000 == 0:
                sys.stderr.write(str(total_lines / 1000) + 'k lines, ' + str(len(gene_by_name)) +
                                 ' genes found\n')
                sys.stderr.flush()

            try:
                _prop_dict = dict((t.strip().split(' ')[0], ' '.join(t.strip().split(' ')[1:]))
                                  for t in props_line.split(';') if t.strip())
            except ValueError:
                sys.stderr.write(format_exc())
                sys.stderr.write(l)

            gene_symbol = _rm_quotes(_prop_dict['gene_name'])
            gene_id = _rm_quotes(_prop_dict['gene_id'])
            gene_biotype = _rm_quotes(_prop_dict['gene_biotype'])
            gene_source = _rm_quotes(_prop_dict['gene_source'])

            # if gene_symbol == 'PTENP1':
            #     sys.stderr.write('PTENP1\n')

            if not ALL_EXONS and gene_biotype not in [
                'protein_coding',
                'nonsense_mediated_decay',
                'non_stop_decay',
                'processed_transcript',
                'polymorphic_pseudogene',
                'sense_intronic',
                'sense_overlapping',
                'antisense',

            ] and not any(b in gene_biotype for b in ['RNA', 'IG_', 'TR_']):
                total_non_coding_genes += 1
                continue

            full_feature_list = ['gene', 'CDS', 'stop_codon', 'exon'] + additional_feature_list
            if ALL_EXONS:
                full_feature_list = ['gene', 'exon']
            # sys.stderr.write('Full feature list: ' + str(full_feature_list) + '\n')
            if feature not in full_feature_list:
                continue

            start, end = int(start) - 1, int(end)

            if int(end) <= int(start):
                sys.stderr.write('Error: start > end: ' + l + '\n')
                continue

            chrom = parse_ensembl_chrom(chrom)
            if not chrom:
                continue

            if feature == 'gene':
                # assert gene_biotype == biotype, 'Gene: gene_biotype "' + gene_biotype + '"
                # do not match biotype "' + biotype + '" for ' + gene_symbol

                gene = Gene(gene_symbol, chrom, start, end, strand,
                            gene_biotype, gene_id, gene_source)

                if gene.name in gene_by_name:
                    prev_gene = gene_by_name[gene.name]

                    if gene.source != prev_gene.source:
                        err('    Duplicated gene in different databases:')
                        err('        This: ' + gene.__repr__())
                        err('        Prev: ' + prev_gene.__repr__())
                        # answer = raw_input('Which one to pick? This (1), prev (2), longest (Enter): ')
                        #
                        # if answer == '1' or answer == '' and gene.end - gene.start >
                        # prev_gene.end - prev_gene.start:
                        #     del gene_by_name[prev_gene.name]
                        #     del gene_by_id[prev_gene.db_id]
                        #
                        # else:
                        #     continue

                        if gene.source == 'ensembl' or prev_gene.source == 'havana':
                            del gene_by_name[prev_gene.name]
                            del gene_by_id[prev_gene.db_id]
                            err('        Picking up this one.')

                        if prev_gene.source == 'ensembl' or gene.source == 'havana':
                            err('        Picking up previous one.')
                            continue

                    else:
                        err('    Duplicated gene in ' + gene.source + ':')
                        err('        ' + gene.__repr__())
                        prev_gene.start = min(prev_gene.start, gene.start)
                        prev_gene.end = max(prev_gene.end, gene.end)
                        prev_gene.feature = 'Multi_Gene'
                        continue

                    sys.stderr.write('\n')

                gene_by_name[gene_symbol] = gene
                gene_by_id[gene_id] = gene

            elif feature in ['CDS', 'stop_codon'] \
                    or feature == 'exon' and ('RNA' in gene_biotype or ALL_EXONS) \
                    or feature in additional_feature_list:
                assert gene_symbol in gene_by_name, 'Error: ' + feature + ' record before gene record ' + \
                        gene_symbol + ', ' + gene_id + '; gene_by_name: ' + str(gene_by_name.keys())
                gene = gene_by_name[gene_symbol]
                if gene.gene_id == gene_id:
                    assert gene_biotype == gene.biotype, feature + ': gene_biotype "' + gene_biotype + \
                         '" do not match biotype "' + gene.biotype + '" for ' + gene_symbol
                    exon = Exon(gene, start, end, gene_biotype, feature)
                    gene.exons.append(exon)

    err()
    err(
        'Processed ' +
        str(total_lines) + ' lines, ' +
        str(total_non_coding_genes) + ' non-coding genes skipped, ' +
        str(len(gene_by_name)) + ' coding genes found')
    err()
    return gene_by_name


if __name__ == '__main__':
    main()
