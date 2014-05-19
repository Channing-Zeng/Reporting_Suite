#!/usr/bin/env python
from itertools import count, izip

import sys
from Bio import SeqIO
from Bio.Seq import Seq


BRCA1_exons_starts = [41196311,41199659,41201137,41203079,41209068,41215349,41215890,41219624,41222944,41226347,41228504,41234420,41242960,41243451,41247862,41249260,41251791,41256138,41256884,41258472,41267742,41276033,41277287]
BRCA1_exons_ends = [41197819,41199720,41201211,41203134,41209152,41215390,41215968,41219712,41223255,41226538,41228631,41234592,41243049,41246877,41247939,41249306,41251897,41256278,41256973,41258550,41267796,41276132,41277500]
BRCA1_exons = [r for r in zip(BRCA1_exons_starts, BRCA1_exons_ends)]

BRCA2_exons_starts = [32889616,32890558,32893213,32899212,32900237,32900378,32900635,32903579,32905055,32906408,32910401,32918694,32920963,32928997,32930564,32931878,32936659,32937315,32944538,32945092,32950806,32953453,32953886,32954143,32968825,32971034,32972298]
BRCA2_exons_ends = [32889804,32890664,32893462,32899321,32900287,32900419,32900750,32903629,32905167,32907524,32915333,32918790,32921033,32929425,32930746,32932066,32936830,32937670,32944694,32945237,32950928,32953652,32954050,32954282,32969070,32971181,32973809]
BRCA2_exons = [r for r in zip(BRCA2_exons_starts, BRCA2_exons_ends)]

class Record():
    def __init__(self, line, gene, ref_aa, aa_pos, alt_aa, ref_nc, nc_pos, alt_nc):
        self.line = line
        self.gene = gene

        self.ref_aa = ref_aa
        self.aa_pos = aa_pos
        self.alt_aa = alt_aa

        self.ref_nc = ref_nc
        self.nc_pos = nc_pos
        self.alt_nc = alt_nc

        self.orig_ref_nc = self.ref_nc
        self.orig_alt_nc = self.alt_nc

        if gene == 'BRCA1':
            self.chr = 'chr17'
            self.ref_nc = str(Seq(str(self.ref_nc)).reverse_complement())
            self.alt_nc = str(Seq(str(self.alt_nc)).reverse_complement())
        if gene == 'BRCA2':
            self.chr = 'chr13'

        self.id = '.'
        self.qual = '.'
        self.filter = '.'
        self.info = ''

    @staticmethod
    def create(l):
        gene, aa_change, nc_change, loc = l.split()

        # self.chr, self.chr_pos = loc.split(':')
        # self.chr_pos = int(self.chr_pos)

        print nc_change
        print nc_change[0:-3]
        ref_aa, aa_pos, alt_aa = aa_change[0], int(aa_change[1:-1]), aa_change[-1]
        ref_nc, nc_pos, alt_nc = nc_change[-3], int(nc_change[0:-3]), nc_change[-1]

        return Record(l, gene, ref_aa, aa_pos, alt_aa, ref_nc, nc_pos, alt_nc)



with open('BRCA_FM.txt') as fm_f,\
     open('chromFa/chr17.fa') as chr17_f,\
     open('chromFa/chr13.fa') as chr13_f,\
     open('BRCA_FM.vcf', 'w') as fm_vcf:
# with open('BIC_BRCA_VUS.txt') as fm_f,\
#      open('chromFa/chr17.fa') as chr17_f,\
#      open('chromFa/chr13.fa') as chr13_f,\
#      open('BIC_BRCA_VUS.vcf', 'w') as fm_vcf:
    chr17_seq = SeqIO.read(chr17_f, 'fasta').seq
    chr13_seq = SeqIO.read(chr13_f, 'fasta').seq

    fm_recs = []
    for line in fm_f:
        if not line.startswith('#'):
            fm_recs.append(Record.create(line))

    # r = Record('BRCA2	S976I	c.2926TC>AT', 'BRCA2', 'S', 976, 'I', 'TC', 2926, 'AT')
    # fm_recs.append(r)
    # r = Record('BRCA2	A3205L	c.9613GC>CT', 'BRCA2', 'A', 3205, 'L', 'GC', 9613, 'CT')
    # fm_recs.append(r)
# BRCA2	S976I	c.2926TC>AT?


    fm_recs = sorted(fm_recs, key=lambda r: (r.gene, r.nc_pos))
    fm_brca1_recs = [rec for rec in fm_recs if rec.gene == 'BRCA1']
    fm_brca2_recs = [rec for rec in fm_recs if rec.gene == 'BRCA2']

    fm_vcf.write('##fileformat=VCFv4.1\n')
    fm_vcf.write('##INFO=<ID=FM_GENE,Number=1,Type=String,Description="Gene name">\n')
    fm_vcf.write('##INFO=<ID=FM_NUC_CHANGE,Number=1,Type=String,Description="R>A">\n')
    fm_vcf.write('##INFO=<ID=FM_AA_CHANGE,Number=1,Type=String,Description="R>A">\n')
    fm_vcf.write('##INFO=<ID=FM_TRANSCRIPT_POS,Number=1,Type=Integer,Description="">\n')
    fm_vcf.write('##INFO=<ID=FM_PROTEIN_POS,Number=1,Type=Integer,Description="">\n')
    fm_vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

    BRCA1_gb_data = SeqIO.read('BRCA1.gb', 'genbank')
    BRCA2_gb_data = SeqIO.read('BRCA2.gb', 'genbank')
    BRCA1_cdna_seq = BRCA1_gb_data.seq
    BRCA2_cdna_seq = BRCA2_gb_data.seq
    BRCA1_cds = next(f for f in BRCA1_gb_data.features if f.type == 'CDS')
    BRCA2_cds = next(f for f in BRCA2_gb_data.features if f.type == 'CDS')

    # Print chromosome fasta
    BRCA1_chr_cds_seq_raw = ''
    BRCA1_rc_exons = []
    BRCA1_raw_positions = []
    for i, (st, en) in enumerate(BRCA1_exons):
        BRCA1_chr_cds_seq_raw += str(chr17_seq[st:en])
        for j in range(st, en):
            BRCA1_raw_positions.append(j)
    BRCA1_positions = BRCA1_raw_positions[::-1][BRCA1_cds.location.start:BRCA1_cds.location.end]
    print BRCA1_raw_positions
    print BRCA1_positions
    BRCA1_chr_cds_seq = Seq(BRCA1_chr_cds_seq_raw).reverse_complement()[BRCA1_cds.location.start:BRCA1_cds.location.end]
    print BRCA1_chr_cds_seq  # chr17_seq[BRCA1_chr_start + cds.location.start:BRCA1_chr_start + cds.location.end].complement()
    print ''.join((str(Seq(chr17_seq[pos]).complement())) for pos in BRCA1_positions)


    BRCA2_chr_cds_seq = ''
    BRCA2_rc_exons = []
    BRCA2_positions = []
    for i, (st, en) in enumerate(BRCA2_exons):
        BRCA2_chr_cds_seq += str(chr13_seq[st:en])
        for j in range(st, en):
            BRCA2_positions.append(j)
    BRCA2_chr_cds_seq = BRCA2_chr_cds_seq[BRCA2_cds.location.start:BRCA2_cds.location.end]
    BRCA2_positions = BRCA2_positions[BRCA2_cds.location.start:BRCA2_cds.location.end]
    print BRCA2_positions
    print ''.join(BRCA2_chr_cds_seq)
    print ''.join((str(chr17_seq[pos]) for pos in BRCA2_positions))


    for cds, cdna_seq, recs, positions, chrs in [
            (BRCA1_cds, BRCA1_cdna_seq, fm_brca1_recs, BRCA1_positions, chr17_seq),
            (BRCA2_cds, BRCA2_cdna_seq, fm_brca2_recs, BRCA2_positions, chr13_seq)
    ]:
        # continue
        # Print chromosome vars
        # vars = ['-' for _ in range(cds.location.start + len(cdna_seq))]
        # for rec in fm_brca1_recs:
        #     vars[cds.location.start + rec.nc_pos - 1] =\
        #         str(Seq(chrs[rec.chr_pos - 1]).reverse_complement())
        # print ''.join(vars[cds.location.start:cds.location.end])


        # Print CDA
        for i, cdna_nuc in izip(count(), cdna_seq[cds.location.start:cds.location.end]):
            sys.stdout.write(cdna_nuc)
        print ''


        # Print variants
        vars = ['-' for _ in cdna_seq]
        for rec in fm_brca1_recs:
            vars[cds.location.start + rec.nc_pos - 1] = rec.ref_nc
        print ''.join(vars[cds.location.start:cds.location.end])


        for rec in recs:
            rec.chr_pos = positions[rec.nc_pos]

        for rec in sorted(recs, key=lambda r: int(r.chr_pos)):
            fm_vcf.write(rec.chr + '\t')
            fm_vcf.write(str(rec.chr_pos) + '\t')
            fm_vcf.write('.\t')
            fm_vcf.write(rec.ref_nc + '\t')
            fm_vcf.write(rec.alt_nc + '\t')
            fm_vcf.write('.\t')
            fm_vcf.write('.\t')
            fm_vcf.write('FM_GENE=%s;' % rec.gene)
            fm_vcf.write('FM_NUC_CHANGE=%s;' % (rec.orig_ref_nc + '>' + rec.orig_alt_nc))
            fm_vcf.write('FM_AA_CHANGE=%s;' % (rec.ref_aa + '>' + rec.alt_aa))
            fm_vcf.write('FM_TRANSCRIPT_POS=%d;' % rec.nc_pos)
            fm_vcf.write('FM_PROTEIN_POS=%d\n' % rec.aa_pos)


# with open('BRCA_FM.txt') as fm_f,\
#      open('chromFa/chr17.fa') as chr17_f,\
#      open('BRCA_FM.vcf', 'a') as fm_vcf:
#     chr17_seq = SeqIO.read(chr17_f, 'fasta').seq
#     print chr17_seq[41276112 + 1]
#
#     fm_recs = []
#     for line in fm_f:
#         fm_recs.append(Record(line))
#     fm_recs = sorted(fm_recs, key=lambda r: (r.gene, r.nc_pos))
#     fm_brca1_recs = [rec for rec in fm_recs if rec.gene == 'BRCA1']
#     fm_brca2_recs = [rec for rec in fm_recs if rec.gene == 'BRCA2']
#
#     fm_vcf.write('##fileformat=VCFv4.1\n')
#     fm_vcf.write('##INFO=<ID=FM_GENE,Number=1,Type=String,Description="Gene name">\n')
#     fm_vcf.write('##INFO=<ID=FM_NUC_CHANGE,Number=1,Type=String,Description="R>A">\n')
#     fm_vcf.write('##INFO=<ID=FM_AA_CHANGE,Number=1,Type=String,Description="R>A">\n')
#     fm_vcf.write('##INFO=<ID=FM_TRANSCRIPT_POS,Number=1,Type=Integer,Description="">\n')
#     fm_vcf.write('##INFO=<ID=FM_PROTEIN_POS,Number=1,Type=Integer,Description="">\n')
#     fm_vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
#
#     gb_data = SeqIO.read('BRCA1.gb', 'genbank')
#     cdna_seq = gb_data.seq
#     cds = next(f for f in gb_data.features if f.type == 'CDS')
#
#     # Print chromosome fasta
#     chr_cds_seq_raw = ''
#     rc_exons = []
#     raw_positions = []
#     for i, (st, en) in enumerate(exons):
#         chr_cds_seq_raw += str(chr17_seq[st:en])
#         for j in range(st, en):
#             raw_positions.append(j)
#
#     print raw_positions
#     positions = raw_positions[::-1][cds.location.start:cds.location.end]
#     print positions
#
#     chr_cds_seq = Seq(chr_cds_seq_raw).reverse_complement()[cds.location.start:cds.location.end]
#     print chr_cds_seq  # chr17_seq[BRCA1_chr_start + cds.location.start:BRCA1_chr_start + cds.location.end].complement()
#     print ''.join((str(Seq(chr17_seq[pos]).complement())) for pos in positions)
#
#     # Print chromosome vars
#     vars = ['-' for _ in range(cds.location.start + len(cdna_seq))]
#     for rec in fm_brca1_recs:
#         vars[cds.location.start + rec.nc_pos - 1] =\
#             str(Seq(chr17_seq[rec.chr_pos - 1]).reverse_complement())
#     print ''.join(vars[cds.location.start:cds.location.end])
#
#
#     # Print CDA
#     for i, cdna_nuc in izip(count(), cdna_seq[cds.location.start:cds.location.end]):
#         sys.stdout.write(cdna_nuc)
#     print ''
#
#
#     # Print variants
#     vars = ['-' for _ in cdna_seq]
#     for rec in fm_brca1_recs:
#         vars[cds.location.start + rec.nc_pos - 1] = rec.ref_nc
#     print ''.join(vars[cds.location.start:cds.location.end])
#
#
#     for rec in fm_brca1_recs:
#         print rec.chr_pos, positions[rec.nc_pos - 2]
#
#     for rec in fm_brca1_recs:
#         fm_vcf.write(rec.chr + '\t')
#         fm_vcf.write(str(rec.chr_pos) + '\t')
#         fm_vcf.write('.\t')
#         fm_vcf.write(str(Seq(str(rec.ref_nc)).reverse_complement()) + '\t')
#         fm_vcf.write(str(Seq(str(rec.alt_nc)).reverse_complement()) + '\t')
#         fm_vcf.write('.\t')
#         fm_vcf.write('.\t')
#         fm_vcf.write('FM_GENE=%s;' % rec.gene)
#         fm_vcf.write('FM_NUC_CHANGE=%s;' % (rec.ref_nc + '>' + rec.alt_nc))
#         fm_vcf.write('FM_AA_CHANGE=%s;' % (rec.ref_aa + '>' + rec.alt_aa))
#         fm_vcf.write('FM_TRANSCRIPT_POS=%d;' % rec.nc_pos)
#         fm_vcf.write('FM_PROTEIN_POS=%d\t' % rec.aa_pos)
#         fm_vcf.write('\n')
#















