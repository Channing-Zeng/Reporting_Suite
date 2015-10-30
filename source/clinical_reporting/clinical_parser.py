from collections import OrderedDict
from itertools import izip
from json import load
from os.path import join

import source
from source import verify_file, info
from source.clinical_reporting.solvebio_mutations import query_mutations
from source.logger import warn, err
from source.reporting.reporting import SampleReport
from source.targetcov.Region import SortableByChrom
from source.targetcov.flag_regions import get_depth_cutoff
from source.targetcov.summarize_targetcov import get_float_val, get_val


ACTIONABLE_GENES_FPATH = join(__file__, '..', 'db', 'broad_db.tsv')


class KeyGene:
    def __init__(self, name, ave_depth=None):
        self.name = name
        self.chrom = None
        self.start = None
        self.end = None
        self.cdss = []
        self.ave_depth = ave_depth
        self.cov_by_threshs = dict()
        self.mutations = []
        self.seq2c_event = None

    def __str__(self):
        return self.name

    def __hash__(self):
        return hash(self.name)


class Mutation(SortableByChrom):
    def __init__(self, chrom, genome):
        SortableByChrom.__init__(self, chrom, genome=genome)
        self.gene = None
        self.transcript = None
        self.codon_change = None
        self.chrom = chrom
        self.pos = None
        self.ref = None
        self.alt = None
        self.depth = None
        self.freq = None
        self.aa_change = None
        self.aa_len = None
        self.eff_type = None
        self.status = None
        self.reason = None
        self.cosmic_ids = []
        self.dbsnp_ids = []
        self.solvebio = None

    def __str__(self):
        return str(self.gene) + ' ' + str(self.chrom) + ':' + \
               str(self.pos) + ' ' + str(self.ref) + '>' + str(self.alt)

    def __hash__(self):
        return hash((self.chrom, self.pos, self.ref, self.alt))

    def get_chrom_key(self):
        return SortableByChrom.get_key(self)

    def get_key(self):
        return SortableByChrom.get_key(self), self.pos, self.ref, self.alt


class Seq2CEvent:
    def __init__(self, gene=None, amp_del=None, fragment=None, ab_log2r=None, log2r=None):
        self.gene = gene
        self.amp_del = amp_del  # Del, Amp
        self.fragment = fragment  # Whole, BP
        self.ab_log2r = ab_log2r
        self.log2r = log2r

    def is_amp(self):
        return self.amp_del == 'Amp'

    def is_del(self):
        return self.amp_del == 'Del'

    def __str__(self):
        return str(self.gene) + ' ' + str(self.log2r if self.log2r is not None else self.ab_log2r)\
               + ' ' + str(self.amp_del) + ' ' + str(self.fragment)

    def __hash__(self):
        return hash((self.gene, self.fragment, self.amp_del))


class CDS:
    def __init__(self):
        self.start = None
        self.end = None

        self.ave_depth = None
        self.cov_by_threshs = dict()


class Patient:
    def __init__(self, gender):
        self.gender = gender


class Target:
    def __init__(self, coverage_percent=None, type_=None, bed_fpath=None):
        self.coverage_percent = coverage_percent
        self.type = type_
        self.bed_fpath = bed_fpath


class ClinicalSampleInfo:
    def __init__(self, cnf):
        self.cnf = cnf
        self.sample = source.BaseSample(cnf.sample, cnf.output_dir,
            targqc_dirpath=cnf.targqc_dirpath, clinical_report_dirpath=cnf.output_dir,
            normal_match=cnf.match_sample_name)
        info('Sample: ' + str(cnf.sample))
        info('Match sample name: ' + str(cnf.match_sample_name))
        info()

        info('Preparing data for a clinical report for AZ 300 key genes ' + str(self.cnf.key_genes) + ', sample ' + self.sample.name)
        self.key_gene_by_name = dict()
        for gene_name in get_key_genes(self.cnf.key_genes):
            self.key_gene_by_name[gene_name] = KeyGene(gene_name)

        info('Parsing target and patient info...')
        self.depth_cutoff = None
        self.patient = Patient(gender=get_gender(self.sample, self.sample.targetcov_json_fpath))
        self.target = Target(
            type_=cnf.target_type,
            coverage_percent=get_target_fraction(self.sample, self.sample.targetcov_json_fpath),
            bed_fpath=self.cnf.bed_fpath)

        info('Parsing TargetCov key genes stats...')
        self.ave_depth = get_ave_coverage(self.sample, self.sample.targetcov_json_fpath)
        self.depth_cutoff = get_depth_cutoff(self.ave_depth, self.cnf.coverage_reports.depth_thresholds)
        self.parse_targetseq_detailed_report()

        info('Parsing actionable genes...')
        self.actionable_genes_dict = parse_broad_actionable()

        info('Parsing mutations...')
        self.total_variants = get_total_variants_number(self.sample, self.cnf.varqc_json_fpath)
        self.mutations = self.parse_mutations(self.cnf.mutations_fpath)
        for mut in self.mutations:
            self.key_gene_by_name[mut.gene.name].mutations.append(mut)

        info('Retrieving SolveBio...')
        self.get_mut_info_from_solvebio()

        info('Parsing Seq2C...')
        self.seq2c_events_by_gene_name = None
        if not self.cnf.seq2c_tsv_fpath:
            warn('No Seq2C results provided by option --seq2c, skipping plotting Seq2C')
        else:
            self.seq2c_events_by_gene_name = self.parse_seq2c_report(self.cnf.seq2c_tsv_fpath)


    def get_mut_info_from_solvebio(self):
        query_mutations(self.cnf, self.mutations)

    def parse_seq2c_report(self, seq2c_tsv_fpath):
        seq2c_events_by_gene_name = dict()
        if not verify_file(seq2c_tsv_fpath, silent=True):
            return None

        with open(seq2c_tsv_fpath) as f_inp:
            for i, l in enumerate(f_inp):
                if i == 0: continue
                fs = l[:-1].split('\t')
                sname, gname = fs[0], fs[1]
                if gname not in self.key_gene_by_name: continue
                if sname != self.sample.name: continue

                sname, gname, chrom, start, end, length, log2r, sig, fragment, amp_del, ab_seg, total_seg, \
                    ab_log2r, log2r_diff, ab_seg_loc, ab_samples, ab_samples_pcnt = fs[:17]

                event = Seq2CEvent(
                    gene=self.key_gene_by_name[gname],
                    fragment=fragment or None,
                    ab_log2r=float(ab_log2r) if ab_log2r else None,
                    log2r=float(log2r) if log2r else None,
                    amp_del=amp_del or None)
                seq2c_events_by_gene_name[gname] = event

        for gn, event in seq2c_events_by_gene_name.items():
            self.key_gene_by_name[gn].seq2c_event = event

        return seq2c_events_by_gene_name

    def parse_targetseq_detailed_report(self):
        info('Preparing coverage stats key gene tables')
        with open(self.sample.targetcov_detailed_tsv) as f_inp:
            for l in f_inp:
                if l.startswith('#'):
                    continue

                fs = l.split('\t')  # Chr	Start	End	Size	Gene	Strand	Feature	Biotype	Min depth	Ave depth	Std dev	W/n 20% of ave depth	1x	5x	10x	25x	50x	100x	500x	1000x	5000x	10000x	50000x
                chrom, start, end, size, symbol, strand, feature, biotype, min_depth, ave_depth, std_dev, wn20pcnt = fs[:12]
                pcnt_val_by_thresh = fs[12:]
                symbol = get_val(symbol)
                if symbol not in self.key_gene_by_name:
                    continue

                chrom = get_val(chrom)
                if not chrom:
                    err('For gene ' + str(symbol) + ', chrom value if empty in line ' + l + ' in ' + self.sample.targetcov_detailed_tsv)
                    continue

                start, end = int(start), int(end)
                ave_depth = get_float_val(ave_depth)
                cov_by_threshs = dict((t, get_float_val(f)) for t, f in izip(self.cnf.coverage_reports.depth_thresholds, pcnt_val_by_thresh))

                if feature in ['Whole-Gene', 'Gene-Exon']:
                    gene = self.key_gene_by_name.get(symbol)
                    gene.chrom = chrom
                    gene.start = start
                    gene.end = end
                    gene.ave_depth = ave_depth
                    gene.cov_by_threshs = cov_by_threshs

                elif feature in ['CDS', 'Exon']:
                    cds = CDS()
                    cds.start = start
                    cds.end = end
                    cds.ave_depth = ave_depth
                    cds.cov_by_threshs = cov_by_threshs
                    gene = self.key_gene_by_name.get(symbol)
                    gene.cdss.append(cds)

        # Cleaning up records that are not found in the target gene panel,
        # so we don't know about them and don't even want to report them
        for gene in self.key_gene_by_name.values():
            if not gene.chrom and gene.name in self.key_gene_by_name:
                del self.key_gene_by_name[gene.name]


    def parse_mutations(self, mutations_fpath):
        mutations = []

        info('Preparing mutations stats for key gene tables')
        info('Checking ' + mutations_fpath)
        if not verify_file(mutations_fpath):
            mut_pass_ending = source.mut_pass_suffix + '.' + source.mut_file_ext
            mut_basename = mutations_fpath.split('.' + mut_pass_ending)[0]
            if self.sample.normal_match:
                mutations_fpath = mut_basename + '.' + source.mut_paired_suffix + '.' + mut_pass_ending
            else:
                mutations_fpath = mut_basename + '.' + source.mut_single_suffix + '.' + mut_pass_ending
            info('Checking ' + mutations_fpath)
            if not verify_file(mutations_fpath):
                err('Cannot find PASSed mutations fpath')
                return []

        info('Reading mutations from ' + mutations_fpath)
        alts_met_before = set()
        with open(mutations_fpath) as f:
            for i, l in enumerate(f):
                if i == 0:
                    continue
                fs = l.strip().split('\t')
                reason = None
                if len(fs) >= 67:
                    sample_name, chrom, start, ids, ref, alt, type_, effect, func, codon_change, aa_change, cdna_change, \
                        aa_len, gname, transcr_biotype, coding, transcript, exon, cosmic_cds_change, cosmic_aa_change, \
                        cosmic_cnt, end, depth, af, bias, pmean, pstd, qual, qstd, sbf, gmaf, vd, clnsif, oddratio, hiaf, \
                        mq, sn, adjaf, nm, shift3, msi, dbsnpbuildid, vtype, status1, paired_pval, paired_oddratiom, \
                        m_depth, m_af, m_vd, m_bias, m_pmean, m_pstd, m_qual, m_qstd, m_hiaf, m_mq, m_sn, m_adjaf, m_nm, \
                        n_sample, n_var, pcnt_sample, ave_af, filt, var_type, var_class, status = fs[:67]  # 67 of them
                    if len(fs) == 68:
                        reason = fs[67]
                else:
                    sample_name, chrom, start, ids, ref, alt, type_, effect, func, codon_change, aa_change, cdna_change, \
                        aa_len, gname, transcr_biotype, coding, transcript, exon, cosmic_cds_change, cosmic_aa_change, \
                        cosmic_cnt, end, depth, af, bias, pmean, pstd, qual, qstd, sbf, gmaf, vd, clnsif, oddratio, hiaf, \
                        mq, sn, adjaf, nm, shift3, msi, dbsnpbuildid, \
                        n_sample, n_var, pcnt_sample, ave_af, filt, var_type, var_class, status = fs[:50]  # 50 of them
                    if len(fs) == 51:
                        reason = fs[50]

                if sample_name == self.sample.name and gname in self.key_gene_by_name:
                    if (chrom, start, ref, alt) in alts_met_before:
                        continue
                    alts_met_before.add((chrom, start, ref, alt))

                    mut = Mutation(chrom=chrom, genome=self.cnf.genome.name)
                    mut.gene = self.key_gene_by_name[gname]
                    mut.transcript = transcript
                    mut.codon_change = codon_change
                    mut.aa_change = aa_change
                    mut.aa_len = aa_len
                    mut.pos = int(start)
                    mut.ref = ref
                    mut.alt = alt
                    mut.depth = int(depth)
                    mut.freq = float(af)
                    mut.dbsnp_ids = [''.join(c for c in id_ if c.isdigit()) for id_ in ids.split(';') if id_.startswith('rs')]
                    mut.cosmic_ids = [''.join(c for c in id_ if c.isdigit()) for id_ in ids.split(';') if id_.startswith('COS')]
                    mut.eff_type = (type_[0] + type_[1:].lower().replace('_', ' ')) if type_ else type_
                    mut.var_type = var_type
                    mut.var_class = var_class
                    mut.status = status
                    if reason:
                        reason = reason.replace('_', ' ')
                        if reason == 'actionable somatic':
                            reason = 'actionable som.'
                        if reason == 'actionable germline':
                            reason = 'actionable germ.'
                        reason = reason.replace('change', 'chg.')
                    mut.reason = reason

                    mutations.append(mut)
        info('Found ' + str(len(mutations)) + ' mutations in key genes')
        return mutations


def parse_broad_actionable():
    act_fpath = verify_file(ACTIONABLE_GENES_FPATH, is_critical=False, description='Actionable genes')
    if act_fpath:
        with open(act_fpath) as f:
            actionable_gene_table = [l.split('\t') for l in f.readlines()]
            return dict((l[0], l[1:]) for l in actionable_gene_table)


def get_key_genes(key_genes_fpath):
    key_genes_fpath = verify_file(key_genes_fpath, is_critical=True, description='300 AZ key genes')
    with open(key_genes_fpath) as f:
        key_gene_names = set([l.strip() for l in f.readlines() if l.strip() != ''])

    return key_gene_names


def get_target_fraction(sample, targqc_json_fpath):
    with open(targqc_json_fpath) as f:
        data = load(f, object_pairs_hook=OrderedDict)
    sr = SampleReport.load(data, sample, None)
    r = sr.find_record(sr.records, 'Percentage of target covered by at least 1 read')
    if not r:
        r = sr.find_record(sr.records, 'Percentage of genome covered by at least 1 read')
    return r.value if r else None


def get_gender(sample, targqc_json_fpath):
    with open(targqc_json_fpath) as f:
        data = load(f, object_pairs_hook=OrderedDict)
    sr = SampleReport.load(data, sample, None)
    r = sr.find_record(sr.records, 'Gender')

    if r:
        if r.value:
            if r.value.startswith('M'):
                return 'Male'
            elif r.value.startswith('F'):
                return 'Female'
            return 'Undetermined'
    return None


def get_ave_coverage(sample, targqc_json_fpath):
    with open(targqc_json_fpath) as f:
        data = load(f, object_pairs_hook=OrderedDict)
    sr = SampleReport.load(data, sample, None)
    r = sr.find_record(sr.records, 'Average target coverage depth')
    if not r:
        r = sr.find_record(sr.records, 'Average genome coverage depth')
    return r.value if r else None


def get_min_coverage(sample, targqc_json_fpath):
    with open(targqc_json_fpath) as f:
        data = load(f, object_pairs_hook=OrderedDict)
    sr = SampleReport.load(data, sample, None)
    r = sr.find_record(sr.records, 'Minimum target coverage depth')
    if not r:
        r = sr.find_record(sr.records, 'Minimum genome coverage depth')
    return r.value if r else None


def get_total_variants_number(sample, varqc_json_fpath):
    with open(varqc_json_fpath) as f:
        data = load(f, object_pairs_hook=OrderedDict)
    sr = SampleReport.load(data, sample, None)
    r = sr.find_record(sr.records, 'Total with rejected')
    return r.value if r else None
