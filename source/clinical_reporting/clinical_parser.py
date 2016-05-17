from collections import OrderedDict
from itertools import izip
from json import load
from os.path import join

import re

import source
from source import verify_file, info
from source.clinical_reporting.known_sv import fusions as known_fusions
from source.file_utils import verify_file, add_suffix, symlink_plus, remove_quotes, adjust_path, verify_dir
from source.clinical_reporting.solvebio_mutations import query_mutations
from source.logger import warn, err, critical
from source.reporting.reporting import SampleReport
from source.targetcov.Region import SortableByChrom
from source.targetcov.bam_and_bed_utils import get_gene_keys
from source.targetcov.flag_regions import get_depth_cutoff
from source.targetcov.summarize_targetcov import get_float_val, get_val
from source.targetcov.Region import get_chrom_order

ACTIONABLE_GENES_FPATH = join(__file__, '..', 'db', 'broad_db.tsv')


class KeyGene:
    def __init__(self, name, chrom=None, ave_depth=None):
        self.name = name
        self.chrom = chrom
        self.key = self.name, self.chrom
        self.start = None
        self.end = None
        self.transcript_id = None
        self.strand = None
        self.cdss = []
        self.ave_depth = ave_depth
        self.cov_by_threshs = dict()
        self.mutations = []
        self.seq2c_events = []
        self.sv_events = set()

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash(self.name)


class Mutation(SortableByChrom):
    def __init__(self, chrom, chrom_ref_order):
        SortableByChrom.__init__(self, chrom, chrom_ref_order)
        self.gene = None
        self.exon = None
        self.transcript = None
        self.codon_change = None
        self.cdna_change = None
        self.chrom = chrom
        self.pos = None
        self.ref = None
        self.alt = None
        self.aa_change = None
        self.aa_len = None
        self.eff_type = None
        self.status = None
        self.reason = None
        self.cosmic_id = None
        self.dbsnp_id = None
        self.solvebio = None

        self.depth = None
        self.freq = None

    def __str__(self):
        return str(self.gene) + ' ' + str(self.chrom) + ':' + \
               str(self.pos) + ' ' + str(self.ref) + '>' + str(self.alt)

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash((self.chrom, self.pos, self.ref, self.alt))

    def get_chrom_key(self):
        return SortableByChrom.get_key(self)

    def get_key(self):
        return SortableByChrom.get_key(self), self.pos, self.ref, self.alt, self.transcript


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

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash((self.gene, self.fragment, self.amp_del))

    def get_key(self):
        return self.gene, self.amp_del, self.log2r

class SVEvent(SortableByChrom):
    class Annotation:
        def __init__(self):
            self.type = None
            self.effect = None
            self.genes = []
            self.transcript = None

            self.known = False
            self.event = None

        def get_key(self):
            if self.event:
                return self.event.get_key(), self.type, self.effect, self.transcript

        @staticmethod
        def parse_annotation(string):
            fs = string.split('|')
            a = SVEvent.Annotation()
            a.type = fs[0]
            a.effect = fs[1]
            genes_val = fs[2]
            if a.type == 'BND':
                if '/' in genes_val:
                    a.genes = sorted(genes_val.split('/'))
                elif genes_val.count('-') == 1:
                    a.genes = sorted(genes_val.split('-'))
                else:
                    return None
            elif genes_val:
                a.genes = [genes_val]
            a.transcript = fs[3]
            a.exon_info = fs[4]
            return a

    @staticmethod
    def parse_sv_event(chr_order, **kwargs):  # caller  sample  chrom  start  end  svtype  known  end_gene  lof  annotation  split_read_support  paired_end_support
        e = SVEvent(chrom=kwargs.get('chrom'), chrom_ref_order=chr_order.get(kwargs.get('chrom')))
        e.caller = kwargs.get('caller')
        e.start = int(kwargs.get('start'))
        e.sample = kwargs.get('sample')
        e.end = int(kwargs.get('end')) if kwargs.get('end') else None

        e.type = None
        e.id = None
        e.mate_id = None
        svt = kwargs.get('svtype')
        if svt:  # BND:MantaBND:12:0:1:0:0:0:1:MantaBND:12:0:1:0:0:0:0 or BND:71_2:71_1 or
            e.type = svt.split(':', 1)[0]
            if e.type == 'BND' and ':' in svt:
                if 'MantaBND' in svt:
                    m = re.match(r'(?P<id1>MantaBND[:0-9]+):(?P<id2>MantaBND[:0-9]+)', svt.split(':', 1)[1])
                else:
                    m = re.match(r'(?P<id1>.+):(?P<id2>.+)', svt.split(':', 1)[1])
                e.id = m.group('id1')
                e.mate_id = m.group('id2')

        e.known_gene_val = kwargs.get('known')
        e.known_gene = e.known_gene_val.split('-with-')[1] if '-with-' in e.known_gene_val else e.known_gene_val
        if kwargs.get('end_gene'):
            e.end_genes = kwargs.get('end_gene').split(',')
        else:
            e.end_genes = []
        e.lof = kwargs.get('lof')
        e.annotations = []
        if kwargs.get('annotation'):
            for s in kwargs.get('annotation').split(','):
                a = SVEvent.Annotation.parse_annotation(s)
                if a:
                    assert a.type == e.type, 'Annotation type and event type does not match: ' + str(e.type) + ', ' + str(a.type)
                    e.annotations.append(a)

        e.split_read_support = kwargs.get('split_read_support').split(',') if kwargs.get('split_read_support') else []
        e.paired_end_support = kwargs.get('paired_end_support').split(',') if kwargs.get('paired_end_support') else []

        return e

        # lof_genes = []
        # if kwargs.get('end_gene'):
        #     for a in kwargs.get('end_gene').split(','):
        #         lof_genes.append(a[1:-1].split('|')[0])
        # with open(adjust_path('~/t.tsv'), 'a') as f:
        #     f.write(str(self.type) + '\t' + kwargs.get('known') + '\t' + kwargs.get('end_gene') + '\t' + ', '.join(lof_genes) + '\n')

    def __init__(self, chrom, chrom_ref_order):
        SortableByChrom.__init__(self, chrom, chrom_ref_order)
        self.caller = None
        self.start = None
        self.sample = None
        self.end = None
        self.type = None
        self.id = None
        self.mate_id = None
        self.known_gene_val = None
        self.known_gene = None
        self.end_genes = []
        self.lof = None
        self.annotations = []
        self.key_annotations = set()
        self.split_read_support = None
        self.paired_end_support = None

    def is_fusion(self):
        return self.type == 'BND'

    def get_possible_fusion_pairs(self):
        if self.type == 'BND':
            return [a.genes for a in self.annotations]

    def is_deletion(self):
        return self.type == 'DEL'

    def is_insertion(self):
        return self.type == 'INS'

    def is_duplication(self):
        return self.type == 'DUP'

    def __str__(self):
        return str(self.chrom) + ':' + str(self.start) + ' ' + str(self.annotations)

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash((self.caller, self.chrom, self.start, self.type, self.id, self.mate_id))

    def get_key(self):
        return self.chrom, self.type  #, tuple(tuple(sorted(a.genes)) for a in self.annotations)

    def get_chrom_key(self):
        return SortableByChrom.get_key(self)


# class FusionEvent(SVEvent):
#     def __init__(self, **kwargs):
#         SVEvent.__init__(self, **kwargs)
#         self.is_know_sv = False
#         self.fusion_pair = False


class CDS:
    def __init__(self):
        self.start = None
        self.end = None

        self.ave_depth = None
        self.cov_by_threshs = dict()


class Patient:
    def __init__(self, gender=None):
        self.gender = gender


class Target:
    def __init__(self, coverage_percent=None, type_=None, bed_fpath=None, targqc_link=None):
        self.coverage_percent = coverage_percent
        self.type = type_
        self.bed_fpath = bed_fpath
        self.targqc_link = targqc_link


def clinical_sample_info_from_bcbio_structure(cnf, bs, sample, is_target2wgs_comparison=False):
    clinical_report_caller = \
        bs.variant_callers.get('vardict') or \
        bs.variant_callers.get('vardict-java')
    if not clinical_report_caller:
        critical('No vardict or vardict-java variant caller in ' + str(bs.variant_callers.keys()))
    vardict_txt_fname = source.mut_fname_template.format(caller_name=clinical_report_caller.name)
    vardict_txt_fpath = join(bs.var_dirpath, vardict_txt_fname)
    mutations_fpath = add_suffix(vardict_txt_fpath, source.mut_pass_suffix)

    return ClinicalExperimentInfo(
        cnf, sample=sample, key_genes=cnf.key_genes,
        target_type=bs.target_type, bed_fpath=bs.bed, mutations_fpath=mutations_fpath, sv_fpath=sample.find_sv_fpath(),
        sv_vcf_fpath=verify_file(cnf.sv_vcf_fpath, is_critical=False) if cnf.sv_vcf_fpath else None,
        varqc_json_fpath=sample.get_varqc_fpath_by_callername(clinical_report_caller.name, ext='.json'),
        seq2c_tsv_fpath=bs.seq2c_fpath, project_name=bs.project_name,
        project_report_path=bs.project_report_html_fpath, targqc_report_path=bs.targqc_summary_fpath,
        is_target2wgs_comparison=is_target2wgs_comparison)


def clinical_sample_info_from_cnf(cnf):
    sample = source.BaseSample(cnf.sample, cnf.output_dir,
        targqc_dirpath=verify_dir(cnf.targqc_dirpath, silent=True),
        clinical_report_dirpath=cnf.output_dir,
        normal_match=cnf.match_sample_name)

    return ClinicalExperimentInfo(
        cnf, sample=sample, key_genes=cnf.key_genes,
        target_type=cnf.target_type,
        bed_fpath=verify_file(cnf.bed_fpath, is_critical=False) if cnf.bed_fpath else None,
        mutations_fpath=verify_file(cnf.mutations_fpath, is_critical=False) if cnf.mutations_fpath else None,
        sv_fpath=verify_file(cnf.sv_fpath, is_critical=False) if cnf.sv_fpath else None,
        sv_vcf_fpath=verify_file(cnf.sv_vcf_fpath, is_critical=False) if cnf.sv_vcf_fpath else None,
        varqc_json_fpath=verify_file(cnf.varqc_json_fpath, is_critical=False) if cnf.varqc_json_fpath else None,
        seq2c_tsv_fpath=verify_file(cnf.seq2c_tsv_fpath, is_critical=False) if cnf.seq2c_tsv_fpath else None,
        project_name=cnf.project_name,
        project_report_path=cnf.project_report_path,
        targqc_report_path=verify_file(cnf.targqc_report_path, silent=False) if cnf.targqc_report_path else None)


class GeneDict(dict):  # supports genes without specified chromosome
    def __init__(self, *args, **kwargs):
        super(GeneDict, self).__init__(**kwargs)
        self.gene_by_name = dict()
        for (gn, ch) in super(GeneDict, self).keys():
            if ch is None:
                self.gene_by_name[gn] = super(GeneDict, self)[(gn, ch)]
                del super(GeneDict, self)[(gn, ch)]

    def __delitem__(self, (gname, chrom)):
        if (gname, chrom) in super(GeneDict, self).keys():
            super(GeneDict, self).__delitem__((gname, chrom))
        if gname in self.gene_by_name:
            del self.gene_by_name[gname]

    def __setitem__(self, (gname, chrom), gene):
        if chrom is None:
            self.gene_by_name[gname] = gene
        else:
            super(GeneDict, self).__setitem__((gname, chrom), gene)
            if gname in self.gene_by_name:
                del self.gene_by_name[gname]

    def __contains__(self, (gname, chrom)):
        return super(GeneDict, self).__contains__((gname, chrom)) or gname in self.gene_by_name

    def __getitem__(self, (gname, chrom)):
        return super(GeneDict, self).get((gname, chrom)) or self.gene_by_name[gname]

    def get(self, (gname, chrom), *args, **kwargs):
        return super(GeneDict, self).get((gname, chrom), *args, **kwargs) or self.gene_by_name.get(gname, *args, **kwargs)

    def keys(self):
        return super(GeneDict, self).keys() + self.gene_by_name.keys()

    def values(self):
        return super(GeneDict, self).values() + self.gene_by_name.values()

    def items(self):
        return super(GeneDict, self).items() + self.gene_by_name.items()

    def __len__(self):
        return len(super(GeneDict, self).items()) + len(self.gene_by_name)


class ClinicalExperimentInfo:
    def __init__(self, cnf, sample, key_genes, target_type,
                 bed_fpath, mutations_fpath, sv_fpath, sv_vcf_fpath, varqc_json_fpath,
                 project_report_path, project_name, seq2c_tsv_fpath=None, targqc_report_path=None,
                 is_target2wgs_comparison=False):
        self.cnf = cnf
        self.sample = sample
        self.project_report_path = project_report_path
        self.project_name = project_name
        self.key_gene_by_name = GeneDict()
        self.key_gene_by_name_chrom = GeneDict()
        self.key_genes_number = None
        self.genes_collection_type = ''
        self.genes_description = ''
        self.key = ''
        self.patient = Patient()
        self.target = Target(type_=target_type, bed_fpath=bed_fpath, targqc_link=targqc_report_path)
        self.ave_depth = None
        self.depth_cutoff = None
        self.region_depth_cutoff = None
        self.actionable_genes_dict = None
        self.total_variants = None
        self.mutations = None
        self.sv_events = None
        self.seq2c_events_by_gene = None
        self.sv_fpath = sv_fpath
        self.sv_vcf_fpath = sv_vcf_fpath
        self.is_target2wgs_comparison = is_target2wgs_comparison

        info('Sample: ' + str(sample.name))
        info('Match sample name: ' + str(sample.normal_match))
        info()

        if not is_target2wgs_comparison:  # use all genes from bed instead of key genes if bed exists and number of genes < 2000
            key_gene_names_chroms, use_custom_panel = get_key_or_target_bed_genes(bed_fpath, key_genes)
        else:
            use_custom_panel = False
            key_gene_names_chroms, _ = get_key_or_target_bed_genes(None, key_genes)

        if use_custom_panel:
            self.genes_collection_type = 'target'
            self.genes_description = 'target genes'
            info('Preparing data for a clinical report for ' + str(len(key_gene_names_chroms)) + ' target genes from ' + str(bed_fpath) + ', sample ' + self.sample.name)
        else:
            self.genes_collection_type = 'key'
            self.genes_description = 'genes that have been previously implicated in various cancers'
            info('Preparing data for a clinical report for AZ 300 key genes ' + str(key_genes) + ', sample ' + self.sample.name)
        info()

        for gene_name, chrom in key_gene_names_chroms:
            self.key_gene_by_name_chrom[(gene_name, chrom)] = KeyGene(gene_name, chrom=chrom)
        self.key_genes_number = len(self.key_gene_by_name_chrom)

        if self.sample.targqc_dirpath and verify_dir(self.sample.targqc_dirpath) \
                and self.sample.targetcov_json_fpath and verify_file(self.sample.targetcov_json_fpath):
            info('Parsing target and patient info from ' + str(self.sample.targetcov_json_fpath))
            self.patient.gender = get_gender(self.sample, self.sample.targetcov_json_fpath)
            self.target.coverage_percent = get_target_fraction(self.sample, self.sample.targetcov_json_fpath)
            info('Parsing TargQC ' + self.genes_collection_type + ' genes stats...')
            self.ave_depth = get_ave_coverage(self.sample, self.sample.targetcov_json_fpath)
            self.depth_cutoff = int(self.ave_depth / 2)
            self.region_depth_cutoff = get_depth_cutoff(self.ave_depth, self.cnf.coverage_reports.depth_thresholds)
            self.sample.targetcov_detailed_tsv = verify_file(self.sample.targetcov_detailed_tsv)
            if self.sample.targetcov_detailed_tsv:
                self.parse_targetseq_detailed_report()
        else:
            warn('No targetcov_json_fpath provided, skipping key genes coverage stats.')
        info()

        info('Parsing actionable genes...')
        self.actionable_genes_dict = parse_broad_actionable()

        if mutations_fpath and verify_file(mutations_fpath):
            info('Parsing mutations from ' + str(mutations_fpath))
            if varqc_json_fpath and verify_file(varqc_json_fpath):
                self.total_variants = get_total_variants_number(self.sample, varqc_json_fpath)
            self.mutations = parse_mutations(self.cnf, self.sample, self.key_gene_by_name_chrom, mutations_fpath, self.genes_collection_type)
            for mut in self.mutations:
                if mut.gene.key in self.key_gene_by_name_chrom:
                    self.key_gene_by_name_chrom[mut.gene.key].mutations.append(mut)
            info('Retrieving SolveBio...')
            self.get_mut_info_from_solvebio()
        else:
            warn('No mutations_fpath provided, skipping mutation stats.')
        info()

        if sv_fpath and verify_file(sv_fpath):
            info('Parsing prioritized SV from ' + str(sv_fpath))
            self.sv_events = self.parse_sv(sv_fpath, self.key_gene_by_name_chrom)
        info()

        if seq2c_tsv_fpath and verify_file(seq2c_tsv_fpath):
            info('Parsing Seq2C from ' + str(seq2c_tsv_fpath))
            self.seq2c_events_by_gene = self.parse_seq2c_report(seq2c_tsv_fpath, self.key_gene_by_name_chrom, self.genes_collection_type)
        else:
            warn('No Seq2C results provided by option --seq2c, skipping plotting Seq2C')
        info()
        info('------')
        info('Done parsing data.')

    def __hash__(self):
        return hash((self.sample.name, self.project_name))

    def get_mut_info_from_solvebio(self):
        query_mutations(self.cnf, self.mutations)

    def parse_sv(self, sv_fpath, key_gene_by_name_chrom):
        info('Parsing prioritized SV events from ' + sv_fpath)
        sv_events = set()
        sv_events_by_gene_name = OrderedDict()

        sorted_known_fusions = [sorted(p) for p in known_fusions['homo_sapiens']]

        chr_order = get_chrom_order(self.cnf)

        with open(sv_fpath) as f:
            header_rows = []
            for i, l in enumerate(f):
                fs = l.strip().split('\t')
                if i == 0:
                    header_rows = fs  # caller  sample  chrom  start  end  svtype  known  end_gene  lof  annotation  split_read_support  paired_end_support
                else:
                    event = SVEvent.parse_sv_event(chr_order, **dict(zip(header_rows, fs)))
                    if event and event.sample == self.sample.name:
                        for annotation in event.annotations:
                            if event.is_fusion() and sorted(annotation.genes) in sorted_known_fusions:
                                info('Found ' + '/'.join(annotation.genes) + ' in known')
                                annotation.known = True

                            for g in annotation.genes:
                                if (g, event.chrom) in key_gene_by_name_chrom:
                                    sv_events_by_gene_name[(g, event.chrom)] = event
                                    sv_events.add(event)
                                    key_gene_by_name_chrom[(g, event.chrom)].sv_events.add(event)
                                    event.key_annotations.add(annotation)
                                    annotation.event = event
        return sv_events

    def parse_seq2c_report(self, seq2c_tsv_fpath, key_gene_by_name_chrom, genes_collection_type):
        seq2c_events_by_gene = dict()
        if not verify_file(seq2c_tsv_fpath, silent=True):
            return None

        with open(seq2c_tsv_fpath) as f_inp:
            for i, l in enumerate(f_inp):
                if i == 0: continue
                fs = l.replace('\n', '').split('\t')
                sname, gname = fs[0], fs[1]
                #if gname not in self.key_gene_by_name: continue
                if sname != self.sample.name: continue
                if 'not_a_gene' in gname: continue

                sname, gname, chrom, start, end, length, log2r, sig, fragment, amp_del, ab_seg, total_seg, \
                    ab_log2r, log2r_diff, ab_seg_loc, ab_samples, ab_samples_pcnt = fs[:17]
                if '_' in chrom:
                    continue
                gene = KeyGene(gname, chrom=chrom)
                gene.chrom, gene.start, gene.end = chrom, int(start), int(end)
                event = Seq2CEvent(
                    gene=gene,
                    fragment=fragment or None,
                    ab_log2r=float(ab_log2r) if ab_log2r else None,
                    log2r=float(log2r) if log2r else None,
                    amp_del=amp_del or None)
                seq2c_events_by_gene[gene] = event

        for gn, event in seq2c_events_by_gene.items():
            if gn.key in self.key_gene_by_name_chrom:
                self.key_gene_by_name_chrom[gn.key].seq2c_events.append(event)
        key_gene_events = sum(1 for gk in self.key_gene_by_name_chrom for e in self.key_gene_by_name_chrom[gk].seq2c_events)

        info('Found ' + str(len(seq2c_events_by_gene.values())) + ' Seq2C events (' + str(key_gene_events) +
             ' in ' + str(len(key_gene_by_name_chrom)) + ' ' + genes_collection_type + ' genes), ' +
             str(sum(1 for e in seq2c_events_by_gene.values() if e.is_amp())) + ' amplifications and ' +
             str(sum(1 for e in seq2c_events_by_gene.values() if e.is_del())) + ' deletions.')

        return seq2c_events_by_gene

    def parse_targetseq_detailed_report(self):
        info('Preparing coverage stats for ' + str(len(self.key_gene_by_name_chrom)) + ' ' + self.genes_collection_type + ' genes')
        if not verify_file(self.sample.targetcov_detailed_tsv):
            return None

        info('Parsing coverage stats from ' + self.sample.targetcov_detailed_tsv)
        depth_thresholds = self.cnf.coverage_reports.depth_thresholds
        with open(self.sample.targetcov_detailed_tsv) as f_inp:
            for l in f_inp:
                if l.startswith('#'):
                    def filter_digits(s):
                        return ''.join(c for c in s if c.isdigit())
                    fs = l.split('\t')
                    if len(fs) > 13:
                        depth_thresholds = [int(filter_digits(d)) for d in fs[13:]]
                    continue

                fs = l.split('\t')  # Chr	Start	End	Size	Gene	Strand	Feature	Biotype	TranscriptID    Min depth	Ave depth	Std dev	W/n 20% of ave depth	1x	5x	10x	25x	50x	100x	500x	1000x	5000x	10000x	50000x
                chrom, start, end, size, symbol, strand, feature, biotype, transcript_id, min_depth, ave_depth, std_dev, wn20pcnt = fs[:13]
                pcnt_val_by_thresh = fs[13:]

                symbol = get_val(symbol)
                if (symbol, chrom) not in self.key_gene_by_name_chrom:
                    continue

                chrom = get_val(chrom)
                if not chrom:
                    err('For gene ' + str(symbol) + ', chrom value is empty in line ' + l + ' in ' + self.sample.targetcov_detailed_tsv)
                    continue

                if start == '.' or end == '.':
                    continue

                start, end = int(start), int(end)
                ave_depth = get_float_val(ave_depth)
                cov_by_threshs = dict((t, get_float_val(f)) for t, f in izip(depth_thresholds, pcnt_val_by_thresh))

                if feature in ['Whole-Gene', 'Gene-Exon']:
                    gene = self.key_gene_by_name_chrom.get((symbol, chrom))
                    if gene:
                        gene.chrom = chrom
                        gene.start = start
                        assert gene.start
                        gene.end = end
                        assert gene.end
                        gene.ave_depth = ave_depth
                        gene.cov_by_threshs = cov_by_threshs
                        gene.transcript_id = transcript_id
                        self.key_gene_by_name_chrom[(symbol, chrom)] = gene

                elif feature in ['CDS', 'Exon']:
                    cds = CDS()
                    cds.start = start
                    cds.end = end
                    cds.ave_depth = ave_depth
                    cds.cov_by_threshs = cov_by_threshs
                    gene = self.key_gene_by_name_chrom.get((symbol, chrom))
                    if gene:
                        gene.cdss.append(cds)
                        gene.strand = strand

                else:
                    gene = self.key_gene_by_name_chrom.get((symbol, chrom))
                    if gene:
                        if gene.chrom is None:
                            gene.chrom = chrom
                        if gene.start is None:
                            gene.start = start
                        if gene.end is None:
                            gene.end = end

        # Cleaning up records that are not found in the target gene panel,
        # so we don't know about them and don't even want to report them
        # for gene in self.key_gene_by_name_chrom.values():
            # if not gene.chrom and (gene.name, gene.chrom) in self.key_gene_by_name_chrom:
            #     del self.key_gene_by_name_chrom[(gene.name, gene.chrom)]
            # if not gene.start or not gene.end:
            #     del self.key_gene_by_name_chrom[(gene.name, gene.chrom)]
        info('Keeping ' + str(len(self.key_gene_by_name_chrom)) + ' ' + self.genes_collection_type + ' genes based on targQC reports')

def parse_mutations(cnf, sample, key_gene_by_name_chrom, mutations_fpath, key_collection_type, for_flagged_report=False):
    mutations = []
    if for_flagged_report:
        info('Preparing mutations stats for flagged regions report')
    else:
        info('Preparing mutations stats for ' + key_collection_type + ' gene tables')
    info('Checking ' + mutations_fpath)
    if not verify_file(mutations_fpath):
        mut_pass_ending = source.mut_pass_suffix + '.' + source.mut_file_ext
        mut_basename = mutations_fpath.split('.' + mut_pass_ending)[0]
        if sample.normal_match:
            mutations_fpath = mut_basename + '.' + source.mut_paired_suffix + '.' + mut_pass_ending
        else:
            mutations_fpath = mut_basename + '.' + source.mut_single_suffix + '.' + mut_pass_ending
        info('Checking ' + mutations_fpath)
        if not verify_file(mutations_fpath):
            err('Cannot find PASSed mutations fpath')
            return []

    with open(cnf.canonical_transcripts) as f:
        canonical_transcripts = [tr.strip() for tr in f]

    info('Reading mutations from ' + mutations_fpath)
    alts_met_before = set()
    sample_col = None
    chr_col = None
    pos_col = None
    ref_col = None
    alt_col = None
    class_col = None
    type_col = None
    allele_freq_col = None
    gene_col = None
    depth_col = None
    codon_chg_col = None
    aa_len_col = None
    aa_chg_col = None
    cdna_chg_col = None
    transcript_col = None
    exon_col = None
    status_col = None
    signif_col = None
    reason_col = None
    ids_col = None
    var_type_col = None

    chr_order = get_chrom_order(cnf)

    with open(mutations_fpath) as f:
        for i, l in enumerate(f):
            if i == 0:
                header = l.strip().split('\t')
                sample_col = header.index('Sample')
                chr_col = header.index('Chr')
                ids_col = header.index('ID')
                pos_col = header.index('Start')
                ref_col = header.index('Ref')
                alt_col = header.index('Alt')
                class_col = header.index('Var_Class')
                type_col = header.index('Type')
                var_type_col = header.index('Var_Type')
                codon_chg_col = header.index('Codon_Change')
                aa_chg_col = header.index('Amino_Acid_Change')
                cdna_chg_col = header.index('cDNA_Change')
                aa_len_col = header.index('Amino_Acid_Length')
                allele_freq_col = header.index('AlleleFreq')
                gene_col = header.index('Gene')
                depth_col = header.index('Depth')
                transcript_col = header.index('Transcript')
                exon_col = header.index('Exon')
                if 'Status' in header:
                    status_col = header.index('Status')
                if 'Significance' in header:
                    signif_col = header.index('Significance')
                else:
                    signif_col = len(header) - header[::-1].index('Status') - 1  # get last index of status
                if 'Reason' in header:
                    reason_col = header.index('Reason')
                continue
            fs = l.replace('\n', '').split('\t')
            sample_name, chrom, start, ref, alt, gname, transcript = \
                fs[sample_col], fs[chr_col], fs[pos_col], fs[ref_col], fs[alt_col], fs[gene_col], fs[transcript_col]
            codon_change, cdna_change, aa_change, aa_len, exon = \
                fs[codon_chg_col], fs[cdna_chg_col], fs[aa_chg_col], fs[aa_len_col], fs[exon_col]
            ids, type_, var_type, var_class = fs[ids_col], fs[type_col], fs[var_type_col], fs[class_col]
            depth, af = fs[depth_col], fs[allele_freq_col]
            status = fs[status_col] if status_col is not None else None
            signif = fs[signif_col] if signif_col is not None else None
            reason = fs[reason_col] if reason_col is not None else None

            if sample_name == sample.name:
                if (chrom, start, ref, alt, transcript) in alts_met_before:
                    continue
                alts_met_before.add((chrom, start, ref, alt, transcript))

                if (gname, chrom) not in key_gene_by_name_chrom:
                    # warn('gene ' + gname + ' at ' + chrom + ' not found in coverage reports, but found in mutation:\n  ' + l)
                    continue

                mut = Mutation(chrom=chrom, chrom_ref_order=chr_order.get(chrom))
                mut.gene = KeyGene(gname, chrom=chrom)
                mut.transcript = transcript
                mut.is_canonical = transcript in canonical_transcripts if canonical_transcripts else True
                mut.codon_change = codon_change
                mut.cdna_change = cdna_change
                mut.aa_change = aa_change
                mut.aa_len = aa_len
                mut.exon = exon
                mut.pos = int(start)
                mut.ref = ref
                mut.alt = alt
                if depth:
                    mut.depth = int(depth)
                else:
                    mut.depth = None
                if af:
                    mut.freq = float(af)
                else:
                    mut.freq = None
                mut.dbsnp_id = next((id_.split(',')[0] for id_ in ids.split(';') if id_.startswith('rs')), None)
                mut.cosmic_id = next((id_.split(',')[0] for id_ in ids.split(';') if id_.startswith('COS')), None)
                mut.eff_type = (type_[0] + type_[1:].lower().replace('_', ' ')) if type_ else type_
                mut.var_type = var_type
                mut.var_class = var_class
                mut.status = status
                mut.signif = signif
                if reason:
                    reason = reason.replace('_', ' ')
                    if reason == 'actionable somatic':
                        reason = 'actionable som.'
                    if reason == 'actionable germline':
                        reason = 'actionable germ.'
                    reason = reason.replace('change', 'chg.')
                mut.reason = reason

                mutations.append(mut)

    info('Found ' + str(len(mutations)) + ' mutations in ' + str(len(key_gene_by_name_chrom)) + ' ' + key_collection_type + ' genes')
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
    r = sr.find_record(sr.records, 'Sex') or sr.find_record(sr.records, 'Gender')

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


def get_key_or_target_bed_genes(bed_fpath, key_genes):
    use_custom_panel = False
    key_gene_names_chroms = None
    if bed_fpath:
        key_gene_names_chroms, gene_names_list = get_gene_keys(bed_fpath)
        if len(key_gene_names_chroms) < 2000:
            use_custom_panel = True
    if not use_custom_panel:
        key_gene_names = get_key_genes(key_genes)
        key_gene_names_chroms = [(gn, None) for gn in key_gene_names]
    key_gene_names_chroms = [(gn, c) for (gn, c) in key_gene_names_chroms if (gn and gn != '.')]
    return key_gene_names_chroms, use_custom_panel
