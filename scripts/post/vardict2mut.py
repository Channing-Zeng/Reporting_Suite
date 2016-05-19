#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc

from collections import defaultdict, OrderedDict
from optparse import OptionParser
from os.path import join
from os.path import exists
import time
import sys
import re
import tabix

from source import verify_file
from source.config import Config, defaults
from source import logger
from source.file_utils import adjust_path, verify_dir
from source.logger import info, critical, err, warn, debug
from source.prepare_args_and_cnf import determine_run_cnf, check_genome_resources, \
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug
from source.prepare_args_and_cnf import determine_sys_cnf
from source.utils import OrderedDefaultDict


def get_args():
    info(' '.join(sys.argv))
    info()
    description = (
        'The program will filter the VarDict output after vcf2txt.pl to '
        'candidate interpretable mutations, somatic or germline.')
    parser = OptionParser(description=description)
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser, threads=1)

    parser.add_option('-o', dest='output_file')
    parser.add_option('--o-all-transcripts', dest='all_transcripts_output_file')
    parser.add_option('--o-fm', dest='fm_output_file')

    parser.add_option('--cohort-freqs', dest='cohort_freqs_fpath')

    parser.add_option('-D', '--min-depth', dest='filt_depth', type='int', help='The minimum total depth')
    parser.add_option('-V', '--min-vd', dest='min_vd', type='int', help='The minimum reads supporting variant')
    parser.add_option('--gmaf', dest='min_gmaf', type='float',
                      help='When the GMAF is greater than specified, it\'s considered common SNP and filtered out.')
    parser.add_option('-f', '--min-freq', dest='min_freq', type='float',
                      help='The minimum allele frequency for regular variants. Default: 0.05')
    parser.add_option('-F', '--min-freq-hs', dest='min_hotspot_freq', type='float',
                      help='The minimum allele frequency hotspot somatic mutations, typically lower then -f. '
                           'Default: 0.01 or half -f, whichever is less')
    parser.add_option('-N', '--keep-utr-intronic', dest='keep_utr_intronic', action='store_true',
                      help='Keep all intronic and UTR in the output, but will be set as "unknown".')

    parser.add_option('-p', '--platform', dest='platform',
                      help='The platform, such as WXS, WGS, RNA-Seq, VALIDATION, etc. No Default. '
                           'Used for output in FM\'s format')

    parser.set_usage('Usage: ' + __file__ + ' vcf2txt_res_fpath [opts] -o output_fpath')

    (opts, args) = parser.parse_args()
    if len(args) < 1:
        critical('Provide the first argument - output from vcf2txt.pl')
    logger.is_debug = opts.debug

    vcf2txt_res_fpath = verify_file(args[0], is_critical=True)

    run_cnf = determine_run_cnf(opts)
    cnf = Config(opts.__dict__, determine_sys_cnf(opts), run_cnf)
    if not cnf.genome:
        critical('Please, specify the --genome option (e.g. --genome hg19)')

    check_genome_resources(cnf)

    if not cnf.output_file:
        critical('Please, specify the output fpath with -o')

    info()

    return cnf, vcf2txt_res_fpath


def main():
    cnf, vcf2txt_res_fpath = get_args()

    info('-' * 70)
    info('Writing to ' + cnf.output_file)
    if cnf.all_transcripts_output_file:
        info('Writing info for all transcripts to ' + cnf.all_transcripts_output_file)
    if cnf.fm_output_file:
        info('Writing in FM format to ' + cnf.fm_output_file)

    f = Filtration(cnf)

    input_f = open(verify_file(vcf2txt_res_fpath))
    output_f = open(adjust_path(cnf.output_file), 'w')
    fm_output_f = open(adjust_path(cnf.fm_output_file), 'w') if cnf.fm_output_file else None
    all_transcripts_output_f = open(adjust_path(cnf.all_transcripts_output_file), 'w') if cnf.all_transcripts_output_file else None

    info()
    info('-' * 70)
    info('Running filtering...')
    f.do_filtering(input_f, output_f, fm_output_f, all_transcripts_output_f)

    input_f.close()
    output_f.close()
    if fm_output_f:
        fm_output_f.close()
    if all_transcripts_output_f:
        all_transcripts_output_f.close()

    info()
    info('Saved to ' + cnf.output_file)


def iter_lines(fpath):
    with open(fpath) as f:
        for l in f:
            l = l.replace('\n', '')
            if not l or l.startswith('#'):
                continue
            yield l


def parse_gene_blacklists(cnf):
    _d = OrderedDict()
    if cnf.variant_filtering.blacklist.genes.published:
        _d['freq mut gene in HGMD'] = 'published/flags_in_hgmd.txt'
        _d['freq mut gene in OMIM'] = 'published/flags_in_omim.txt'
        _d['freq mut gene'] = 'published/flags.txt'
        _d['incidentalome gene'] = 'published/incidentalome.txt'
        _d['mutSigCV gene'] = 'published/mutsigcv.txt'
    if cnf.variant_filtering.blacklist.genes.low_complexity:
        _d['low complexity gene'] = 'low_complexity/low_complexity_entire_gene.txt'
    if cnf.variant_filtering.blacklist.genes.repetitive_single_exome:
        _d['repetitive single exon gene'] = 'low_complexity/repetitive_single_exon_gene.txt'
    if cnf.variant_filtering.blacklist.genes.abnormal_gc:
        _d['low GC gene'] = 'low_complexity/low_gc.txt'
        _d['high GC gene'] = 'low_complexity/high_gc.txt'
    if cnf.variant_filtering.blacklist.genes.too_many_cosmic_mutations:
        _d['gene with too many COSMIC mutations'] = 'low_complexity/too_many_cosmic_mutations.txt'

    d = OrderedDefaultDict(dict)
    for reason, fn in _d.items():
        fpath = verify_file(join(cnf.crapomedir, fn), description=reason + ' blacklist genes file', is_critical=True)
        for l in iter_lines(fpath):
            fs = l.split('\t')
            gene_name = l.split('\t')[0]
            meta_info = l.split('\t')[1] if len(fs) == 2 else ''
            d[reason][gene_name] = meta_info
    return d


def load_region_blacklists(cnf):
    _d = OrderedDict()
    if cnf.variant_filtering.blacklist.regions.low_complexity:
        _d['low complexity region'] = 'low_complexity.bed.gz'
    if cnf.variant_filtering.blacklist.regions.abnormal_gc:
        _d['low GC region'] = 'low_gc.bed.gz'
        _d['high GC region'] = 'high_gc.bed.gz'
    if cnf.variant_filtering.blacklist.regions.hengs_universal_mask:
        _d['Hengs mask'] = 'heng_um75-hs37d5.bed.gz'
    if cnf.variant_filtering.blacklist.regions.repeats:
        _d['repetitive region'] = 'repeats.bed.gz'

    d = OrderedDict()
    for reason, fn in _d.items():
        fpath = verify_file(join(cnf.genome.tricky_regions, fn), description=reason + ' tricky regions file', is_critical=True)
        # d[reason] = build_interval_tree(fpath)
        d[reason] = tabix.open(fpath)
    return d


# def build_interval_tree(bed_fpath):
#     # dictionary mapping chromosome names to interval trees
#     info('Building interval tree for ' + bed_fpath + ': ', ending='')
#     genome = defaultdict(IntervalTree)
#
#     # parse the BED file and build the interval trees
#     for l in iter_lines(bed_fpath):
#         fs = l.split('\t')[:3]
#         chrom = fs[0]
#         if chrom not in genome:
#             info(chrom, ending=' ', print_date=False)
#         start = int(fs[1])
#         end = int(fs[2])
#         genome[chrom].addi(start, end)
#     info('', print_date=False)
#
#     return genome


class Filtration:
    statuses = ['', 'known', 'likely', 'unknown', 'incidentalome']  # Tier 1, 2, 3, 4
    sensitization_aa_changes = {'EGFR-T790M': 'TKI'}

    def __init__(self, cnf):
        self.canonical_reject_counter = OrderedDefaultDict(int)
        self.no_transcript_reject_counter = OrderedDefaultDict(int)
        self.all_reject_counter = OrderedDefaultDict(int)

        self.canonical_counter = OrderedDefaultDict(int)
        self.no_transcript_counter = OrderedDefaultDict(int)
        self.all_counter = OrderedDefaultDict(int)

        self.canonical_blacklist_counter = OrderedDefaultDict(int)
        self.no_transcript_blacklist_counter = OrderedDefaultDict(int)
        self.gene_blacklist_counter = OrderedDefaultDict(int)
        self.region_blacklist_counter = OrderedDefaultDict(int)

        cnf.genome.compendia_ms7_hotspot    = verify_file(cnf.genome.compendia_ms7_hotspot, 'compendia_ms7_hotspot')
        cnf.genome.actionable               = verify_file(cnf.genome.actionable, 'actionable')
        cnf.genome.filter_common_snp        = verify_file(cnf.genome.filter_common_snp, 'filter_common_snp')
        cnf.genome.filter_common_artifacts  = verify_file(cnf.genome.filter_common_artifacts, 'filter_common_artifacts')
        cnf.genome.splice                   = verify_file(cnf.genome.splice, 'splice')
        cnf.suppressors                     = verify_file(cnf.suppressors, 'suppressors')
        cnf.oncogenes                       = verify_file(cnf.oncogenes, 'oncogenes')
        cnf.ruledir                         = verify_dir(cnf.ruledir, 'ruledir')
        cnf.snpeffect_export_polymorphic    = verify_file(cnf.snpeffect_export_polymorphic, 'snpeffect_export_polymorphic')
        cnf.actionable_hotspot              = verify_file(cnf.actionable_hotspot, 'actionable_hotspot')
        cnf.specific_mutations              = verify_file(cnf.specific_mutations, 'specific_mutations')
        cnf.last_critical_aa                = verify_file(cnf.last_critical_aa, 'last_critical_aa')
        cnf.crapomedir                      = verify_dir(cnf.crapomedir, 'crapomedir')
        if not all([cnf.genome.compendia_ms7_hotspot,
                    cnf.genome.actionable,
                    cnf.genome.filter_common_snp,
                    cnf.genome.filter_common_artifacts,
                    cnf.genome.splice,
                    cnf.suppressors,
                    cnf.oncogenes,
                    cnf.ruledir,
                    cnf.snpeffect_export_polymorphic,
                    cnf.actionable_hotspot,
                    cnf.specific_mutations,
                    cnf.last_critical_aa,
                    cnf.crapomedir
                    ]):
            critical('Critical: some of the required files are not found or empty (see above)')

        self.suppressors = parse_genes_list(adjust_path(cnf.suppressors))
        self.oncogenes = parse_genes_list(adjust_path(cnf.oncogenes))

        self.reg_exp_sample = cnf.reg_exp_sample
        self.platform = cnf.platform

        canon_tr_fpath = verify_file(cnf.canonical_transcripts, is_critical=True)
        info('Using canonical transcripts from ' + canon_tr_fpath)
        with open(canon_tr_fpath) as f:
            self.canonical_transcripts = [tr.strip().split('.')[0] for tr in f]

        self.max_ratio = cnf.max_ratio or cnf.variant_filtering.max_ratio
        self.max_sample_cnt = cnf.max_sample_cnt or cnf.variant_filtering.max_sample_cnt
        self.min_freq = cnf.min_freq or cnf.variant_filtering.min_freq_vardict2mut
        self.min_hotspot_freq = cnf.min_hotspot_freq or cnf.variant_filtering.min_hotspot_freq
        if self.min_hotspot_freq is None or self.min_hotspot_freq == 'default':
            self.min_hotspot_freq = min(0.01, self.min_freq / 2)
        self.filt_depth = cnf.variant_filtering.filt_depth
        self.min_vd = cnf.variant_filtering.min_vd
        self.min_gmaf = cnf.variant_filtering.min_gmaf

        # self.freq_in_sample_by_vark = dict()
        # if cnf.cohort_freqs_fpath:
        #     cohort_freqs_fpath = verify_file(cnf.cohort_freqs_fpath, is_critical=True)
        #     with open(cohort_freqs_fpath) as f:
        #         for l in f:
        #             fs = l.replace('\n', '').split()
        #             self.freq_in_sample_by_vark[fs[0]] = float(fs[1])

        info('Parsing filtering data...')
        self.tp53_groups = {'Group 1': parse_mut_tp53(join(cnf.ruledir, 'Rules', 'DNE.txt')),
                            'Group 2': parse_mut_tp53(join(cnf.ruledir, 'Rules', 'TA0-25.txt')),
                            'Group 3': parse_mut_tp53(join(cnf.ruledir, 'Rules', 'TA25-50_SOM_10x.txt'))}

        self.splice_positions_by_gene = defaultdict(set)
        for l in iter_lines(cnf.genome.splice):
            pos, g = l.split('\t')
            self.splice_positions_by_gene[g].add(pos)

        self.last_critical_aa_pos_by_gene = dict()
        for l in iter_lines(cnf.last_critical_aa):
            g, aa_pos, transcript = l.split('\t')
            self.last_critical_aa_pos_by_gene[g] = int(aa_pos)

        self.filter_snp = set()
        for l in iter_lines(cnf.genome.filter_common_snp):
            fields = l.split('\t')
            self.filter_snp.add('-'.join(fields[1:5]))

        self.snpeff_snp = set()
        self.snpeff_snp_rsids = set()
        for l in iter_lines(cnf.snpeffect_export_polymorphic):
            fields = l.split('\t')
            snpeff_aachg = fields[2]
            snpeff_rsid = fields[5]
            if len(fields) > 11 and fields[11]:
                snpeff_gene = fields[11]
                self.snpeff_snp.add('-'.join([snpeff_gene, snpeff_aachg]))
            elif snpeff_rsid != '-':
                self.snpeff_snp_rsids.add(snpeff_rsid)

        self.filter_artifacts = set()
        for l in iter_lines(cnf.genome.filter_common_artifacts):
            fields = l.split('\t')
            self.filter_artifacts.add('-'.join(fields[1:5]))

        self.actionable_hotspots = defaultdict(set)
        for l in iter_lines(cnf.actionable_hotspot):
            fields = l.split('\t')
            self.actionable_hotspots[fields[0]].add(fields[1])

        self.act_somatic = dict()
        self.act_germline = set()

        self.rules = defaultdict(lambda: defaultdict(list))
        # inframe_del = 'inframe-del'
        # inframe_ins = 'inframe-ins'
        # self.rules[inframe_del] = {}
        # self.rules[inframe_ins] = {}
        for l in iter_lines(cnf.genome.actionable):
            fields = l.split('\t')
            if fields[7] == 'germline':
                key = '-'.join(fields[1:5])
                self.act_germline.add(key)
            elif fields[7] == 'somatic':
                change = fields[8].strip()
                if fields[6] == 'rule':
                    if fields[4] == '*' and len(fields[3]) == 1:
                        key = '-'.join(fields[1:4])
                        self.act_somatic[key] = change
                    else:
                        indel_type = fields[5]
                        gene = fields[0]
                        chrom = fields[1]
                        start = int(fields[2])
                        end = int(fields[3])
                        n = int(fields[4])
                        self.rules[indel_type][gene].append([chrom, start, end, n, change])
                    # elif fields[5] == inframe_del:
                    #     self.rules[inframe_del].setdefault(fields[0], []).append([fields[1]] + [int (f) for f in fields[2:5]])
                    # elif fields[5] == inframe_ins:
                    #     self.rules[inframe_ins].setdefault(fields[0], []).append([fields[1]] + [int (f) for f in fields[2:5]])

                else:
                    key = '-'.join(fields[1:5])
                    self.act_somatic[key] = change

        self.hotspot_nucleotides = set()
        self.hotspot_proteins = set()
        for l in iter_lines(cnf.genome.compendia_ms7_hotspot):
            fields = l.split('\t')
            if fields[5].startswith('g.'):
                continue
            self.hotspot_nucleotides.add('-'.join(fields[1:5]))
            if not fields[6]:
                continue
            self.hotspot_proteins.add('-'.join([fields[0], fields[6]]))

        info('Parsing gene blacklists...')
        self.gene_blacklists_by_reason = parse_gene_blacklists(cnf)
        for r in self.gene_blacklists_by_reason.keys():
            self.gene_blacklist_counter[r] = 0
        self.gene_to_soft_filter = list(iter_lines(join(cnf.crapomedir, 'soft_filter.txt')))

        info('Parsing region blacklists...')
        self.region_blacklists_by_reason = load_region_blacklists(cnf)
        for r in self.region_blacklists_by_reason.keys():
            self.region_blacklist_counter[r] = 0

        info('Parsing spreadsheat with actionable rules...')
        self.tier_by_specific_mutations, \
        self.genes_with_generic_rules, \
        self.tier_by_type_by_region_by_gene, \
        self.sensitizations_by_gene, \
        self.spec_transcripts_by_aachg = parse_specific_mutations(cnf.specific_mutations)

        if not all([self.rules, self.splice_positions_by_gene, self.act_somatic, self.act_germline, self.actionable_hotspots]):
            if not self.rules:
                err('No rules, cannot proceed')
            if not self.splice_positions_by_gene:
                err('No tp53_positions, cannot proceed')
            if not self.act_somatic:
                err('No act_somatic, cannot proceed')
            if not self.act_germline:
                err('No act_germline, cannot proceed')
            if not self.actionable_hotspots:
                err('No actionable_hotspots, cannot proceed')

        self.status = None
        self.reason_by_status = None

    def update_status(self, new_status, new_reason, force=False):
        if isinstance(new_reason, list):
            for r in new_reason:
                self.reason_by_status[new_status].add(r)
        else:
            self.reason_by_status[new_status].add(new_reason)

        if not force and Filtration.statuses.index(new_status) > Filtration.statuses.index(self.status):
            return self.status
        else:
            self.status = new_status
        return self.status

    def check_by_var_class(self, var_class, cosmic_aachg, cosmic_counts):
        if var_class == 'ClnSNP':
            self.update_status('likely', 'clin_SNP')
        if var_class == 'dbSNP_del':
            self.update_status('likely', 'dbSNP_del')
        if var_class == 'ClnSNP_known':
            self.update_status('known', 'clin_SNP_known')
        if var_class == 'ClnSNP_unknown':
            self.update_status('unknown', 'clin_SNP_unknown')
        if var_class == 'COSMIC':
            if cosmic_counts:
                for c in cosmic_counts:
                    if c >= 5:
                        self.update_status('likely', 'COSMIC_5+')
            if cosmic_aachg is not None and cosmic_aachg.startswith('p.'):
                cosmic_aachg = cosmic_aachg[2:]
        return cosmic_aachg

    def check_by_type(self, var_type, aa_chg, cdna_chg, effect):
        if 'FRAME_SHIFT' in var_type or 'FRAMESHIFT' in var_type:
            self.update_status('likely', 'frame_shift')
        elif stop_gain_pattern.match(aa_chg) or 'STOP_GAIN' in var_type:
            self.update_status('likely', 'stop_gained')
        elif 'START_LOST' in var_type and effect == 'HIGH':
            self.update_status('likely', 'start_lost')

        if 'SPLICE' in var_type and ('ACCEPTOR' in var_type or 'SPLICE_DONOR' in var_type):
            self.update_status('likely', 'splice_site')
            aa_chg = 'splice'
        elif not aa_chg and 'SPLICE' in var_type and 'REGION_VARIANT' not in var_type:
            if cdna_chg:
                cdna_pos = None
                m = re.match('.*\d+\+(\d+)', cdna_chg).groups()
                if m:
                    cdna_pos = int(m[0])
                else:
                    m = re.match('.*\d+-(\d+)[^_]\S+$', cdna_chg).groups()
                    if m:
                        cdna_pos = int(m[0])
                if cdna_pos is not None and cdna_pos <= 2:
                    self.update_status('likely', 'splice_site')
                    aa_chg = 'splice'
            else:  # No cDNA_Change, for earlier version compatibility
                self.update_status('likely', 'splice_site')
                aa_chg = 'splice'
        return aa_chg

    aa_chg_trim_pattern = re.compile('^([A-Z]\d+)[A-Z]$')
    def check_actionable(self, chrom, pos, ref, alt, gene, aa_chg):
        key = '-'.join([chrom, pos, ref, alt])
        if key in self.act_somatic:
            return aa_chg, self.update_status('known', 'act_somatic')
        if key in self.act_germline:
            return aa_chg, self.update_status('known', 'act_germline')
        if len(ref) == 1 and len(ref) == len(alt):
            key = '-'.join([chrom, pos, ref])
            if key in self.act_somatic:
                return aa_chg, self.update_status('known', 'act_somatic')

        if Filtration.aa_chg_trim_pattern.match(aa_chg) and gene in self.actionable_hotspots:
            act_hotspots = self.actionable_hotspots[gene]
            if aa_chg in act_hotspots:
                return aa_chg, self.update_status('known', 'act_hotspot')
            aa_chg_trim = re.findall(Filtration.aa_chg_trim_pattern, aa_chg)[0]
            if aa_chg_trim in act_hotspots:
                return aa_chg, self.update_status('known', 'act_hotspot')

        if gene == 'TP53':
            tp53_group = classify_tp53(aa_chg, pos, ref, alt, self.splice_positions_by_gene[gene], self.tp53_groups)
            if tp53_group != 'NA':
                return aa_chg, self.update_status('known', 'act_somatic_tp53_group_' + tp53_group.split(' ')[1])

        if gene in self.rules['inframe-del'] and len(ref) > len(alt) and (len(ref) - len(alt)) % 3 == 0:
            for r in self.rules['inframe-del'][gene]:
                if r[0] == chrom and r[1] <= int(pos) <= r[2] and len(ref) - len(alt) >= r[3]:
                    aa_chg = r[4]
                    return aa_chg, self.update_status('known', 'act_somatic_inframe_del')

        elif gene in self.rules['inframe-ins'] and len(ref) < len(alt) and (len(alt) - len(ref)) % 3 == 0:
            for r in self.rules['inframe-ins'][gene]:
                if r[0] == chrom and r[1] <= int(pos) <= r[2] and len(alt) - len(ref) >= r[3]:
                    aa_chg = r[4]
                    return aa_chg, self.update_status('known', 'act_somatic_inframe_ins')

        elif gene in self.rules['indel'] and len(ref) != len(alt):
            for r in self.rules['indel'][gene]:
                if r[0] == chrom and r[1] <= int(pos) <= r[2] and len(alt) - len(ref) >= r[3]:
                    aa_chg = r[4]
                    return aa_chg, self.update_status('known', 'act_somatic_indel')

        elif gene in self.rules['del'] and len(ref) > len(alt):
            for r in self.rules['del'][gene]:
                if r[0] == chrom and r[1] <= int(pos) <= r[2] and len(ref) - len(alt) >= r[3]:
                    aa_chg = r[4]
                    return aa_chg, self.update_status('known', 'act_somatic_del')

        elif gene in self.rules['ins'] and len(ref) < len(alt):
            for r in self.rules['ins'][gene]:
                if r[0] == chrom and r[1] <= int(pos) <= r[2] and len(alt) - len(ref) >= r[3]:
                    aa_chg = r[4]
                    return aa_chg, self.update_status('known', 'act_somatic_ins')

        return aa_chg, False

    def check_rob_hedley_actionable(self, gene, aa_chg, effect, region, transcript):
        if aa_chg:
            gene_aachg = '-'.join([gene, aa_chg])
            if transcript and gene_aachg in self.spec_transcripts_by_aachg:
                if self.spec_transcripts_by_aachg[gene_aachg] != transcript:
                    return None
            if gene_aachg in self.tier_by_specific_mutations:
                tier = self.tier_by_specific_mutations[gene_aachg]
                self.update_status(Filtration.statuses[tier], 'actionable')
                return True

        if region and effect in ['HIGH', 'MODERATE']:
            codon = re.sub('[^0-9]', '', aa_chg)
            gene_codon_chg = '-'.join([gene, region, codon])
            if gene_codon_chg in self.tier_by_specific_mutations:
                tier = self.tier_by_specific_mutations[gene_codon_chg]
                # status, reasons = self.update_status(status, reasons, statuses[tier], 'manually_curated_codon_' + codon + '_in_exon_' + region)
                self.update_status(Filtration.statuses[tier], 'actionable_codon_' + codon + '_in_exon_' + region)
                return True

    def check_by_general_rules(self, var_type, aa_chg, is_lof, gene):
        # if 'splice_site' in reasons:
        #     status, reasons = self.update_status(status, reasons, 'known', ['general_rules'] + reasons)
        if is_lof:
            self.update_status('known', 'act_lof' + '_of_gene_' + gene)
            # status, reasons = self.update_status(status, reasons, 'known', ['lof_in_gene_' + gene])
        elif 'EXON_LOSS' in var_type or 'EXON_DELETED' in var_type:
            self.update_status('known', 'act_exon_loss' + '_in_gene_' + gene)
            # status, reasons = self.update_status(status, reasons, 'known', ['exon_loss_in_gene_' + gene])
        elif Filtration.statuses.index(self.status) <= 2:
            info(str(gene) + ' ' + str(aa_chg) + ' is in general rules, but does not alter protein function.'
                 ' Keeping status as ' + str(self.status))
        #     status, reasons = update_status(status, reasons, 'unlikely', 'but_not_alter_protein_function', force=True)

    def check_by_mut_type(self, cdna_chg, region, types_by_region, gene):
        if region in types_by_region:
            for type_ in types_by_region[region]:
                if type_ in cdna_chg:
                    tier = types_by_region[region][type_]
                    self.update_status(Filtration.statuses[tier], 'act_' + type_ + '_in_gene_' + gene)

    # def print_mutations_for_one_gene(self, output_f, fm_output_f, all_transcripts_output_f, cur_gene_mutations, gene, sensitizations):
    #     for cur_gene_mut_info in cur_gene_mutations:
    #         fields, status, reasons, gene_aachg, fm_data, is_canonical, no_transcript = cur_gene_mut_info
    #         for sensitization, sens_or_res in self.sensitizations_by_gene[gene]:
    #             sensitization_aa_chg = Filtration.sensitization_aa_changes.keys()[Filtration.sensitization_aa_changes.values().index(sensitization)]
    #
    #             if sensitization in sensitizations:
    #                 for s in Filtration.statuses:
    #                     self.reason_by_status[s].add(sensitization + '_' + sens_or_res)
    #                 # if status == 'likely':
    #                 # self.update_status('known', reasons + [aa_chg + '_' + sens_or_res], force=True)
    #                 # elif status in ['unlikely', 'unknown']:
    #                 #     self.update_status('likely', [aa_chg + '_' + sens_or_res], force=True)
    #
    #             # if sensitization not in sensitizations and status == 'known':
    #             #     self.update_status('likely', reasons + [sensitization_aa_chg + '_required'], force=True)
    #             # elif sensitization not in sensitizations and status == 'likely':
    #             #     self.update_status('unlikely', sensitization_aa_chg + '_required', force=True)
    #         self.print_mutation(output_f, fm_output_f, all_transcripts_output_f, status, reasons, fields, is_canonical, no_transcript, fm_data=fm_data)

    def print_mutation(self, output_f, fm_output_f, all_transcripts_output_f, status, reasons, fields, is_canonical, no_transcript, fm_data):
        self.apply_counter('lines_written', is_canonical, no_transcript)
        self.apply_counter(status, is_canonical, no_transcript)

        if is_canonical or no_transcript:
            if fm_data and fm_output_f:
                sample, platform, gene, pos, cosm_aa_chg, aa_chg, cdna_chg, chrom, depth, allele_freq = fm_data
                fm_output_f.write('\t'.join([sample, platform, 'short-variant', gene, status, aa_chg, cdna_chg,
                                       chrom + ':' + pos, str(depth), str(allele_freq * 100),
                                       '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-'])  + '\n')
            output_f.write('\t'.join(fields + [status]) + ('\t' + ', '.join(reasons) + '\n'))
        if all_transcripts_output_f:
            all_transcripts_output_f.write('\t'.join(fields + [status]) + ('\t' + ', '.join(reasons) + '\n'))
            # if status != fields[-2] or ','.join(reasons) != fields[-1]:
            #     out_f.write('\t'.join(fields[1:15] + fields[-3:] + [status]) + ('\t' + ','.join(reasons) + '\n'))

    def apply_counter(self, reason, is_canonical, no_transcript):
        self.all_counter[reason] += 1
        if is_canonical:
            self.canonical_counter[reason] += 1
        if no_transcript:
            self.no_transcript_counter[reason] += 1

    def apply_reject_counter(self, reason, is_canonical, no_transcript):
        self.all_reject_counter[reason] += 1
        if is_canonical:
            self.canonical_reject_counter[reason] += 1
        if no_transcript:
            self.no_transcript_reject_counter[reason] += 1

    def apply_gene_blacklist_counter(self, reason):
        self.gene_blacklist_counter[reason] += 1

    def apply_region_blacklist_counter(self, reason):
        self.region_blacklist_counter[reason] += 1

    def check_blacklist_genes(self, gene_name, aa_pos=None):
        reasons = []
        for reason, data in self.gene_blacklists_by_reason.items():
            if gene_name in data:
                meta_info = data[gene_name]
                if meta_info == '':
                    reasons.append(reason)
                elif aa_pos is not None:
                    fs = meta_info.split(':')  # regions in form of :232, 553:, 42:111
                    if fs[0] and aa_pos >= int(fs[0]):
                        reasons.append(reason)
                    if fs[1] and aa_pos < int(fs[1]):
                        reasons.append(reason)
        return reasons

    def check_blacklist_regions(self, chrom, start, end):
        reasons = []
        for reason, tb in self.region_blacklists_by_reason.items():
            try:
                records = list(tb.query(chrom, start, end))
            except tabix.TabixError:
                pass
            else:
                if records:
                    reasons.append(reason)
        return reasons

    def do_filtering(self, input_f, output_f, fm_output_f=None, all_transcripts_output_f=None):
        pass_col = None
        sample_col = None
        chr_col = None
        pos_col = None
        ref_col = None
        alt_col = None
        class_col = None
        type_col = None
        func_col = None
        gene_code_col = None
        allele_freq_col = None
        gene_col = None
        depth_col = None
        vd_col = None
        aa_chg_col = None
        cosmaachg_col = None
        cosmcnt_col = None
        msicol = None
        cdna_chg_col = None
        gene_coding_col = None
        transcript_col = None
        effect_col = None
        exon_col = None
        status_col = None
        reason_col = None
        lof_col = None

        headers = []

        sensitizations = []
        cur_gene_mutations = []
        prev_gene = ''

        # platform_regexp = re.compile('-\d\d[_-]([^_\d]+?)$')  # TODO: fix pattern
        platform_regexp = re.compile('[_-]([^_\d]+?)$')  # TODO: fix pattern
        sample_regexp = re.compile(self.reg_exp_sample) if self.reg_exp_sample else None

        if fm_output_f:
            fm_output_f.write('SAMPLE ID\tANALYSIS FILE LOCATION\tVARIANT-TYPE\tGENE\tSOMATIC STATUS/FUNCTIONAL IMPACT\tSV-PROTEIN-CHANGE\tSV-CDS-CHANGE\tSV-GENOME-POSITION\tSV-COVERAGE\tSV-PERCENT-READS\tCNA-COPY-NUMBER\tCNA-EXONS\tCNA-RATIO\tCNA-TYPE\tREARR-GENE1\tREARR-GENE2\tREARR-DESCRIPTION\tREARR-IN-FRAME?\tREARR-POS1\tREARR-POS2\tREARR-NUMBER-OF-READS\n')
        for i, l in enumerate(input_f):
            l = l.replace('\n', '')
            if not l:
                continue
            if i == 0:
                headers = l.split('\t')
                pass_col = headers.index('PASS')
                sample_col = headers.index('Sample')
                chr_col = headers.index('Chr')
                pos_col = headers.index('Start')
                ref_col = headers.index('Ref')
                alt_col = headers.index('Alt')
                class_col = headers.index('Var_Class')
                type_col = headers.index('Type')
                effect_col = headers.index('Effect')
                func_col = headers.index('Functional_Class')
                gene_code_col = headers.index('Gene_Coding')
                allele_freq_col = headers.index('AlleleFreq')
                gene_col = headers.index('Gene')
                depth_col = headers.index('Depth')
                vd_col = headers.index('VD')
                aa_chg_col = headers.index('Amino_Acid_Change')
                cosmaachg_col = headers.index('COSMIC_AA_Change')
                cosmcnt_col = headers.index('COSMIC_Cnt') if 'COSMIC_Cnt' in headers else None
                msicol = headers.index('MSI')
                cdna_chg_col = headers.index('cDNA_Change')
                gene_coding_col = headers.index('Gene_Coding')
                transcript_col = headers.index('Transcript')
                exon_col = headers.index('Exon')
                try:
                    reason_col = headers.index('Reason')
                except ValueError:
                    reason_col = None
                else:
                    status_col = reason_col - 1
                try:
                    lof_col = headers.index('LOF')
                except ValueError:
                    lof_col = None
                if not status_col and not reason_col:
                    l += '\tSignificance\tReason'
                output_f.write(l + '\n')
                if all_transcripts_output_f:
                    all_transcripts_output_f.write(l + '\n')
                continue
            fields = l.split('\t')
            if len(fields) < len(headers):
                critical('Error: len of line ' + str(i) + ' is ' + str(len(fields)) + ', which is less than the len of header (' + str(len(headers)) + ')')

            no_transcript = True
            is_canonical = False
            if fields[transcript_col] and fields[gene_coding_col] == 'transcript':
                no_transcript = False
                if fields[transcript_col].split('.')[0] in self.canonical_transcripts:
                    is_canonical = True
            if not all_transcripts_output_f:
                if not is_canonical and not no_transcript:
                    self.apply_reject_counter('not canonical transcript', True, no_transcript)
                    continue

            if fields[pass_col] != 'TRUE':
                self.apply_reject_counter('PASS=False', is_canonical, no_transcript)
                continue

            if reason_col:
                fields = fields[:-1]
            if status_col:
                fields = fields[:-1]

            sample, chrom, pos, ref, alt, aa_chg, cosm_aa_chg, gene, depth = \
                fields[sample_col], fields[chr_col], fields[pos_col], fields[ref_col], \
                fields[alt_col], fields[aa_chg_col], fields[cosmaachg_col], fields[gene_col], \
                float(fields[depth_col])

            # gene_aachg = '-'.join([gene, aa_chg])
            if 'chr' not in chrom: chrom = 'chr' + chrom
            key = '-'.join([chrom, pos, ref, alt])
            allele_freq = float(fields[allele_freq_col])

            var_class, var_type, fclass, gene_coding, effect, cdna_chg, transcript = \
                fields[class_col], fields[type_col], fields[func_col], fields[gene_code_col], \
                fields[effect_col], fields[cdna_chg_col], fields[transcript_col]
            var_type = var_type.upper()

            # debug(chrom + ':' + pos + ' ' + ref + '>' + alt + ' ' + fields[func_col] +
            #       ' ' + aa_chg + ' ' + fields[headers.index('MSI')] + ' ' + str(allele_freq))

            if var_type.startswith('PROTEIN_PROTEIN_CONTACT'):
                self.apply_reject_counter('PROTEIN_PROTEIN_CONTACT', is_canonical, no_transcript)
                continue

            region = ''
            if fields[exon_col]:
                region = fields[exon_col].split('/')[0]
                if 'intron' in var_type:
                    region = 'intron' + region

            is_lof = fields[lof_col]

            self.status = 'unknown'
            self.reason_by_status = {k: set() for k in Filtration.statuses}

            #################################
            # Checking actionable mutations #
            #################################
            aa_chg, is_act = self.check_actionable(chrom, pos, ref, alt, gene, aa_chg)
            fields[aa_chg_col] = aa_chg
            is_act = self.check_rob_hedley_actionable(gene, aa_chg, effect, region, transcript) or is_act

            if not is_act:
                if gene in self.tier_by_type_by_region_by_gene:
                    self.check_by_mut_type(cdna_chg, region, self.tier_by_type_by_region_by_gene[gene], gene)
                elif gene in self.genes_with_generic_rules:  # it can be only either in tier_by_type_by_region_by_gene or genes_with_generic_rules
                    self.check_by_general_rules(var_type, aa_chg, is_lof, gene)

            if is_lof:
                if gene in self.oncogenes:
                    if is_act:
                        warn(aa_chg + ' in ' + gene + ': LOF in a oncogoene, and actionable')
                    for s in Filtration.statuses:
                        self.reason_by_status[s].add('oncogene_lof')
                    # self.update_status('unlikely', reasons + ['oncogene_lof'], force=True)
                    # info('gene ' + gene + ' is an oncogene, and mutation is LOF. Updating status from ' + status + ' to unlikely')
                elif gene in self.suppressors:
                    # if not is_act:
                    #     warn((aa_chg if aa_chg else effect) + ' in ' + gene + ': LOF in a suppressor, but not actionable')
                    # is_act = True
                    self.update_status('likely', 'suppressor_lof')
                    # for s in Filtration.statuses:
                    #     self.reason_by_status[s].add('suppressor_lof')
                    # info('gene ' + gene + ' is a suppressor, and mutation is LOF. Making actionable. Status was ' + self.status)
                    # self.update_status(status, reasons, 'known', 'lof_in_suppressor')

            if not is_act:
                if key in self.filter_snp:
                    self.apply_reject_counter('not act and in filter_common_snp', is_canonical, no_transcript)
                    continue
                if '-'.join([gene, aa_chg]) in self.snpeff_snp:
                    self.apply_reject_counter('not act and in snpeff_snp', is_canonical, no_transcript)
                    continue
                if key in self.filter_artifacts and allele_freq < 0.2:
                    self.apply_reject_counter('not act and in filter_artifacts and AF < 0.5', is_canonical, no_transcript)
                    continue
                gmaf = fields[headers.index('GMAF')]
                if gmaf and all(not g or float(g) == 0 or float(g) > self.min_gmaf for g in gmaf.split(',')):
                    self.apply_reject_counter('not act and all GMAF > ' + str(self.min_gmaf) + ' or zero', is_canonical, no_transcript)
                    continue

            if depth < self.filt_depth:
                self.apply_reject_counter('depth < ' + str(self.filt_depth) + ' (filt_depth)', is_canonical, no_transcript)
                continue
            if fields[vd_col] < self.min_vd and allele_freq >= 0.5:
                self.apply_reject_counter('VD < ' + str(self.min_vd) + ' (min_vd) and AF >= 0.5', is_canonical, no_transcript)
                continue

            snps = re.findall(r'rs\d+', fields[3])
            if any(snp in self.snpeff_snp_rsids for snp in snps):
                self.apply_reject_counter('snp in snpeffect_export_polymorphic', is_canonical, no_transcript)
                continue

            platform = re.findall(platform_regexp, sample)[0] if platform_regexp.match(sample) else ''
            if platform.lower() not in [p.lower() for p in ['WXS', 'RNA-Seq', 'VALIDATION', 'WGS']]:
                platform = ''
            if self.platform:
                platform = self.platform

            if sample_regexp:
                if sample_regexp.match(sample):
                    sample = re.findall(self.reg_exp_sample, sample)[0]
                else:
                    self.apply_reject_counter('sample not matching ' + self.reg_exp_sample, is_canonical, no_transcript)
                    continue

            # Filter low AF MSI
            if abs(len(ref) - len(alt)) == 1:
                msi = float(fields[msicol])
                msi_fail = any([
                    msi <=  7 and allele_freq < 0.05,
                    msi ==  8 and allele_freq < 0.07,
                    msi ==  9 and allele_freq < 0.125,
                    msi == 10 and allele_freq < 0.175,
                    msi == 11 and allele_freq < 0.25,
                    msi == 12 and allele_freq < 0.3,
                    msi >  12 and allele_freq < 0.35])
                if msi_fail:
                    self.apply_reject_counter('MSI fail', is_canonical, no_transcript)
                    continue

            cosmic_counts = map(int, fields[cosmcnt_col].split()) if cosmcnt_col is not None else None
            cosm_aa_chg = self.check_by_var_class(var_class, cosm_aa_chg, cosmic_counts)
            aa_chg = self.check_by_type(var_type, aa_chg, cdna_chg, effect)

            if is_hotspot_nt(chrom, pos, ref, alt, self.hotspot_nucleotides):
                self.update_status('likely', 'hotspot_nucl_change')
            elif is_hotspot_prot(gene, aa_chg, self.hotspot_proteins):
                self.update_status('likely', 'hotspot_AA_change')

            if is_act:
                if self.min_hotspot_freq is not None and allele_freq < self.min_hotspot_freq:
                    self.apply_reject_counter('act and AF < ' + str(self.min_hotspot_freq) + ' (min_hotspot_freq)', is_canonical, no_transcript)
                    continue
            else:
                if var_type.startswith('SYNONYMOUS') or fclass.upper() == 'SILENT':
                    # Discarding any dbSNP silent mutation.
                    # Caveat: any silent mutation with entries in both dbSNP and COSMIC will be filtered,
                    # so it might filter out some somatic mutation. But keeping those will increase silent
                    # mutation nearly 10 times. Just another evidence how noisy COSMIC is.
                    if var_class == 'dbSNP' or any(f.startswith('rs') for f in fields[headers.index('ID')].split(';')):
                        self.apply_reject_counter('SYNONYMOUS and dbSNP', is_canonical, no_transcript)
                        continue
                    self.update_status('unknown', 'silent')
                if var_type.startswith('INTRON') and self.status == 'unknown':
                    self.apply_reject_counter('not act and unknown and in INTRON', is_canonical, no_transcript)
                    continue
                if 'SPLICE' in var_type and not aa_chg and self.status == 'unknown':
                    self.apply_reject_counter('not act and SPLICE and no aa_ch\g and unknown', is_canonical, no_transcript)
                    continue
                if self.min_freq and allele_freq < self.min_freq:
                    self.apply_reject_counter('not act and AF < ' + str(self.min_freq) + ' (min_freq)', is_canonical, no_transcript)
                    continue

            if self.status != 'known' and not is_act:
                if var_type.startswith('UPSTREAM'):
                    self.apply_reject_counter('not known and UPSTREAM', is_canonical, no_transcript)
                    continue
                if var_type.startswith('DOWNSTREAM'):
                    self.apply_reject_counter('not known and DOWNSTREAM', is_canonical, no_transcript)
                    continue
                if var_type.startswith('INTERGENIC'):
                    self.apply_reject_counter('not known and INTERGENIC', is_canonical, no_transcript)
                    continue
                if var_type.startswith('INTRAGENIC'):
                    self.apply_reject_counter('not known and INTRAGENIC', is_canonical, no_transcript)
                    continue
                if 'UTR_' in var_type and 'CODON' not in var_type:
                    self.apply_reject_counter('not known and not UTR_/CODON', is_canonical, no_transcript)
                    continue
                if 'NON_CODING' in gene_coding.upper():
                    self.apply_reject_counter('not known and NON_CODING', is_canonical, no_transcript)
                    continue
                if fclass.upper().startswith('NON_CODING'):
                    self.apply_reject_counter('not known and fclass=NON_CODING', is_canonical, no_transcript)
                    continue
                if var_class == 'dbSNP':
                    self.apply_reject_counter('not known and dbSNP', is_canonical, no_transcript)
                    continue

            # Ignore any variants that occur after last known critical amino acid
            aa_chg_pos_pattern = re.compile('^[A-Z](\d+).*')
            aa_pos = None
            if aa_chg_pos_pattern.match(aa_chg):
                aa_pos = int(aa_chg_pos_pattern.findall(aa_chg)[0])
                if gene in self.last_critical_aa_pos_by_gene and aa_pos >= self.last_critical_aa_pos_by_gene[gene]:
                    self.apply_reject_counter('variants occurs after last known critical amino acid', is_canonical, no_transcript)
                    continue

            if self.status != 'known' and not is_act:
                bl_gene_reasons = self.check_blacklist_genes(gene, aa_pos)
                bl_region_reasons = self.check_blacklist_regions(chrom=chrom, start=int(pos) - 1, end=int(pos) - 1 + len(ref))
                if bl_gene_reasons or bl_region_reasons:
                    if self.status == 'unknown' and 'silent' in self.reason_by_status[self.status]:
                        self.apply_reject_counter('blacklist and silent', is_canonical, no_transcript)
                        continue

                    self.apply_gene_blacklist_counter(', '.join(bl_gene_reasons + bl_region_reasons))
                    # if gene in self.gene_to_soft_filter:
                    #     self.update_status('unknown', 'blacklist gene', force=True)
                    # else:
                    if bl_gene_reasons or bl_region_reasons:
                        self.update_status('incidentalome', bl_gene_reasons + bl_region_reasons, force=True)
                    if bl_gene_reasons:
                        self.apply_reject_counter('gene blacklist', is_canonical, no_transcript)
                    elif bl_region_reasons:
                        self.apply_reject_counter('region blacklist', is_canonical, no_transcript)

                    # if float(fields[pcnt_sample_col]) > self.max_ratio:
                    # if self.freq_in_sample_by_vark:
                    #     vark = ':'.join([chrom, pos, ref, alt])
                    #     if vark in self.freq_in_sample_by_vark:
                    #         cohort_freq = self.freq_in_sample_by_vark[vark]
                    #         if cohort_freq > self.max_ratio:
                    #     self.filter_reject_counter['not known and Pcnt_sample > max_ratio (' + str(self.max_ratio) + ')'] += 1
                    #     continue

            # if gene in self.sensitizations_by_gene and (prev_gene == gene or not cur_gene_mutations):
            #     if gene_aachg in Filtration.sensitization_aa_changes:
            #         sensitization = Filtration.sensitization_aa_changes[gene_aachg]
            #         sensitizations.append(sensitization)
            #     fm_data = [sample, platform, prev_gene, pos, aa_chg_col, cdna_chg_col, chr_col, depth, allele_freq]
            #     cur_gene_mutations.append([fields, self.status, self.reason_by_status[self.status], gene_aachg, fm_data, is_canonical, no_transcript])
            # else:
            #     if cur_gene_mutations:
            #         self.print_mutations_for_one_gene(output_f, fm_output_f, all_transcripts_output_f, cur_gene_mutations, prev_gene, sensitizations)
            #         cur_gene_mutations = []
            #         sensitizations = []
            # prev_gene = gene
            #
            # if gene not in self.sensitizations_by_gene:  # sens/res mutations in spec. gene are written in print_mutations_for_one_gene
            self.print_mutation(output_f, fm_output_f, all_transcripts_output_f,
                self.status, self.reason_by_status[self.status], fields, is_canonical, no_transcript,
                fm_data=[sample, platform, gene, pos, cosm_aa_chg, aa_chg, cdna_chg, chrom, depth, allele_freq])

        info('Done.')

        counters = [['All statistics', self.all_counter, self.all_reject_counter, self.gene_blacklist_counter, self.region_blacklist_counter]]
        # if all_transcripts_output_f:
        #     counters.extend([
        #         ['No transcript', self.no_transcript_counter, self.no_transcript_reject_counter, self.no_transcript_blacklist_counter],
        #         ['Canonical', self.canonical_counter, self.canonical_reject_counter, self.canonical_blacklist_counter]])

        for title, counter, reject_counter, gene_blacklist_counter, region_blacklist_counter in counters:
            info(title + ':')
            info('    Written ' + str(counter['lines_written']) + ' lines')
            info('    Set known: ' + str(counter['known']))
            info('    Set likely: ' + str(counter['likely']))
            info('    Kept unknown: ' + str(counter['unknown']))
            info('    Incidentalome: ' + str(counter['incidentalome']))
            for reason, count in gene_blacklist_counter.items():
                info('        ' + str(count) + ' ' + reason)
            info('    Dropped: ' + str(sum(reject_counter.values())))
            for reason, count in reject_counter.items():
                info('        ' + str(count) + ' ' + reason)
            # info('    Region blacklist: ' + str(reject_counter['region blacklist']))
            # for reason, count in region_blacklist_counter.items():
            #     info('        ' + str(count) + ' ' + reason)
            info()


def parse_mut_tp53(mut_fpath):
    mut_tp53 = set()
    if verify_file(mut_fpath):
        with open(mut_fpath) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
                line = l.split('\t')
                if not line[19] or 'p.' not in line[19]:
                    continue
                prot = line[19].replace('p.', '')
                mut_tp53.add(prot)

    return mut_tp53


def is_hotspot_nt(chr, pos, ref, alt, hotspot_nucleotides):
    if len(ref) > len(alt) and alt != '-':
        ref = ref[1:]
        if len(alt) > 1:
            alt = alt[1:]
        else:
            alt = '-'
    elif len(alt) > len(ref) and ref != "-":
        alt = alt[1:]
        if len(ref) > 1:
            ref = ref[1:]
        else:
            ref = '-'
    key = '-'.join([chr, pos, ref, alt])
    return key in hotspot_nucleotides


def is_hotspot_prot(gene, aa_chg, hotspot_proteins):
    aa_chg = aa_chg.replace('p.', '')
    if not aa_chg: return False
    key = '-'.join([gene, aa_chg])
    return key in hotspot_proteins


stop_gain_pattern = re.compile('^[A-Z]+\d+\*')
fs_pattern = re.compile('^[A-Z]+(\d+)fs')
aa_snp_chg_pattern = re.compile('^[A-Z]\d+[A-Z]$')

def classify_tp53(aa_chg, pos, ref, alt, tp53_positions, tp53_groups):
    ref = ref.replace(' ', '')
    alt = alt.replace(' ', '')
    pos = pos.replace(' ', '')
    aa_chg = aa_chg.replace(' ', '')
    if pos in tp53_positions and len(ref) == 1 and len(alt) == 1:
        return 'Group 6'
    aa_chg = aa_chg.replace('p.', '')
    aa_num = 0
    if aa_chg:
        aa_num = int(re.sub('[^0-9]', '', aa_chg))
    if aa_snp_chg_pattern.match(aa_chg):
        if aa_chg in tp53_groups['Group 1']:
            return 'Group 1'
        if aa_chg in tp53_groups['Group 2']:
            return 'Group 2'
        if aa_chg in tp53_groups['Group 3']:
            return 'Group 3'
    elif stop_gain_pattern.match(aa_chg):
        if aa_num < 359:
            return 'Group 4'
    elif fs_pattern.match(aa_chg):
        if aa_num < 359:
            return 'Group 5'
    return 'NA'


def parse_specific_mutations(specific_mut_fpath):
    genes_with_generic_rules = set()
    tier_by_specific_mutations = dict()
    tier_by_type_by_region_by_gene = defaultdict(dict)
    spec_transcripts_by_aachg = defaultdict()
    dependent_mutations_by_gene = defaultdict(set)  # when other mutation is required

    with open(specific_mut_fpath) as f:
        for i, l in enumerate(f):
            if i == 0:
                continue
            l = l.replace('\n', '')
            if not l:
                continue
            line = l.split('\t')
            gene = line[0].upper()
            regions = re.findall(r'\d+', line[1])
            if '-' in line[1]:
                for region_num in range(int(regions[0]) + 1, int(regions[1])):
                    regions.append(str(region_num))
            if 'intron' in line[1]:
                regions = ['intron' + region for region in regions]
            for index in range(2, len(line) - 1):
                if line[index]:
                    mut = line[index]
                    tier = index - 1
                    if mut == 'generic':
                        genes_with_generic_rules.add(gene)
                    elif 'types' in mut:
                        types = mut.split(':')[1].split(',')
                        for region in regions:
                            tier_by_type_by_region_by_gene[gene][region] = dict()
                            for type in types:
                                tier_by_type_by_region_by_gene[gene][region][type] = tier
                    else:
                        mutations = []
                        if 'codon' in mut:
                            codons = re.findall(r'\d+', mut)
                            if '-' in mut and len(codons) == 2:
                                codons = range(int(codons[0]), int(codons[1]) + 1)
                            for region in regions:
                                for codon in codons:
                                    tier_by_specific_mutations['-'.join([gene, region, str(codon)])] = tier
                                    mutations.append('-'.join([gene, region, str(codon)]))
                        elif 'sens' in mut or 'res' in mut:
                            pattern = re.compile('\((\D+)\s+\D+\)')
                            sensitization = re.findall(pattern, mut)[0]  # like TKI
                            prot_chg = mut.split()[0].strip().replace('p.', '')
                            mutations = ['-'.join([gene, prot_chg])]
                            tier_by_specific_mutations['-'.join([gene, prot_chg])] = tier
                            dependent_mutations_by_gene[gene].add((sensitization, 'sens' if 'sens' in mut else 'res'))
                        else:
                            prot_chg = line[index].replace('p.', '').strip()
                            mutations = ['-'.join([gene, prot_chg])]
                            tier_by_specific_mutations['-'.join([gene, mut])] = tier
                        if 'NM' in line[-1] and mutations:
                            for mut in mutations:
                                spec_transcripts_by_aachg[mut] = line[-1].strip()

    return tier_by_specific_mutations, genes_with_generic_rules, \
           tier_by_type_by_region_by_gene, dependent_mutations_by_gene, spec_transcripts_by_aachg


# def is_loss_of_function(reasons, is_lof=None):
#     if is_lof is not None:
#         return is_lof
    # lof_reasons = ['frame_shift', 'stop_gained', 'start_lost', 'splice_site']
    # return any(reason in lof_reasons for reason in reasons)


def parse_genes_list(fpath):
    genes = []
    if fpath and verify_file(fpath):
        genes = [line.strip() for line in open(fpath)]
    return genes


if __name__ == '__main__':
    main()
