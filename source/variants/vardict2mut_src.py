from collections import defaultdict, OrderedDict
from os.path import join, abspath, dirname
import sys
import re
import tabix

from source import verify_file
from source.file_utils import adjust_path, verify_dir, adjust_system_path
from source.logger import info, critical, err, warn, debug
from source.tools_from_cnf import get_system_path
from source.utils import OrderedDefaultDict


def iter_lines(fpath):
    with open(fpath) as f:
        for l in f:
            l = l.replace('\n', '')
            if not l or l.startswith('#'):
                continue
            yield l


def _read_list(reason, fpath):
    gene_d = {}
    fpath = verify_file(fpath, description=reason + ' blacklist genes file', is_critical=True)
    for l in iter_lines(fpath):
        fs = l.split('\t')
        gene_name = l.split('\t')[0]
        meta_info = l.split('\t')[1] if len(fs) == 2 else ''
        gene_d[gene_name] = meta_info
    return gene_d

def parse_gene_blacklists(cnf):
    _d = OrderedDict()
    if 'published' in cnf.variant_filtering.blacklist.genes:
        _d['freq mut gene in HGMD'] = 'published/flags_in_hgmd.txt'
        _d['freq mut gene in OMIM'] = 'published/flags_in_omim.txt'
        _d['freq mut gene'] = 'published/flags.txt'
        _d['incidentalome gene'] = 'published/incidentalome.txt'
        _d['mutSigCV gene'] = 'published/mutsigcv.txt'
    if 'low_complexity' in cnf.variant_filtering.blacklist.genes:
        _d['low complexity gene'] = 'low_complexity/low_complexity_entire_gene.txt'
    if 'repetitive_single_exome' in cnf.variant_filtering.blacklist.genes:
        _d['repetitive single exon gene'] = 'low_complexity/repetitive_single_exon_gene.txt'
    if 'abnormal_gc' in cnf.variant_filtering.blacklist.genes:
        _d['low GC gene'] = 'low_complexity/low_gc.txt'
        _d['high GC gene'] = 'low_complexity/high_gc.txt'
    if 'too_many_cosmic_mutations' in cnf.variant_filtering.blacklist.genes:
        _d['gene with too many COSMIC mutations'] = 'low_complexity/too_many_cosmic_mutations.txt'
    _d['hardfilter'] = 'blacklist.txt'

    d = OrderedDefaultDict(dict)
    for reason, fn in _d.items():
        d[reason] = _read_list(reason, join(cnf.incidentalome_dir, fn))

    return d


def load_region_blacklists(cnf):
    d = OrderedDict()
    for region_type in cnf.variant_filtering.blacklist.regions:
        fpath = verify_file(join(cnf.genome.tricky_regions, 'new', region_type + '.bed.gz'),
                            description=region_type + ' tricky regions file', is_critical=True)
        # d[reason] = build_interval_tree(fpath)
        reason = region_type.replace('_', ' ').replace('heng', 'Heng\'s').replace('lt51bp', '< 51bp')
        if 'gc' in reason:
            reason = reason.replace('to', '-').replace('gc', 'GC ') + '%'
        if 'complexity' in reason:
            reason = reason.replace('to', '-')
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


class Rule:
    def __init__(self, gene, chrom=None, start=None, end=None, length=None, ref=None,
                 required_inframe=None, indel_type=None, change=None, action=None):
        self.gene = gene
        self.chrom = chrom
        self.start = start
        self.end = end
        self.length = length
        self.ref = ref
        self.required_inframe = required_inframe
        self.indel_type = indel_type
        self.change = change
        self.action = action


class Filtration:
    statuses = ['', 'known', 'likely', 'unknown', 'uncallable regions']  # Tier 1, 2, 3, 4
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

        cnf.genome.compendia_ms7_hotspot    = verify_file(adjust_system_path(cnf.genome.compendia_ms7_hotspot), 'compendia_ms7_hotspot')
        cnf.genome.actionable               = verify_file(adjust_system_path(cnf.genome.actionable), 'actionable')
        cnf.genome.filter_common_snp        = verify_file(adjust_system_path(cnf.genome.filter_common_snp), 'filter_common_snp')
        cnf.genome.filter_common_artifacts  = verify_file(adjust_system_path(cnf.genome.filter_common_artifacts), 'filter_common_artifacts')
        cnf.genome.splice                   = verify_file(adjust_system_path(cnf.genome.splice), 'splice')
        cnf.suppressors                     = verify_file(adjust_system_path(cnf.suppressors), 'suppressors')
        cnf.oncogenes                       = verify_file(adjust_system_path(cnf.oncogenes), 'oncogenes')
        cnf.ruledir                         = verify_dir(adjust_system_path(cnf.ruledir), 'ruledir')
        cnf.snpeffect_export_polymorphic    = verify_file(adjust_system_path(cnf.snpeffect_export_polymorphic), 'snpeffect_export_polymorphic')
        cnf.actionable_hotspot              = verify_file(adjust_system_path(cnf.actionable_hotspot), 'actionable_hotspot')
        cnf.specific_mutations              = verify_file(adjust_system_path(cnf.specific_mutations), 'specific_mutations')
        cnf.last_critical_aa                = verify_file(adjust_system_path(cnf.last_critical_aa), 'last_critical_aa')
        cnf.incidentalome_dir               = verify_dir(adjust_system_path(cnf.incidentalome_dir), 'incidentalome')
        cnf.genome.tricky_regions           = verify_dir(cnf.genome.tricky_regions, 'tricky regions')
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
                    cnf.incidentalome_dir,
                    cnf.genome.tricky_regions,
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

        c = cnf.variant_filtering
        self.max_ratio = cnf.max_ratio or c.max_ratio
        self.max_sample_cnt = cnf.max_sample_cnt or c.max_sample_cnt

        self.min_freq = cnf.min_freq or c.min_freq  # for all variants
        self.act_min_freq = cnf.act_min_freq or c.act_min_freq
        if self.act_min_freq is None or self.act_min_freq == 'default':
            self.act_min_freq = min(0.01, self.min_freq / 2)
        self.germline_min_freq = c.germline_min_freq

        self.filt_depth = c.filt_depth
        self.min_vd = c.min_vd
        self.min_gmaf = c.min_gmaf

        # self.freq_in_sample_by_vark = dict()
        # if cnf.cohort_freqs_fpath:
        #     cohort_freqs_fpath = verify_file(cnf.cohort_freqs_fpath, is_critical=True)
        #     with open(cohort_freqs_fpath) as f:
        #         for l in f:
        #             fs = l.replace('\n', '').split()
        #             self.freq_in_sample_by_vark[fs[0]] = float(fs[1])

        info('Parsing filtering data...')
        self.tp53_groups = {'Group 1': parse_mut_tp53(join(cnf.ruledir, 'DNE.txt')),
                            'Group 2': parse_mut_tp53(join(cnf.ruledir, 'TA0-25.txt')),
                            'Group 3': parse_mut_tp53(join(cnf.ruledir, 'TA25-50_SOM_10x.txt'))}

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
        self.filter_rules_by_gene = defaultdict(list)
        for l in iter_lines(cnf.genome.filter_common_artifacts):
            fields = l.split('\t')
            gene, chrom, start, ref = fields[:4]
            if fields[5] == 'rule':
                action = fields[4]
                rule = Rule(gene, chrom=chrom, start=start, ref=ref, action=action)
                self.filter_rules_by_gene[gene].append(rule)
            else:
                alt = fields[4]
                self.filter_artifacts.add('-'.join([chrom, start, ref, alt]))

        self.actionable_hotspot_by_gene = defaultdict(dict)
        self.common_snps_by_gene = defaultdict(set)
        with open(cnf.actionable_hotspot) as f:
            for l in f:
                l = l.replace('\n', '')
                if not l or l.startswith('##'):
                    continue
                fields = l.split('\t')
                gene = fields[0]
                prot_change = fields[1]
                if gene.startswith('#'):  # VUS, No special treatment for now
                    gene = gene[1:]
                elif gene.startswith('^'):
                    gene = gene[1:]
                    self.common_snps_by_gene[gene].add(prot_change)
                else:
                    is_somatic = fields[2] == 'somatic'
                    self.actionable_hotspot_by_gene[gene][prot_change] = 'somatic' if is_somatic else 'germline'

        self.act_somatic = dict()
        self.act_germline = set()
        self.rules = defaultdict(list)
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
                        indel_type = ''
                        if 'indel' in fields[5]: indel_type = 'indel'
                        elif 'ins' in fields[5]: indel_type = 'ins'
                        elif 'del' in fields[5]: indel_type = 'del'
                        rule = Rule(gene=fields[0],
                                    chrom=fields[1],
                                    start=int(fields[2]),
                                    end=int(fields[3]),
                                    length=int(fields[4]),
                                    required_inframe='inframe' in fields[5],
                                    indel_type=indel_type,
                                    change=change)
                        self.rules[rule.gene].append(rule)
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
        self.gene_to_soft_filter = list(iter_lines(join(cnf.incidentalome_dir, 'soft_filter.txt')))

        info('Parsing region blacklists...')
        self.region_blacklists_by_reason = load_region_blacklists(cnf)
        for r in self.region_blacklists_by_reason.keys():
            self.region_blacklist_counter[r] = 0

        info('Parsing spreadsheat with actionable rules...')
        self.tier_by_specific_mutations, \
        self.genes_with_generic_rules, \
        self.tier_by_type_by_region_by_gene, \
        self.sensitizations_by_gene, \
        self.specific_transcripts_by_aachg = parse_specific_mutations(cnf.specific_mutations)

        if not all([self.rules, self.splice_positions_by_gene, self.act_somatic, self.act_germline, self.actionable_hotspot_by_gene]):
            if not self.rules:
                err('No rules, cannot proceed')
            if not self.splice_positions_by_gene:
                err('No tp53_positions, cannot proceed')
            if not self.act_somatic:
                err('No act_somatic, cannot proceed')
            if not self.act_germline:
                err('No act_germline, cannot proceed')
            if not self.actionable_hotspot_by_gene:
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

    platform_regexp = re.compile('[_-]([^_\d]+?)$')  # TODO: fix pattern
    def parse_platform(self, sample):
        platform = re.findall(Filtration.platform_regexp, sample)[0] \
            if Filtration.platform_regexp.match(sample) else ''
        if platform.lower() not in [p.lower() for p in ['WXS', 'RNA-Seq', 'VALIDATION', 'WGS']]:
            platform = ''
        if self.platform:
            platform = self.platform
        return platform

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

    def check_by_effect(self, var_type, aa_chg, cdna_chg, effect):
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
    def check_actionable(self, chrom, pos, ref, alt, gene, aa_chg, cosm_aa_chg, af, clnsig):
        change_len = len(alt) - len(ref)

        key = '-'.join([chrom, str(pos), ref, alt])
        if key in self.act_somatic:
            self.update_status('known', 'act_somatic')
            return 'somatic'
        if key in self.act_germline and af >= self.germline_min_freq:
            self.update_status('known', 'act_germline')
            return 'germline'

        if len(ref) == 1 and change_len == 0:  # SNP
            key = '-'.join([chrom, str(pos), ref])
            if key in self.act_somatic:
                self.update_status('known', 'act_somatic')
                return 'somatic'

        if gene in self.actionable_hotspot_by_gene and Filtration.aa_chg_trim_pattern.match(aa_chg):
            act_hotspot_by_aa_chg = self.actionable_hotspot_by_gene[gene]
            status = act_hotspot_by_aa_chg.get(aa_chg)
            if status is not None:
                if status == 'somatic' and af > self.act_min_freq:
                    self.update_status('known', 'act_hotspot_somatic')
                    return 'somatic'
                elif status == 'germline' and af > self.germline_min_freq:
                    self.update_status('known', 'act_hotspot_germline')
                    return 'germline'
            aa_chg_trim = re.findall(Filtration.aa_chg_trim_pattern, aa_chg)[0]
            status = act_hotspot_by_aa_chg.get(aa_chg_trim)
            if status is not None:
                self.update_status('known', 'act_hotspot_' + status)
                return status

        if gene == 'TP53':
            tp53_group = self.classify_tp53(aa_chg, pos, ref, alt)
            if tp53_group is not None:
                self.update_status('known', 'act_somatic_tp53_group_' + str(tp53_group))
                return 'somatic'

        if gene in self.rules:
            for r in self.rules[gene]:
                if change_len >= r.length and r.start <= pos <= r.end:
                    if r.required_inframe and change_len % 3 != 0:
                        continue
                    if any([r.indel_type == 'ins' and change_len > 0,
                            r.indel_type == 'del' and change_len < 0,
                            r.indel_type == 'indel' and change_len != 0]):
                        self.update_status('known', 'act_somatic_' + r.aa_chg)
                        return 'somatic'
        return None

    def classify_tp53(self, aa_chg, pos, ref, alt):
        aa_chg = aa_chg.replace(' ', '')
        if str(pos) in self.splice_positions_by_gene['TP53'] and len(ref) == 1 and len(alt) == 1:
            return 6
        aa_chg = aa_chg.replace('p.', '')
        aa_num = 0
        if aa_chg:
            aa_num = int(re.sub('[^0-9]', '', aa_chg))
        if aa_snp_chg_pattern.match(aa_chg):
            for i in [1, 2, 3]:
                if aa_chg in self.tp53_groups['Group ' + str(i)]:
                    return i
        elif stop_gain_pattern.match(aa_chg):
            if aa_num < 359:
                return 4
        elif fs_pattern.match(aa_chg):
            if aa_num < 359:
                return 5
        return None

    def check_rob_hedley_actionable(self, gene, aa_chg, effect, region, transcript):
        if aa_chg:
            gene_aachg = '-'.join([gene, aa_chg])
            if transcript and gene_aachg in self.specific_transcripts_by_aachg \
                and self.specific_transcripts_by_aachg[gene_aachg] != transcript:
                    return None
            if gene_aachg in self.tier_by_specific_mutations:
                tier = self.tier_by_specific_mutations[gene_aachg]
                if tier == 1:
                    self.update_status(Filtration.statuses[tier], 'actionable')
                else:
                    self.update_status(Filtration.statuses[tier], 'tier2')
                return True

        if region and effect in ['HIGH', 'MODERATE']:
            codon = re.sub('[^0-9]', '', aa_chg)
            gene_codon_chg = '-'.join([gene, region, codon])
            if gene_codon_chg in self.tier_by_specific_mutations:
                tier = self.tier_by_specific_mutations[gene_codon_chg]
                # status, reasons = self.update_status(status, reasons, statuses[tier], 'manually_curated_codon_' + codon + '_in_exon_' + region)
                if tier == 1:
                    self.update_status(Filtration.statuses[tier], 'tier2_codon_' + codon + '_in_exon_' + region)
                else:
                    self.update_status(Filtration.statuses[tier], 'actionable_codon_' + codon + '_in_exon_' + region)
                return True

    def check_by_general_rules(self, var_type, is_lof, gene):
        if gene in self.genes_with_generic_rules:
            # if 'splice_site' in reasons:
            #     status, reasons = self.update_status(status, reasons, 'known', ['general_rules'] + reasons)
            if is_lof:
                self.update_status('known', 'act_lof' + '_of_gene_' + gene)
                return True
            elif 'EXON_LOSS' in var_type or 'EXON_DELETED' in var_type:
                self.update_status('known', 'act_exon_loss' + '_in_gene_' + gene)
                return True
            else:
                return False
                # Is in general rules, but does not alter protein function

    def check_by_type_and_region(self, cdna_chg, region, gene):
        types_by_region = self.tier_by_type_by_region_by_gene.get(gene)
        if types_by_region:
            for type_ in types_by_region.get(region, []):
                if type_ in cdna_chg:
                    tier = types_by_region[region][type_]
                    self.update_status(Filtration.statuses[tier], 'act_' + type_ + '_in_gene_' + gene)
                    return True
        return False

    def fails_filters(self, chrom, pos, ref, alt, gene, aa_chg):
        for r in self.filter_rules_by_gene.get(gene, []):
            if r.action == 'ignore' and chrom == r.chrom and r.start <= pos <= r.end:
                return 'ignore'
        return None

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

    def print_mutation(self, output_f, fm_output_f, all_transcripts_output_f, status, reasons, blacklisted_reasons,
                       fields, is_canonical, no_transcript, fm_data):
        self.apply_counter('lines_written', is_canonical, no_transcript)
        self.apply_counter(status, is_canonical, no_transcript)

        if is_canonical or no_transcript:
            if fm_data and fm_output_f:
                sample, platform, gene, pos, cosm_aa_chg, aa_chg, cdna_chg, chrom, depth, allele_freq = fm_data
                fm_output_f.write('\t'.join([sample, platform, 'short-variant', gene, status, aa_chg, cdna_chg,
                                       chrom + ':' + pos, str(depth), str(allele_freq * 100),
                                       '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-'])  + '\n')
            output_f.write('\t'.join(fields + [status] + [', '.join(reasons)] + [', '.join(blacklisted_reasons)]) + '\n')
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

    def reject_mutation(self, reason, is_canonical, no_transcript, rejected_output_f, status, fields):
        self.apply_reject_counter(reason, is_canonical, no_transcript)
        if rejected_output_f:
            self.print_rejected_mutation(rejected_output_f, status, reason, fields)

    def print_rejected_mutation(self, output_f, status, reason, fields):
        output_f.write('\t'.join(fields + [status] + [reason]) + '\n')

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

    def do_filtering(self, input_f, output_f, fm_output_f=None, all_transcripts_output_f=None, rejected_output_f=None):
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
        msi_col = None
        cdna_chg_col = None
        gene_coding_col = None
        transcript_col = None
        effect_col = None
        exon_col = None
        lof_col = None
        clnsig_col = None

        prev_status_col = None
        prev_reason_col = None
        prev_incidentalome_col = None

        headers = []

        sensitizations = []
        cur_gene_mutations = []
        prev_gene = ''

        sample_regexp = re.compile(self.reg_exp_sample) if self.reg_exp_sample else None
        aa_chg_pos_regexp = re.compile('^[A-Z](\d+).*')

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
                msi_col = headers.index('MSI')
                cdna_chg_col = headers.index('cDNA_Change')
                gene_coding_col = headers.index('Gene_Coding')
                transcript_col = headers.index('Transcript')
                exon_col = headers.index('Exon')
                clnsig_col = headers.index('CLNSIG')
                lof_col = headers.index('LOF')

                if 'Significance' in headers:
                    prev_status_col = headers.index('Significance')
                    prev_reason_col = headers.index('Reason')
                    new_headers = headers + ['NewSignificance', 'NewReason', 'Incidentalome']
                    if 'Incidentalome' in headers:
                        prev_incidentalome_col = headers.index('Incidentalome')
                        new_headers = headers[:-1] + ['NewSignificance', 'NewReason', headers[-1]]
                else:
                    new_headers = headers + ['Significance', 'Reason', 'Incidentalome']

                l = '\t'.join(new_headers) + '\n'
                output_f.write(l)
                if all_transcripts_output_f:
                    all_transcripts_output_f.write(l)
                if rejected_output_f:
                    header = '\t'.join(new_headers[:-1]) + '\n'
                    rejected_output_f.write(header)
                continue

            fields = l.split('\t')
            if len(fields) < len(headers):
                critical('Error: len of line ' + str(i) + ' is ' + str(len(fields)) + ', which is less than the len of header (' + str(len(headers)) + ')')

            self.status = 'unknown'
            no_transcript = True
            is_canonical = False
            if fields[transcript_col] and fields[gene_coding_col] == 'transcript':
                no_transcript = False
                if fields[transcript_col].split('.')[0] in self.canonical_transcripts:
                    is_canonical = True
            if not all_transcripts_output_f:
                if not is_canonical and not no_transcript:
                    self.reject_mutation('not canonical transcript', True, no_transcript, rejected_output_f, self.status, fields)
                    continue

            if fields[pass_col] != 'TRUE':
                self.reject_mutation('PASS=False', is_canonical, no_transcript, rejected_output_f, self.status, fields)
                continue

            if prev_incidentalome_col:
                fields = fields[:-1]

            sample, chrom, pos, ref, alt, aa_chg, cosm_aa_chg, gene, depth = \
                fields[sample_col], fields[chr_col], int(fields[pos_col]), fields[ref_col], \
                fields[alt_col], fields[aa_chg_col], fields[cosmaachg_col], fields[gene_col], \
                float(fields[depth_col])

            if pos == 120611964:  #161514542
                pass

            # gene_aachg = '-'.join([gene, aa_chg])
            if 'chr' not in chrom: chrom = 'chr' + chrom
            nt_chg_key = '-'.join([chrom, str(pos), ref, alt])
            af = float(fields[allele_freq_col])

            var_class, var_type, fclass, gene_coding, effect, cdna_chg, transcript = \
                fields[class_col], fields[type_col], fields[func_col], fields[gene_code_col], \
                fields[effect_col], fields[cdna_chg_col], fields[transcript_col]
            var_type = var_type.upper()

            if var_type.startswith('PROTEIN_PROTEIN_CONTACT'):
                self.reject_mutation('PROTEIN_PROTEIN_CONTACT', is_canonical, no_transcript, None, self.status, fields)
                continue

            if depth < self.filt_depth:
                self.reject_mutation('depth < ' + str(self.filt_depth) + ' (filt_depth)', is_canonical, no_transcript,
                                     rejected_output_f, self.status, fields)
                continue
            if fields[vd_col] < self.min_vd and af >= 0.5:
                self.reject_mutation('VD < ' + str(self.min_vd) + ' (min_vd) and AF >= 0.5', is_canonical, no_transcript,
                                     rejected_output_f, self.status, fields)
                continue

            region = ''
            if fields[exon_col]:
                region = fields[exon_col].split('/')[0]
                if 'intron' in var_type:
                    region = 'intron' + region

            is_lof = fields[lof_col]

            self.reason_by_status = {k: set() for k in Filtration.statuses}

            #################################
            # Checking actionable mutations #
            #################################
            actionability = \
                self.check_actionable(chrom, pos, ref, alt, gene, aa_chg, cosm_aa_chg, af, fields[clnsig_col]) or \
                self.check_rob_hedley_actionable(gene, aa_chg, effect, region, transcript) or \
                self.check_by_type_and_region(cdna_chg, region, gene) or \
                self.check_by_general_rules(var_type, is_lof, gene)

            fail_reason = self.fails_filters(chrom, pos, ref, alt, gene, aa_chg)
            if fail_reason:
                self.reject_mutation(fail_reason, is_canonical, no_transcript, rejected_output_f, self.status, fields)
                continue

            if gene in self.common_snps_by_gene and aa_chg in self.common_snps_by_gene[gene]:
                self.reject_mutation('common SNP', is_canonical, no_transcript, rejected_output_f, self.status, fields)
                continue

            if is_lof:
                if gene in self.oncogenes:
                    for s in Filtration.statuses:
                        self.reason_by_status[s].add('oncogene_lof')
                elif gene in self.suppressors:
                    self.update_status('likely', 'suppressor_lof')

            if not actionability:
                if nt_chg_key in self.filter_snp:
                    self.reject_mutation('not act and in filter_common_snp', is_canonical, no_transcript,
                                         rejected_output_f, self.status, fields)
                    continue
                if nt_chg_key in self.filter_artifacts and af < 0.35:
                    self.reject_mutation('not act and in filter_artifacts and AF < 0.35', is_canonical, no_transcript,
                                         rejected_output_f, self.status, fields)
                    continue
                gmaf = fields[headers.index('GMAF')]
                if gmaf and all(not g or float(g) > self.min_gmaf for g in gmaf.split(',')):
                    self.reject_mutation('not act and all GMAF > ' + str(self.min_gmaf), is_canonical, no_transcript,
                                         rejected_output_f, self.status, fields)
                    continue
                clncheck = check_clnsig(fields[clnsig_col])
                if clncheck == 'dbSNP':  # Even if it's COSMIC in status, it's going to be filtered in case of low ClinVar significance
                    self.reject_mutation('clnsig dbSNP', is_canonical, no_transcript, rejected_output_f, self.status, fields)
                    continue
                if '-'.join([gene, aa_chg]) in self.snpeff_snp and clncheck != 'ClnSNP_known':
                    self.reject_mutation('not act and not ClnSNP_known and in snpeff_snp', is_canonical, no_transcript,
                                         rejected_output_f, self.status, fields)
                    continue

            snps = re.findall(r'rs\d+', fields[3])
            if any(snp in self.snpeff_snp_rsids for snp in snps):
                self.reject_mutation('snp in snpeffect_export_polymorphic', is_canonical, no_transcript,
                                     rejected_output_f, self.status, fields)
                continue

            if sample_regexp:
                if sample_regexp.match(sample):
                    sample = re.findall(self.reg_exp_sample, sample)[0]
                else:
                    continue

            # Filter low AF MSI
            msi = float(fields[msi_col])
            if abs(len(ref) - len(alt)) == 1 and msi > 3:
                msi_fail = any([
                    msi <=  7 and af < 0.03,
                    msi ==  8 and af < 0.06,
                    msi ==  9 and af < 0.125,
                    msi == 10 and af < 0.175,
                    msi == 11 and af < 0.25,
                    msi == 12 and af < 0.3,
                    msi >  12 and af < 0.35])
                if msi_fail:
                    self.reject_mutation('MSI fail', is_canonical, no_transcript, rejected_output_f, self.status, fields)
                    continue

            cosmic_counts = map(int, fields[cosmcnt_col].split()) if cosmcnt_col is not None else None
            cosm_aa_chg = self.check_by_var_class(var_class, cosm_aa_chg, cosmic_counts)
            aa_chg = self.check_by_effect(var_type, aa_chg, cdna_chg, effect)

            if is_hotspot_nt(chrom, pos, ref, alt, self.hotspot_nucleotides):
                self.update_status('likely', 'hotspot_nucl_change')
            elif is_hotspot_prot(gene, aa_chg, self.hotspot_proteins):
                self.update_status('likely', 'hotspot_AA_change')

            if actionability:
                if actionability == 'germline' and af < self.germline_min_freq:
                    self.reject_mutation('act germline and AF < ' + str(self.act_min_freq), is_canonical, no_transcript,
                                         rejected_output_f, self.status, fields)
                    continue
                if af < self.act_min_freq:
                    self.reject_mutation('act somatic and AF < ' + str(self.germline_min_freq), is_canonical, no_transcript,
                                         rejected_output_f, self.status, fields)
                    continue
            else:
                if var_type.startswith('SYNONYMOUS') or fclass.upper() == 'SILENT':
                    # Discarding any dbSNP silent mutation.
                    # Caveat: any silent mutation with entries in both dbSNP and COSMIC will be filtered,
                    # so it might filter out some somatic mutation. But keeping those will increase silent
                    # mutation nearly 10 times. Just another evidence how noisy COSMIC is.
                    if var_class == 'dbSNP' or any(f.startswith('rs') for f in fields[headers.index('ID')].split(';')):
                        self.reject_mutation('SYNONYMOUS and dbSNP', is_canonical, no_transcript,
                                             rejected_output_f, self.status, fields)
                        continue
                    self.update_status('unknown', 'silent')
                if var_type.startswith('INTRON') and self.status == 'unknown':
                    self.reject_mutation('not act and unknown and in INTRON', is_canonical, no_transcript,
                                         rejected_output_f, self.status, fields)
                    continue
                if 'SPLICE' in var_type and not aa_chg and self.status == 'unknown':
                    self.reject_mutation('not act and SPLICE and no aa_ch\g and unknown', is_canonical, no_transcript,
                                         rejected_output_f, self.status, fields)
                    continue
                if self.min_freq and af < self.min_freq:
                    self.reject_mutation('not act and AF < ' + str(self.min_freq) + ' (min_freq)', is_canonical, no_transcript,
                                         rejected_output_f, self.status, fields)
                    continue

            if not actionability and self.status != 'known':
                if var_type.startswith('UPSTREAM'):
                    self.reject_mutation('not known and UPSTREAM', is_canonical, no_transcript, rejected_output_f, self.status, fields)
                    continue
                if var_type.startswith('DOWNSTREAM'):
                    self.reject_mutation('not known and DOWNSTREAM', is_canonical, no_transcript, rejected_output_f, self.status, fields)
                    continue
                if var_type.startswith('INTERGENIC'):
                    self.reject_mutation('not known and INTERGENIC', is_canonical, no_transcript, rejected_output_f, self.status, fields)
                    continue
                if var_type.startswith('INTRAGENIC'):
                    self.reject_mutation('not known and INTRAGENIC', is_canonical, no_transcript, rejected_output_f, self.status, fields)
                    continue
                if 'UTR_' in var_type and 'CODON' not in var_type:
                    self.reject_mutation('not known and not UTR_/CODON', is_canonical, no_transcript, rejected_output_f, self.status, fields)
                    continue
                if 'NON_CODING' in gene_coding.upper():
                    self.reject_mutation('not known and NON_CODING', is_canonical, no_transcript, rejected_output_f, self.status, fields)
                    continue
                if fclass.upper().startswith('NON_CODING'):
                    self.reject_mutation('not known and fclass=NON_CODING', is_canonical, no_transcript, rejected_output_f, self.status, fields)
                    continue
                if var_class == 'dbSNP':
                    self.reject_mutation('not known and dbSNP', is_canonical, no_transcript, rejected_output_f, self.status, fields)
                    continue

            # Ignore any variants that occur after last known critical amino acid
            aa_pos = None
            if aa_chg_pos_regexp.match(aa_chg):
                aa_pos = int(aa_chg_pos_regexp.findall(aa_chg)[0])
            if aa_pos is not None:
                if gene in self.last_critical_aa_pos_by_gene and aa_pos >= self.last_critical_aa_pos_by_gene[gene]:
                    self.reject_mutation('variants occurs after last known critical amino acid', is_canonical, no_transcript,
                                         rejected_output_f, self.status, fields)
                    continue

            # if not actionability and self.status != 'known':
            bl_gene_reasons = self.check_blacklist_genes(gene, aa_pos)
            bl_region_reasons = self.check_blacklist_regions(chrom=chrom, start=pos - 1, end=pos - 1 + len(ref))
            blacklisted_reasons = bl_gene_reasons + bl_region_reasons
            if bl_gene_reasons or bl_region_reasons:
                # if self.status == 'unknown' and 'silent' in self.reason_by_status[self.status]:
                #     self.reject_mutation('blacklist and silent', is_canonical, no_transcript, rejected_output_f, self.status, fields)
                #     continue
                self.apply_gene_blacklist_counter(', '.join(bl_gene_reasons + bl_region_reasons))
            if 'hardfilter' in bl_gene_reasons and not actionability and self.status != 'known':
                self.reject_mutation('blacklist gene', is_canonical, no_transcript, rejected_output_f, self.status, fields)
                continue

            # if bl_region_reasons and not actionability and self.status != 'known' and gene not in self.gene_to_soft_filter:
            #     self.update_status('uncallable regions', bl_region_reasons, force=True)

            # if not is_act:
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
                self.status, self.reason_by_status[self.status], blacklisted_reasons, fields, is_canonical, no_transcript,
                fm_data=[sample, self.parse_platform(sample), gene, str(pos), cosm_aa_chg, aa_chg, cdna_chg, chrom, depth, af])

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
            info('    In uncallable regions: ' + str(counter['uncallable regions']))
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


def is_hotspot_nt(chrom, pos, ref, alt, hotspot_nucleotides):
    if len(ref) > len(alt) and alt != '-':
        ref = ref[1:]
        if len(alt) > 1:
            alt = alt[1:]
        else:
            alt = '-'
    elif len(alt) > len(ref) and ref != '-':
        alt = alt[1:]
        if len(ref) > 1:
            ref = ref[1:]
        else:
            ref = '-'
    key = '-'.join([chrom, str(pos), ref, alt])
    return key in hotspot_nucleotides


def is_hotspot_prot(gene, aa_chg, hotspot_proteins):
    aa_chg = aa_chg.replace('p.', '')
    if not aa_chg: return False
    key = '-'.join([gene, aa_chg])
    return key in hotspot_proteins


stop_gain_pattern = re.compile('^[A-Z]+\d+\*')
fs_pattern = re.compile('^[A-Z]+(\d+)fs')
aa_snp_chg_pattern = re.compile('^[A-Z]\d+[A-Z]$')


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


def check_clnsig(clnsig):
    """
    :param clnsig: a |-separate values of variant clinical significance from ClinVar:
            0 - uncertain significance,
            1 - not provided,
            2 - benign,
            3 - likely benign,
            4 - likely pathogenic,
            5 - pathogenic,
            6 - drug response,
            7 - histocompatibility,
            255 - other
    :return: variant class:
            ClnSNP_known - mostly (likely) pathogenic
            ClnSNP_unknown - mostly uncertain/other
            dbSNP - mostly (likely) benign
    """
    if not clnsig:
        return None
    flag255 = 0
    flagno = 0
    flagyes = 0
    flags = 0
    for cl in re.split('\||,', clnsig):
        cl = int(cl)
        if 4 <= cl <= 6:  # likely pathogenic, pathogenic, drug response
            flagyes += 1
        if 2 <= cl <= 3:  # benign, likely benign
            flagno += 1
        if cl == 255 or cl == 0:  # uncertain, other
            flag255 += 1
        if cl == 1:  # not provided
            flags += 1

    if flagyes:
        if flagyes > 1:
            return 'ClnSNP_known'  # twice pathogenic
        if flagyes > 0 and flagno == 0:
            return 'ClnSNP_known'  # pathogenic and not benign
        if flags and flagyes >= flagno and flagno <= 1 and flagyes / flags >= 0.5:
            return 'ClnSNP_known'
        else:
            return 'ClnSNP_unknown'

    if flagno > 1 and flagno >= flag255:  # benign is > 2 and > "not provided"
        return 'dbSNP'
    if flag255 > 0:
        return 'ClnSNP_unknown'  # Keep unknown significant variants
    return 'dbSNP'
