#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc

from collections import defaultdict
from optparse import OptionParser
from os.path import join
from os.path import exists
import time
import sys
import re

from source import verify_file
from source.config import Config, defaults
from source import logger
from source.file_utils import adjust_path, verify_dir
from source.logger import info, critical, err, warn
from source.prepare_args_and_cnf import determine_run_cnf, check_genome_resources, \
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug
from source.prepare_args_and_cnf import determine_sys_cnf


aa_chg_pattern = re.compile('^([A-Z]\d+)[A-Z]$')


def get_args():
    info(' '.join(sys.argv))
    info()
    description = (
        'The program will filter the VarDict output after vcf2txt.pl to '
        'candidate interpretable mutations, somatic or germline.')
    parser = OptionParser(description=description)
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser, threads=1)

    parser.add_option('-o', dest='output_file')
    parser.add_option('--cohort-freqs', dest='cohort_freqs_fpath')
    parser.add_option('-D', '--min-depth', dest='filt_depth', type='int', help='The minimum total depth')
    parser.add_option('-V', '--min-vd', dest='min_vd', type='int', help='The minimum reads supporting variant')
    parser.add_option('-f', '--min-freq', dest='min_freq', type='float',
                      help='The minimum allele frequency for regular variants. Default: 0.05')
    parser.add_option('-F', '--min-freq-hs', dest='min_hotspot_freq', type='float',
                      help='The minimum allele frequency hotspot somatic mutations, typically lower then -f. '
                           'Default: 0.01 or half -f, whichever is less')
    parser.add_option('-R', '--max-rate', dest='max_ratio', type='float',
                      help=('If a variant is present in > [max_ratio] fraction of samples, it\'s deemed not a mutation.\n' +
                            '[default: 1.0, or no filtering].\n'
                            'Use with caution. '
                            'It\'ll filter even if it\'s in COSMIC, unless if actionable.'
                            'Don\'t use it if the sample is homogeneous. Use only in heterogeneous samples.'))
    parser.add_option('-N', '--keep-utr-intronic', dest='keep_utr_intronic', action='store_true',
                      help='Keep all intronic and UTR in the output, but will be set as "unknown".')
    parser.add_option('-r', '--keep-max-rate', dest='keep_max_ratio', action='store_true',
                      help='Keep only those variants satisfying -R option. The option is meant to find '
                           're-occuring variants or artifacts.')

    parser.add_option('-M', '--fm', dest='is_output_fm', action='store_true', help='Output in FM\'s format')
    parser.add_option('-p', '--platform', dest='platform',
                      help='The platform, such as WXS, WGS, RNA-Seq, VALIDATION, etc. No Default. '
                           'Used when output is in FM\'s format (-M option)')

    parser.set_usage('Usage: ' + __file__ + ' vcf2txt_res_fpath [opts] -o output_fpath')

    (opts, args) = parser.parse_args()
    if len(args) < 1:
        critical('Provide the first argument - output from vcf2txt.pl')
    logger.is_debug = opts.debug

    vcf2txt_res_fpath = verify_file(args[0])

    run_cnf = determine_run_cnf(opts)
    cnf = Config(opts.__dict__, determine_sys_cnf(opts), run_cnf)
    if not cnf.genome:
        critical('Please, specify the --genome option (e.g. --genome hg19)')

    check_genome_resources(cnf)

    if not cnf.output_file:
        critical('Please, specify the output fpath with -o')

    info()

    return cnf, vcf2txt_res_fpath, adjust_path(cnf.output_file)


def main():
    cnf, vcf2txt_res_fpath, cnf.out_fpath = get_args()

    info('-' * 70)
    info('Writing to ' + cnf.out_fpath)

    f = Filtration(cnf)
    f.do_filtering(vcf2txt_res_fpath, cnf.out_fpath)

    info()
    info('Saved to ' + cnf.out_fpath)


class Filtration:
    statuses = ['', 'known', 'likely', 'unlikely', 'unknown']  # Tier 1, 2, 3, 3
    sensitization_aa_changes = {'EGFR-T790M': 'TKI'}

    def __init__(self, cnf):
        cnf.genome.compendia_ms7_hotspot    = verify_file(cnf.genome.compendia_ms7_hotspot, 'compendia_ms7_hotspot')
        cnf.genome.actionable               = verify_file(cnf.genome.actionable, 'actionable')
        cnf.genome.filter_common_snp        = verify_file(cnf.genome.filter_common_snp, 'filter_common_snp')
        cnf.genome.filter_common_artifacts  = verify_file(cnf.genome.filter_common_artifacts, 'filter_common_artifacts')
        cnf.suppressors                     = verify_file(cnf.suppressors, 'suppressors')
        cnf.oncogenes                       = verify_file(cnf.oncogenes, 'oncogenes')
        cnf.ruledir                         = verify_dir(cnf.ruledir, 'ruledir')
        cnf.snpeffect_export_polymorphic    = verify_file(cnf.snpeffect_export_polymorphic, 'snpeffect_export_polymorphic')
        cnf.actionable_hotspot              = verify_file(cnf.actionable_hotspot, 'actionable_hotspot')
        cnf.specific_mutations              = verify_file(cnf.specific_mutations, 'specific_mutations')
        if not all([cnf.genome.compendia_ms7_hotspot,
                    cnf.genome.actionable,
                    cnf.genome.filter_common_snp,
                    cnf.genome.filter_common_artifacts,
                    cnf.suppressors,
                    cnf.oncogenes,
                    cnf.ruledir,
                    cnf.snpeffect_export_polymorphic,
                    cnf.actionable_hotspot,
                    cnf.specific_mutations]):
            critical('Critical: some of the required files are not found or empty (see above)')

        self.suppressors = parse_genes_list(adjust_path(cnf.suppressors))
        self.oncogenes = parse_genes_list(adjust_path(cnf.oncogenes))

        self.reg_exp_sample = cnf.reg_exp_sample
        self.is_output_fm = cnf.is_output_fm
        self.platform = cnf.platform

        with open(verify_file(cnf.canonical_transcripts or cnf.snpeff_transcripts)) as f:
            self.canonical_transcripts = [tr.strip() for tr in f]

        self.min_freq = cnf.min_freq or cnf.variant_filtering.min_freq_vardict2mut
        self.min_hotspot_freq = cnf.min_hotspot_freq or cnf.variant_filtering.min_hotspot_freq
        if self.min_hotspot_freq is None or self.min_hotspot_freq == 'default':
            self.min_hotspot_freq = min(0.01, self.min_freq / 2)
        self.filt_depth = cnf.variant_filtering.filt_depth
        self.max_ratio = cnf.max_ratio or cnf.variant_filtering.max_ratio_vardict2mut
        self.min_vd = cnf.variant_filtering.min_vd

        # self.freq_in_sample_by_vark = dict()
        # if cnf.cohort_freqs_fpath:
        #     cohort_freqs_fpath = verify_file(cnf.cohort_freqs_fpath, is_critical=True)
        #     with open(cohort_freqs_fpath) as f:
        #         for l in f:
        #             fs = l.replace('\n', '').split()
        #             self.freq_in_sample_by_vark[fs[0]] = float(fs[1])

        self.tp53_positions = []
        self.tp53_groups = dict()
        self.tp53_groups = {'Group 1': parse_mut_tp53(join(cnf.ruledir, 'Rules', 'DNE.txt')),
                            'Group 2': parse_mut_tp53(join(cnf.ruledir, 'Rules', 'TA0-25.txt')),
                            'Group 3': parse_mut_tp53(join(cnf.ruledir, 'Rules', 'TA25-50_SOM_10x.txt'))}
        if cnf.genome.name.startswith('hg38'):
            self.tp53_positions = (
                '7670716 7670717 7673533 7673534 7673609 7673610 7673699 7673700 7673838 7673839 7674179 '
                '7674180 7674291 7674292 7674857 7674858 7674972 7674973 7675051 7675052 7675237 7675238 '
                '7675992 7675993 7676273 7676274 7676380 7676381 7676404 7676405 7676519 7676520 7670715 '
                '7673608 7673837 7674290 7674971 7675236 7676272 7676403 7673535 7673701 7674181 7674859 '
                '7675053 7675994 7676382 7676521').split()
        elif cnf.genome.name.startswith('hg19') or cnf.genome.name.startswith('GRCh37'):
            self.tp53_positions = (
                '7574034 7574035 7576851 7576852 7576927 7576928 7577017 7577018 7577156 7577157 7577497 '
                '7577498 7577609 7577610 7578175 7578176 7578290 7578291 7578369 7578370 7578555 7578556 '
                '7579310 7579311 7579591 7579592 7579698 7579699 7579722 7579723 7579837 7579838 7574033 '
                '7576926 7577155 7577608 7578289 7578554 7579590 7579721 7576853 7577019 7577499 7578177 '
                '7578371 7579312 7579700 7579839').split()

        # Set up common SNP filter
        self.filter_snp = set()
        with open(cnf.genome.filter_common_snp) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
                fields = l.split('\t')
                self.filter_snp.add('-'.join(fields[1:5]))

        self.snpeff_snp = set()
        self.snpeff_snp_rsids = set()
        with open(cnf.snpeffect_export_polymorphic) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
                fields = l.split('\t')
                snpeff_aachg = fields[2]
                snpeff_rsid = fields[5]
                if len(fields) > 11 and fields[11]:
                    snpeff_gene = fields[11]
                    self.snpeff_snp.add('-'.join([snpeff_gene, snpeff_aachg]))
                elif snpeff_rsid != '-':
                    self.snpeff_snp_rsids.add(snpeff_rsid)

        self.filter_artifacts = set()
        with open(cnf.genome.filter_common_artifacts) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
                fields = l.split('\t')
                self.filter_artifacts.add('-'.join(fields[1:5]))

        self.actionable_hotspots = defaultdict(set)
        with open(cnf.actionable_hotspot) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
                fields = l.split('\t')
                self.actionable_hotspots[fields[0]].add(fields[1])

        self.act_somatic = dict()
        self.act_germline = set()

        self.rules = defaultdict(lambda: defaultdict(list))
        # inframe_del = 'inframe-del'
        # inframe_ins = 'inframe-ins'
        # self.rules[inframe_del] = {}
        # self.rules[inframe_ins] = {}
        with open(cnf.genome.actionable) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
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
        with open(cnf.genome.compendia_ms7_hotspot) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
                fields = l.split('\t')
                if fields[5].startswith('g.'):
                    continue
                self.hotspot_nucleotides.add('-'.join(fields[1:5]))
                if not fields[6]:
                    continue
                self.hotspot_proteins.add('-'.join([fields[0], fields[6]]))

        self.tier_by_specific_mutations, \
        self.genes_with_generic_rules, \
        self.tier_by_type_by_region_by_gene, \
        self.sensitizations_by_gene, \
        self.spec_transcripts_by_aachg = parse_specific_mutations(cnf.specific_mutations)

        if not all([self.rules, self.tp53_positions, self.tp53_groups, self.act_somatic, self.act_germline, self.actionable_hotspots]):
            if not self.rules:
                err('No rules, cannot proceed')
            if not self.tp53_groups:
                err('No tp53_groups, cannot proceed')
            if not self.tp53_positions:
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

        if aa_chg_pattern.match(aa_chg) and gene in self.actionable_hotspots:
            act_hotspots = self.actionable_hotspots[gene]
            if aa_chg in act_hotspots:
                return aa_chg, self.update_status('known', 'act_hotspot')
            aa_chg_trim = re.findall(aa_chg_pattern, aa_chg)[0]
            if aa_chg_trim in act_hotspots:
                return aa_chg, self.update_status('known', 'act_hotspot')

        if gene == 'TP53':
            tp53_group = classify_tp53(aa_chg, pos, ref, alt, self.tp53_positions, self.tp53_groups)
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
                    return None, None
            if gene_aachg in self.tier_by_specific_mutations:
                tier = self.tier_by_specific_mutations[gene_aachg]
                self.update_status(Filtration.statuses[tier], 'actionable')
                # status, reasons = self.update_status(status, reasons, statuses[tier], 'actionable')

        if region and effect in ['HIGH', 'MODERATE']:
            codon = re.sub('[^0-9]', '', aa_chg)
            gene_codon_chg = '-'.join([gene, region, codon])
            if gene_codon_chg in self.tier_by_specific_mutations:
                tier = self.tier_by_specific_mutations[gene_codon_chg]
                # status, reasons = self.update_status(status, reasons, statuses[tier], 'manually_curated_codon_' + codon + '_in_exon_' + region)
                self.update_status(Filtration.statuses[tier], 'actionable_codon_' + codon + '_in_exon_' + region)

    def check_by_general_rules(self, var_type, aa_chg, is_lof, gene):
        # if 'splice_site' in reasons:
        #     status, reasons = self.update_status(status, reasons, 'known', ['general_rules'] + reasons)
        if is_lof:
            self.update_status('known', 'lof_of_gene_' + gene)
            # status, reasons = self.update_status(status, reasons, 'known', ['lof_in_gene_' + gene])
        elif 'EXON_LOSS' in var_type or 'EXON_DELETED' in var_type:
            self.update_status('known', 'exon_loss_in_gene_' + gene)
            # status, reasons = self.update_status(status, reasons, 'known', ['exon_loss_in_gene_' + gene])
        elif self.status != 'unlikely' and self.status != 'unknown':
            info(str(gene) + ' ' + str(aa_chg) + ' is in general rules, but does not alter protein function.'
                 ' Keeping status as ' + str(self.status))
        #     status, reasons = update_status(status, reasons, 'unlikely', 'but_not_alter_protein_function', force=True)

    def check_by_mut_type(self, cdna_chg, region, types_by_region, gene):
        if region in types_by_region:
            for type_ in types_by_region[region]:
                if type_ in cdna_chg:
                    tier = types_by_region[region][type_]
                    self.update_status(Filtration.statuses[tier], type_ + '_in_gene_' + gene)

    def print_mutations_for_one_gene(self, out_f, cur_gene_mutations, gene, sensitizations):
        for cur_gene_mut_info in cur_gene_mutations:
            fields, status, reasons, gene_aachg, fm_data = cur_gene_mut_info
            for sensitization, sens_or_res in self.sensitizations_by_gene[gene]:
                sensitization_aa_chg = Filtration.sensitization_aa_changes.keys()[Filtration.sensitization_aa_changes.values().index(sensitization)]

                if sensitization in sensitizations:
                    for s in Filtration.statuses:
                        self.reason_by_status[s].add(sensitization + '_' + sens_or_res)
                    # if status == 'likely':
                    # self.update_status('known', reasons + [aa_chg + '_' + sens_or_res], force=True)
                    # elif status in ['unlikely', 'unknown']:
                    #     self.update_status('likely', [aa_chg + '_' + sens_or_res], force=True)

                # if sensitization not in sensitizations and status == 'known':
                #     self.update_status('likely', reasons + [sensitization_aa_chg + '_required'], force=True)
                # elif sensitization not in sensitizations and status == 'likely':
                #     self.update_status('unlikely', sensitization_aa_chg + '_required', force=True)
            self.print_mutation(out_f, status, reasons, fields, fm_data)

    def print_mutation(self, out_f, status, reasons, fields, fm_data=None):
        self.lines_written += 1
        if status == 'unknown': self.unknown_count += 1
        if status == 'unlikely': self.unlikely_count += 1
        if status == 'likely': self.likely_count += 1
        if status == 'known': self.known_count += 1

        if fm_data and self.is_output_fm:
            sample, platform, gene, pos, cosm_aa_chg, aa_chg, cdna_chg, chrom, depth, allele_freq = fm_data
            out_f.write('\t'.join([sample, platform, 'short-variant', gene, status, aa_chg, cdna_chg,
                                   chrom + ':' + pos, str(depth), str(allele_freq * 100),
                                   '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-'])  + '\n')
        else:
            out_f.write('\t'.join(fields + [status]) + ('\t' + ', '.join(reasons) + '\n'))
            # if status != fields[-2] or ','.join(reasons) != fields[-1]:
            #     out_f.write('\t'.join(fields[1:15] + fields[-3:] + [status]) + ('\t' + ','.join(reasons) + '\n'))

    def do_filtering(self, vcf2txt_res_fpath, out_fpath):
        self.lines_written = 0
        self.unknown_count = 0
        self.unlikely_count = 0
        self.likely_count = 0
        self.known_count = 0
        self.filter_reject_counter = defaultdict(int)

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
        transcript_col = None
        effect_col = None
        exon_col = None
        pcnt_sample_col = None
        status_col = None
        reason_col = None
        lof_col = None

        lines_written = 0
        header = ''

        sensitizations = []
        cur_gene_mutations = []
        prev_gene = ''

        # platform_regexp = re.compile('-\d\d[_-]([^_\d]+?)$')  # TODO: fix pattern
        platform_regexp = re.compile('[_-]([^_\d]+?)$')  # TODO: fix pattern
        sample_regexp = re.compile(self.reg_exp_sample) if self.reg_exp_sample else None

        with open(verify_file(vcf2txt_res_fpath)) as f, open(out_fpath, 'w') as out_f:
            if self.is_output_fm:
                out_f.write('SAMPLE ID\tANALYSIS FILE LOCATION\tVARIANT-TYPE\tGENE\tSOMATIC STATUS/FUNCTIONAL IMPACT\tSV-PROTEIN-CHANGE\tSV-CDS-CHANGE\tSV-GENOME-POSITION\tSV-COVERAGE\tSV-PERCENT-READS\tCNA-COPY-NUMBER\tCNA-EXONS\tCNA-RATIO\tCNA-TYPE\tREARR-GENE1\tREARR-GENE2\tREARR-DESCRIPTION\tREARR-IN-FRAME?\tREARR-POS1\tREARR-POS2\tREARR-NUMBER-OF-READS\n')
            for i, l in enumerate(f):
                l = l.replace('\n', '')
                if not l:
                    continue
                if i == 0:
                    header = l.split('\t')
                    pass_col = header.index('PASS')
                    sample_col = header.index('Sample')
                    chr_col = header.index('Chr')
                    pos_col = header.index('Start')
                    ref_col = header.index('Ref')
                    alt_col = header.index('Alt')
                    class_col = header.index('Var_Class')
                    type_col = header.index('Type')
                    effect_col = header.index('Effect')
                    func_col = header.index('Functional_Class')
                    gene_code_col = header.index('Gene_Coding')
                    allele_freq_col = header.index('AlleleFreq')
                    gene_col = header.index('Gene')
                    depth_col = header.index('Depth')
                    vd_col = header.index('VD')
                    aa_chg_col = header.index('Amino_Acid_Change')
                    cosmaachg_col = header.index('COSMIC_AA_Change')
                    cosmcnt_col = header.index('COSMIC_Cnt') if 'COSMIC_Cnt' in header else None
                    msicol = header.index('MSI')
                    cdna_chg_col = header.index('cDNA_Change')
                    transcript_col = header.index('Transcript')
                    exon_col = header.index('Exon')
                    pcnt_sample_col = header.index('Pcnt_sample')
                    try:
                        reason_col = header.index('Reason')
                    except ValueError:
                        reason_col = None
                    else:
                        status_col = reason_col - 1
                    try:
                        lof_col = header.index('LOF')
                    except ValueError:
                        lof_col = None
                    if not self.is_output_fm:
                        if not status_col and not reason_col:
                            l += '\tSignificance\tReason'
                        out_f.write(l + '\n')
                    continue
                fields = l.split('\t')
                if len(fields) < len(header):
                    critical('Error: len of line ' + str(i) + ' is ' + str(len(fields)) + ', which is less than the len of header (' + str(len(header)) + ')')

                if fields[pass_col] != 'TRUE':
                    self.filter_reject_counter['PASS=False'] += 1
                    continue
                if fields[transcript_col] not in self.canonical_transcripts:
                    self.filter_reject_counter['non-canonical transcript'] += 1
                    continue

                if reason_col:
                    fields = fields[:-1]
                if status_col:
                    fields = fields[:-1]

                sample, chrom, pos, ref, alt, aa_chg, gene, depth = \
                    fields[sample_col], fields[chr_col], fields[pos_col], fields[ref_col], \
                    fields[alt_col], fields[aa_chg_col], fields[gene_col], float(fields[depth_col])

                # gene_aachg = '-'.join([gene, aa_chg])
                if 'chr' not in chrom: chrom = 'chr' + chrom
                key = '-'.join([chrom, pos, ref, alt])
                allele_freq = float(fields[allele_freq_col])

                var_class, var_type, fclass, gene_coding, effect, cdna_chg, transcript = \
                    fields[class_col], fields[type_col], fields[func_col], fields[gene_code_col], \
                    fields[effect_col], fields[cdna_chg_col], fields[transcript_col]
                var_type = var_type.upper()

                region = ''
                if fields[exon_col]:
                    region = fields[exon_col].split('/')[0]
                    if 'intron' in var_type:
                        region = 'intron' + region

                is_lof = fields[lof_col]

                self.status = 'unknown'
                self.reason_by_status = {k: set() for k in Filtration.statuses}

                if pos == '46924425':
                    pass

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
                        if not is_act:
                            warn(aa_chg + ' in ' + gene + ': LOF in a suppressor, but not actionable')
                        # is_act = True
                        self.update_status('likely', 'suppressor_lof')
                        # for s in Filtration.statuses:
                        #     self.reason_by_status[s].add('suppressor_lof')
                        # info('gene ' + gene + ' is a suppressor, and mutation is LOF. Making actionable. Status was ' + self.status)
                        # self.update_status(status, reasons, 'known', 'lof_in_suppressor')

                if not is_act:
                    if key in self.filter_snp:
                        self.filter_reject_counter['not act and in filter_common_snp'] += 1
                        continue
                    if '-'.join([gene, aa_chg]) in self.snpeff_snp:
                        self.filter_reject_counter['not act and in snpeff_snp'] += 1
                        continue
                    if key in self.filter_artifacts and allele_freq < 0.2:
                        self.filter_reject_counter['not act and in filter_artifacts and AF < 0.5'] += 1
                        continue

                if depth < self.filt_depth:
                    self.filter_reject_counter['depth < ' + str(self.filt_depth) + ' (filt_depth)'] += 1
                    continue
                if fields[vd_col] < self.min_vd:
                    self.filter_reject_counter['VD < ' + str(self.min_vd) + ' (min_vd)'] += 1
                    continue

                snps = re.findall(r'rs\d+', fields[3])
                if any(snp in self.snpeff_snp_rsids for snp in snps):
                    self.filter_reject_counter['snp in snpeffect_export_polymorphic'] += 1
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
                        self.filter_reject_counter['sample not matching ' + self.reg_exp_sample] += 1
                        continue

                if var_type.startswith('PROTEIN_PROTEIN_CONTACT'):
                    self.filter_reject_counter['PROTEIN_PROTEIN_CONTACT'] += 1
                    continue

                # Filter low AF MSI
                if abs(len(ref) - len(alt)) == 1:
                    msi = int(fields[msicol])
                    msi_fail = any([
                        msi <=  7 and allele_freq < 0.05,
                        msi ==  8 and allele_freq < 0.07,
                        msi ==  9 and allele_freq < 0.125,
                        msi == 10 and allele_freq < 0.175,
                        msi == 11 and allele_freq < 0.25,
                        msi == 12 and allele_freq < 0.3,
                        msi >  12 and allele_freq < 0.35])
                    if msi_fail:
                        self.filter_reject_counter['MSI fail'] += 1
                        continue

                cosmic_counts = map(int, fields[cosmcnt_col].split()) if cosmcnt_col is not None else None
                cosm_aa_chg = fields[cosmaachg_col]
                cosm_aa_chg = self.check_by_var_class(var_class, cosm_aa_chg, cosmic_counts)
                aa_chg = self.check_by_type(var_type, aa_chg, cdna_chg, effect)

                if is_hotspot_nt(chrom, pos, ref, alt, self.hotspot_nucleotides):
                    self.update_status('likely', 'hotspot_nucl_change')
                elif is_hotspot_prot(gene, aa_chg, self.hotspot_proteins):
                    self.update_status('likely', 'hotspot_AA_change')

                if is_act:
                    if self.min_hotspot_freq is not None and allele_freq < self.min_hotspot_freq:
                        self.filter_reject_counter['act and AF < ' + str(self.min_hotspot_freq) + ' (min_hotspot_freq)'] += 1
                        continue
                    # if allele_freq < 0.2 and key in self.act_germline:
                    #     self.filter_reject_counter['act and AF < 0.2 and is act_germline'] += 1
                    #     continue
                else:
                    if self.min_freq and allele_freq < self.min_freq:
                        self.filter_reject_counter['not act and AF < ' + str(self.min_freq) + ' (min_freq)'] += 1
                        continue
                    if var_type.startswith('INTRON') and self.status == 'unknown':
                        self.filter_reject_counter['not act and unknown and in INTRON'] += 1
                        continue
                    if var_type.startswith('SYNONYMOUS'):
                        self.filter_reject_counter['not act and SYNONYMOUS'] += 1
                        continue
                    if fclass.upper() == 'SILENT':
                        self.filter_reject_counter['not act and SILENT'] += 1
                        continue
                    if 'SPLICE' in var_type and not aa_chg and self.status == 'unknown':
                        self.filter_reject_counter['not act and SPLICE and no aa_chg and unknown'] += 1
                        continue

                    if self.status != 'known':
                        if var_type.startswith('UPSTREAM'):
                            self.filter_reject_counter['not known and UPSTREAM'] += 1
                            continue
                        if var_type.startswith('DOWNSTREAM'):
                            self.filter_reject_counter['not known and DOWNSTREAM'] += 1
                            continue
                        if var_type.startswith('INTERGENIC'):
                            self.filter_reject_counter['not known and INTERGENIC'] += 1
                            continue
                        if var_type.startswith('INTRAGENIC'):
                            self.filter_reject_counter['not known and INTRAGENIC'] += 1
                            continue
                        if 'UTR_' in var_type and 'CODON' not in var_type:
                            self.filter_reject_counter['not known and not UTR_/CODON'] += 1
                            continue
                        if 'NON_CODING' in gene_coding.upper():
                            self.filter_reject_counter['not known and NON_CODING'] += 1
                            continue
                        if fclass.upper().startswith('NON_CODING'):
                            self.filter_reject_counter['not known and NON_CODING'] += 1
                            continue
                        if var_class == 'dbSNP':
                            self.filter_reject_counter['not known and dbSNP'] += 1
                            continue
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
                #     cur_gene_mutations.append([fields, self.status, self.reason_by_status[self.status], gene_aachg, fm_data])
                # else:
                #     if cur_gene_mutations:
                #         self.print_mutations_for_one_gene(out_f, cur_gene_mutations, prev_gene, sensitizations)
                #         cur_gene_mutations = []
                #         sensitizations = []
                # prev_gene = gene
                #
                # if gene not in self.sensitizations_by_gene:  # sens/res mutations in spec. gene are written in print_mutations_for_one_gene
                self.print_mutation(out_f, self.status, self.reason_by_status[self.status], fields,
                    fm_data=[sample, platform, gene, pos, cosm_aa_chg, aa_chg, cdna_chg, chrom, depth, allele_freq])
                lines_written += 1

        info()
        info('Written ' + str(lines_written) + ' lines')

        info()
        info('Filtering stats:')
        info('  Dropped: ' + str(sum(self.filter_reject_counter.values())))
        for reason, count in self.filter_reject_counter.items():
            info('      ' + str(count) + ' ' + reason)
        info('  Kept unknown: ' + str(self.unknown_count))
        info('  Set unlikely: ' + str(self.unlikely_count))
        info('  Set likely: ' + str(self.likely_count))
        info('  Set known: ' + str(self.known_count))


def parse_mut_tp53(mut_fpath):
    mut_tp53 = set()
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


def classify_tp53(aa_chg, pos, ref, alt, tp53_pos, tp53_groups):
    ref = ref.replace(' ', '')
    alt = alt.replace(' ', '')
    pos = pos.replace(' ', '')
    aa_chg = aa_chg.replace(' ', '')
    if pos in tp53_pos and len(ref) == 1 and len(alt) == 1:
        return 'Group 6'
    aa_chg = aa_chg.replace('p.', '')
    aa_num = 0
    if aa_chg:
        aa_num = int(re.sub('[^0-9]', '', aa_chg))
    if fs_pattern.match(aa_chg):
        if aa_num < 359:
            return 'Group 5'
    elif stop_gain_pattern.match(aa_chg):
        if aa_num < 359:
            return 'Group 4'
    elif aa_chg_pattern.match(aa_chg):
        if aa_chg in tp53_groups['Group 1']:
            return 'Group 1'
        if aa_chg in tp53_groups['Group 2']:
            return 'Group 2'
        if aa_chg in tp53_groups['Group 3']:
            return 'Group 3'
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
