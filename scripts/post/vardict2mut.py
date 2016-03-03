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

from source import info, verify_file
from source.config import Config, defaults
from source import logger
from source.file_utils import adjust_path, verify_dir
from source.logger import critical
from source.prepare_args_and_cnf import determine_run_cnf, check_genome_resources, \
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug
from source.prepare_args_and_cnf import determine_sys_cnf


statuses = ['', 'known', 'likely', 'unlikely', 'unknown']  # Tier 1, 2, 3, 3
sensitization_aa_changes = {'EGFR-T790M': 'TKI'}

SIMULATE_OLD_VARDICT2MUT = False


def get_args():
    info(' '.join(sys.argv))
    info()
    description = (
        'The program will filter the VarDict output after vcf2txt.pl to '
        'candidate interpretable mutations, somatic or germline.')
    parser = OptionParser(description=description)
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser)

    parser.add_option('-o', dest='output_file')
    parser.add_option('-D', '--min-depth', dest='filt_depth', type='int', help='The minimum total depth')
    parser.add_option('-V', '--min-vd', dest='min_vd', type='int', help='The minimum reads supporting variant')
    parser.add_option('-f', '--min-freq', dest='min_freq', type='float',
                      help='The minimum allele frequency for regular variants. Default: 0.05')
    parser.add_option('-F', '--min-freq-hs', dest='min_hotspot_freq', type='float',
                      help='The minimum allele frequency hotspot somatic mutations, typically lower then -f. '
                           'Default: 0.01 or half -f, whichever is less')
    parser.add_option('-R', '--max-rate', dest='max_rate', type='float',
                      help=('If a variant is present in > [double] fraction of samples, it\'s deemed not a mutation.\n' +
                            '[default: 1.0, or no filtering].\n'
                            'Use with caution.  It\'ll filter even if it\'s in COSMIC, unless if actionable.'
                            'Don\'t use it if the sample is homogeneous. Use only in heterogeneous samples.'))
    parser.add_option('-M', '--fm', dest='is_output_fm', action='store_true', help='Output in FM\'s format')

    parser.set_usage('Usage: ' + __file__ + ' vcf2txt_res_fpath [opts] > output_fpath')

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

    set_filtering_params(cnf)

    info()

    return cnf, vcf2txt_res_fpath, adjust_path(cnf.output_file)


def main():
    cnf, vcf2txt_res_fpath, cnf.out_fpath = get_args()

    info('-' * 70)
    info('Writing to ' + cnf.out_fpath)

    do_filtering(cnf, vcf2txt_res_fpath, cnf.out_fpath)

    info('Saved to ' + cnf.out_fpath)


def do_filtering(cnf, vcf2txt_res_fpath, out_fpath):
    suppressors = parse_genes_list(adjust_path(cnf.suppressors))
    oncogenes = parse_genes_list(adjust_path(cnf.oncogenes))

    tp53_positions = []
    tp53_groups = dict()
    if cnf.ruledir:
        cnf.ruledir = verify_dir(cnf.ruledir, is_critical=True)
        tp53_groups = {'Group 1': parse_mut_tp53(join(cnf.ruledir, 'Rules', 'DNE.txt')),
                       'Group 2': parse_mut_tp53(join(cnf.ruledir, 'Rules', 'TA0-25.txt')),
                       'Group 3': parse_mut_tp53(join(cnf.ruledir, 'Rules', 'TA25-50_SOM_10x.txt'))}
        if cnf.genome.name.startswith('hg38'):
            tp53_positions = (
                '7670716 7670717 7673533 7673534 7673609 7673610 7673699 7673700 7673838 7673839 7674179 '
                '7674180 7674291 7674292 7674857 7674858 7674972 7674973 7675051 7675052 7675237 7675238 '
                '7675992 7675993 7676273 7676274 7676380 7676381 7676404 7676405 7676519 7676520 7670715 '
                '7673608 7673837 7674290 7674971 7675236 7676272 7676403 7673535 7673701 7674181 7674859 '
                '7675053 7675994 7676382 7676521').split()
        elif cnf.genome.name.startswith('hg19') or cnf.genome.name.startswith('GRCh37'):
            tp53_positions = (
                '7574034 7574035 7576851 7576852 7576927 7576928 7577017 7577018 7577156 7577157 7577497 '
                '7577498 7577609 7577610 7578175 7578176 7578290 7578291 7578369 7578370 7578555 7578556 '
                '7579310 7579311 7579591 7579592 7579698 7579699 7579722 7579723 7579837 7579838 7574033 '
                '7576926 7577155 7577608 7578289 7578554 7579590 7579721 7576853 7577019 7577499 7578177 '
                '7578371 7579312 7579700 7579839').split()

    # Set up common SNP filter
    filter_snp = set()
    if cnf.genome.filter_common_snp:
        with open(adjust_path(cnf.genome.filter_common_snp)) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
                fields = l.split('\t')
                filter_snp.add('-'.join(fields[1:5]))

    snpeff_snp = set()
    snpeff_snp_ids = set()
    if cnf.snpeffect_export_polymorphic:
        with open(adjust_path(cnf.snpeffect_export_polymorphic)) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
                fields = l.split('\t')
                if len(fields) > 11 and fields[11]:
                    snpeff_snp.add('-'.join([fields[11], fields[2]]))
                elif fields[5] != '-':
                    snpeff_snp_ids.add(fields[5])

    filter_artifacts = set()
    if cnf.genome.filter_common_artifacts:
        with open(adjust_path(cnf.genome.filter_common_artifacts)) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
                fields = l.split('\t')
                filter_artifacts.add('-'.join(fields[1:5]))

    actionable_hotspots = defaultdict(set)
    if cnf.actionable_hotspot:
        with open(adjust_path(cnf.actionable_hotspot)) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
                fields = l.split('\t')
                actionable_hotspots[fields[0]].add(fields[1])

    act_somatic = set()
    act_germline = set()

    rules = {}
    inframe_del = 'inframe-del'
    inframe_ins = 'inframe-ins'
    rules[inframe_del] = {}
    rules[inframe_ins] = {}
    if cnf.genome.actionable:
        with open(adjust_path(cnf.genome.actionable)) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
                fields = l.split('\t')
                if fields[7] == 'germline':
                    key = '-'.join(fields[1:5])
                    act_germline.add(key)
                elif fields[7] == 'somatic':
                    if fields[6] == 'rule':
                        if fields[4] == '*' and len(fields[3]) == 1:
                            key = '-'.join(fields[1:4])
                            act_somatic.add(key)
                        elif fields[5] == inframe_del:
                            rules[inframe_del].setdefault(fields[0], []).append([fields[1]] + [int (f) for f in fields[2:5]])
                        elif fields[5] == inframe_ins:
                            rules[inframe_ins].setdefault(fields[0], []).append([fields[1]] + [int (f) for f in fields[2:5]])
                    else:
                        key = '-'.join(fields[1:5])
                        act_somatic.add(key)
    hotspot_nucleotides = set()
    hotspot_proteins = set()
    if cnf.genome.compendia_ms7_hotspot:
        with open(adjust_path(cnf.genome.compendia_ms7_hotspot)) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
                fields = l.split('\t')
                hotspot_nucleotides.add('-'.join(fields[1:5]))
                hotspot_proteins.add('-'.join([fields[0], fields[6]]))

    specific_mutations, genes_with_generic_rules, genes_with_spec_types, \
        genes_with_dependent_mutations, specific_transcripts = parse_specific_mutations(adjust_path(cnf.specific_mutations))

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

    platform_regexp = re.compile('-\d\d[_-]([^_\d]+?)$')  # TODO: fix pattern
    sample_regexp = re.compile(cnf.reg_exp_sample) if cnf.reg_exp_sample else None

    filter_matches_counter = defaultdict(int)

    with open(verify_file(vcf2txt_res_fpath)) as f, open(out_fpath, 'w') as out_f:
        if cnf.is_output_fm:
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
                if not cnf.is_output_fm:
                    if not status_col and not reason_col:
                        l += '\tSignificance\tReason'
                    out_f.write(l + '\n')
                continue
            fields = l.split('\t')
            if len(fields) < len(header):
                critical('Error: len of line ' + str(i) + ' is ' + str(len(fields)) + ', which is less than the len of header (' + str(len(header)) + ')')
            if fields[pass_col] != 'TRUE':
                filter_matches_counter['PASS=False'] += 1
                continue
            if reason_col:
                fields = fields[:-1]
            if status_col:
                fields = fields[:-1]

            sample, chr, pos, ref, alt, aa_chg, gene, depth = \
                fields[sample_col], fields[chr_col], fields[pos_col], fields[ref_col], \
                fields[alt_col], fields[aa_chg_col], fields[gene_col], float(fields[depth_col])

            gene_aachg = '-'.join([gene, aa_chg])
            if 'chr' not in chr:
                chr = 'chr' + chr
            key = '-'.join([chr, pos, ref, alt])
            allele_freq = float(fields[allele_freq_col])

            is_act = False
            if all([rules, act_somatic, act_germline, actionable_hotspots, tp53_positions, tp53_groups]):
                is_act = is_actionable(chr, pos, ref, alt, gene, aa_chg, rules, act_somatic,
                                       act_germline, actionable_hotspots, tp53_positions, tp53_groups)
            if not is_act:
                if key in filter_snp:
                    filter_matches_counter['not act and in filter_common_snp'] += 1
                    continue
                if '-'.join([gene, aa_chg]) in snpeff_snp:
                    filter_matches_counter['not act and in snpeff_snp'] += 1
                    continue
                if key in filter_artifacts and allele_freq < 0.5:
                    filter_matches_counter['not act and in filter_artifacts and AF < 0.5'] += 1
                    continue

            if depth < cnf.variant_filtering.filt_depth:
                filter_matches_counter['depth < ' + str(cnf.variant_filtering.filt_depth) + ' (filt_depth)'] += 1
                continue
            if fields[vd_col] < cnf.variant_filtering.min_vd:
                filter_matches_counter['VD < ' + str(cnf.variant_filtering.min_vd) + ' (min_vd)'] += 1
                continue

            snps = re.findall(r'rs\d+', fields[3])
            if any(snp in snpeff_snp_ids for snp in snps):
                filter_matches_counter['snp in snpeffect_export_polymorphic'] += 1
                continue

            platform = re.findall(platform_regexp, sample)[0] if platform_regexp.match(sample) else ''
            if sample_regexp:
                if sample_regexp.match(sample):
                    sample = re.findall(cnf.reg_exp_sample, sample)[0]
                else:
                    filter_matches_counter['sample not matching ' + cnf.reg_exp_sample] += 1
                    continue
            status = 'unknown'
            reasons = []
            var_class, var_type, fclass, gene_coding, effect, cdna_chg, transcript = \
                fields[class_col], fields[type_col], fields[func_col], fields[gene_code_col], fields[effect_col], fields[cdna_chg_col], fields[transcript_col]
            var_type = var_type.upper()

            if var_type.startswith('PROTEIN_PROTEIN_CONTACT'):
                filter_matches_counter['PROTEIN_PROTEIN_CONTACT'] += 1
                continue

            status, reasons = check_by_var_class(var_class, status, reasons, fields, header)
            status, reasons = check_by_type(var_type, status, reasons, aa_chg, effect)

            if is_hotspot_nt(chr, pos, ref, alt, hotspot_nucleotides):
                status, reasons = update_status(status, reasons, 'likely', 'hotspot_nucl_change')
            elif is_hotspot_prot(gene, aa_chg, hotspot_proteins):
                status, reasons = update_status(status, reasons, 'likely', 'hotspot_AA_change')
            if key in act_somatic:
                status, reasons = update_status(status, reasons, 'known', 'actionable_somatic')
            elif key in act_germline:
                status, reasons = update_status(status, reasons, 'known', 'actionable_germline')
            elif is_act:
                status, reasons = update_status(status, reasons, 'known', 'actionable')

            region = ''
            if fields[exon_col]:
                region = fields[exon_col].split('/')[0]
                if 'intron' in var_type:
                    region = 'intron' + region

            is_specific_mutation = False
            if not SIMULATE_OLD_VARDICT2MUT:
                status, reasons, is_specific_mutation = check_for_specific_mutation(specific_mutations, gene, aa_chg, effect,
                                                                                    region, status, reasons, transcript, specific_transcripts)

            if not SIMULATE_OLD_VARDICT2MUT and not is_specific_mutation:
                is_lof = fields[lof_col] if lof_col else None
                if not SIMULATE_OLD_VARDICT2MUT and status != 'known' and gene in genes_with_generic_rules:
                    status, reasons = check_by_general_rules(var_type, status, reasons, aa_chg, is_lof)
                if not SIMULATE_OLD_VARDICT2MUT and status != 'known' and gene in genes_with_spec_types:
                    status, reasons = check_by_mut_type(cdna_chg, status, reasons, region, genes_with_spec_types[gene])

                if is_loss_of_function(reasons, is_lof):
                    if gene in oncogenes:
                        info('gene ' + gene + ' is an oncogene, and mutation is LOF. Updating status from ' + status + ' to unlikely')
                        status, reasons = update_status(status, reasons, 'unlikely', reasons + ['oncogene'], force=True)
                    elif gene in suppressors:
                        is_act = True
                        info('gene ' + gene + ' is a suppressor, and mutation is LOF. Updating status from ' + status + ' to known')
                        status, reasons = update_status(status, reasons, 'known', 'lof_in_suppressor')

                if is_act:
                    if allele_freq < cnf.variant_filtering.min_hotspot_freq:
                        filter_matches_counter['act and AF < ' + str(cnf.variant_filtering.min_hotspot_freq) + ' (min_hotspot_freq)'] += 1
                        continue
                    if allele_freq < 0.2 and key in act_germline:
                        filter_matches_counter['act and AF < 0.2 and is act_germline'] += 1
                        continue
                else:
                    if cnf.variant_filtering.min_freq_vardict2mut and allele_freq < cnf.variant_filtering.min_freq_vardict2mut:
                        filter_matches_counter['not act and AF < ' + str(cnf.variant_filtering.min_freq_vardict2mut) + ' (min_freq_vardict2mut)'] += 1
                        continue
                    if var_type.startswith('INTRON'):
                        filter_matches_counter['not act and in INTRON'] += 1
                        continue
                    if var_type.startswith('SYNONYMOUS'):
                        filter_matches_counter['not act and SYNONYMOUS'] += 1
                        continue
                    if fclass.upper() == 'SILENT':
                        filter_matches_counter['not act and SILENT'] += 1
                        continue
                    if var_type == 'SPLICE_REGION_VARIANT' and not aa_chg:
                        filter_matches_counter['not act and SPLICE_REGION_VARIANT'] += 1
                        continue

                if status != 'known':
                    if var_type.startswith('UPSTREAM'):
                        filter_matches_counter['not known and UPSTREAM'] += 1
                        continue
                    if var_type.startswith('DOWNSTREAM'):
                        filter_matches_counter['not known and DOWNSTREAM'] += 1
                        continue
                    if var_type.startswith('INTERGENIC'):
                        filter_matches_counter['not known and INTERGENIC'] += 1
                        continue
                    if var_type.startswith('INTRAGENIC'):
                        filter_matches_counter['not known and INTRAGENIC'] += 1
                        continue
                    if 'UTR_' in var_type and 'CODON' not in var_type:
                        filter_matches_counter['not known and not UTR_/CODON'] += 1
                        continue
                    if 'NON_CODING' in gene_coding.upper():
                        filter_matches_counter['not known and NON_CODING'] += 1
                        continue
                    if fclass.upper().startswith('NON_CODING'):
                        filter_matches_counter['not known and NON_CODING'] += 1
                        continue
                    if var_class == 'dbSNP':
                        filter_matches_counter['not known and dbSNP'] += 1
                        continue
                    if float(fields[pcnt_sample_col]) > cnf.variant_filtering.max_ratio:
                        filter_matches_counter['not known and Pcnt_sample > variant_filtering (' + str(cnf.variant_filtering.max_ratio) + ')'] += 1
                        continue

            if not SIMULATE_OLD_VARDICT2MUT:
                if gene in genes_with_dependent_mutations and (prev_gene == gene or not cur_gene_mutations):
                    if gene_aachg in sensitization_aa_changes:
                        sens_mut = sensitization_aa_changes[gene_aachg]
                        sensitizations.append(sens_mut)
                    fm_data = [sample, platform, prev_gene, pos, aa_chg_col, cdna_chg_col, chr_col, depth, allele_freq]
                    cur_gene_mutations.append([fields, status, reasons, gene_aachg, fm_data])
                else:
                    if cur_gene_mutations:
                        lines_written = print_mutations_for_one_gene(out_f, lines_written,
                            cur_gene_mutations, prev_gene, genes_with_dependent_mutations, sensitizations,
                            is_output_fm=cnf.is_output_fm)
                        cur_gene_mutations = []
                        sensitizations = []
                prev_gene = gene

            if gene not in genes_with_dependent_mutations:  # sens/res mutations in spec. gene are written in print_mutations_for_one_gene
                lines_written = print_mutation(out_f, lines_written, status, reasons, fields,
                    is_output_fm=cnf.is_output_fm, fm_data=[sample, platform, gene, pos, aa_chg_col, cdna_chg_col, chr_col, depth, allele_freq])

    info()
    info('Written ' + str(lines_written) + ' lines')

    info()
    info('Filtered stats:')
    for reason, count in filter_matches_counter.items():
        info('  ' + str(count) + ' ' + reason)
    info()


aa_chg_pattern = re.compile('^([A-Z]\d+)[A-Z]$')


def is_actionable(chr, pos, ref, alt, gene, aa_chg, rules, act_som, act_germ, act_hotspots, tp53_pos, tp53_groups):
    key = '-'.join([chr, pos, ref, alt])
    if key in act_som:
        return True
    if key in act_germ:
        return 'germline'
    if len(ref) == 1 and len(ref) == len(alt):
        key = '-'.join([chr, pos, ref])
        if key in act_som:
            return True
    if aa_chg_pattern.match(aa_chg) and gene in act_hotspots:
        act_hotspots = act_hotspots[gene]
        if aa_chg in act_hotspots:
            return True
        aa_chg_trim = re.findall(aa_chg_pattern, aa_chg)[0]
        if aa_chg_trim in act_hotspots:
            return True
    if gene == 'TP53':
        tp53_group = classify_tp53(aa_chg, pos, ref, alt, tp53_pos, tp53_groups)
        if tp53_group != 'NA':
            return 'somatic'
    if gene in rules['inframe-del'] and len(ref) > len(alt) and (len(ref) - len(alt)) % 3 == 0:
        for r in rules['inframe-del'][gene]:
            if r[0] == chr and r[1] <= int(pos) <= r[2] and len(ref) - len(alt) >= r[3]:
                return 'somatic'
    elif gene in rules['inframe-ins'] and len(ref) < len(alt) and (len(alt) - len(ref)) % 3 == 0:
        for r in rules['inframe-ins'][gene]:
            if r[0] == chr and r[1] <= int(pos) <= r[2] and len(alt) - len(ref) >= r[3]:
                return 'somatic'
    return False


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


def update_status(cur_status, reasons, new_status, new_reason, force=False):
    if not force and statuses.index(new_status) > statuses.index(cur_status):
        return cur_status, reasons
    if cur_status != new_status:
        if isinstance(new_reason, list):
            reasons = new_reason
        else:
            reasons = [new_reason]
    else:
        if isinstance(new_reason, list):
            reasons.extend(new_reason)
        else:
            reasons.append(new_reason)

    return new_status, reasons


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


def is_hotspot_prot(gene, pchg, hotspot_proteins):
    pchg = pchg.replace('p.', '')
    key = '-'.join([gene, pchg])
    return key in hotspot_proteins


def check_by_var_class(var_class, status, reasons, line, header):
    if var_class == 'ClnSNP':
        status, reasons = update_status(status, reasons, 'likely', 'clin_SNP')
    elif var_class == 'dbSNP_del':
        status, reasons = update_status(status, reasons, 'likely', 'dbSNP_del')
    elif var_class == 'ClnSNP_known':
        status, reasons = update_status(status, reasons, 'known', 'clin_SNP_known')
    elif var_class == 'ClnSNP_unknown':
        status, reasons = update_status(status, reasons, 'unknown', 'clin_SNP_unknown')
    elif var_class == 'COSMIC':
        if 'COSMIC_Cnt' in header:
            if sum([int(l) for l in line[header.index('COSMIC_Cnt')].split(',') if l]) >= 5:
                status, reasons = update_status(status, reasons, 'likely', 'COSMIC_5+')
        else:
            status, reasons = update_status(status, reasons, 'likely', 'COSMIC')
    return status, reasons


stop_gain_pattern = re.compile('^[A-Z]+\d+\*')


def check_by_type(var_type, status, reasons, aa_chg, effect):
    if 'FRAME_SHIFT' in var_type or 'FRAMESHIFT' in var_type:
        status, reasons = update_status(status, reasons, 'likely', 'frame_shift')
    elif stop_gain_pattern.match(aa_chg) or 'STOP_GAIN' in var_type:
        status, reasons = update_status(status, reasons, 'likely', 'stop_gained')
    elif 'START_LOST' in var_type and effect == 'HIGH':
        status, reasons = update_status(status, reasons, 'likely', 'start_lost')
    elif not aa_chg and 'SPLICE' in var_type and 'REGION_VARIANT' not in var_type:
        status, reasons = update_status(status, reasons, 'likely', 'splice_site')
    elif 'SPLICE_DONOR' in var_type or 'SPLICE_ACCEPTOR' in var_type:
        status, reasons = update_status(status, reasons, 'likely', 'splice_site')
    return status, reasons


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
    specific_mutations = {}
    genes_with_generic_rules = set()
    genes_with_spec_types = defaultdict(dict)
    specific_transcripts = defaultdict()

    genes_with_dependent_mutations = defaultdict(set) # other mutation is required
    with open(specific_mut_fpath) as f:
        for l in f:
            l = l.strip()
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
                            genes_with_spec_types[gene][region] = {}
                            for type in types:
                                genes_with_spec_types[gene][region][type] = tier
                    else:
                        mutations = []
                        if 'codon' in mut:
                            codons = re.findall(r'\d+', mut)
                            if '-' in mut and len(codons) == 2:
                                codons = range(int(codons[0]), int(codons[1]) + 1)
                            for region in regions:
                                for codon in codons:
                                    specific_mutations['-'.join([gene, region, str(codon)])] = tier
                                    mutations.append('-'.join([gene, region, str(codon)]))
                        elif 'sens' in mut or 'res' in mut:
                            pattern = re.compile('\((\D+)\s+\D+\)')
                            dependent_mutation = re.findall(pattern, mut)[0]
                            prot_chg = mut.split()[0].strip().replace('p.', '')
                            mutations = ['-'.join([gene, prot_chg])]
                            specific_mutations['-'.join([gene, prot_chg])] = tier
                            genes_with_dependent_mutations[gene].add(dependent_mutation)
                        else:
                            prot_chg = line[index].replace('p.', '').strip()
                            mutations = ['-'.join([gene, prot_chg])]
                            specific_mutations['-'.join([gene, mut])] = tier
                        if 'NM' in line[-1] and mutations:
                            for mut in mutations:
                                specific_transcripts[mut] = line[-1].strip()

    return specific_mutations, genes_with_generic_rules,  genes_with_spec_types, genes_with_dependent_mutations, specific_transcripts


def check_for_specific_mutation(specific_mutations, gene, aa_chg, effect, region, status, reasons, transcript, specific_transcripts):
    if aa_chg:
        gene_aachg = '-'.join([gene, aa_chg])
        if transcript and gene_aachg in specific_transcripts:
            if specific_transcripts[gene_aachg] != transcript:
                return status, reasons, False
        if gene_aachg in specific_mutations:
            tier = specific_mutations[gene_aachg]
            status, reasons = update_status(status, reasons, statuses[tier], 'manually_curated')
            return status, reasons, True
    if region and effect in ['HIGH', 'MODERATE']:
        codon = re.sub('[^0-9]', '', aa_chg)
        gene_codon_chg = '-'.join([gene, region, codon])
        if gene_codon_chg in specific_mutations:
            tier = specific_mutations[gene_codon_chg]
            status, reasons = update_status(status, reasons, statuses[tier], 'manually_curated')
            return status, reasons, True

    return status, reasons, False


def check_by_general_rules(var_type, status, reasons, aa_chg, is_lof):
    if 'splice_site' in reasons:
        status, reasons = update_status(status, reasons, 'known', 'actionable')
    elif is_loss_of_function(reasons, is_lof):
        status, reasons = update_status(status, reasons, 'known', 'actionable')
    elif 'EXON_LOSS' in var_type or 'EXON_DELETED' in var_type:
        status, reasons = update_status(status, reasons, 'known', 'actionable')
    elif status != 'unlikely':
        status, reasons = update_status(status, reasons, 'unlikely', 'not_alter_protein_function', force=True)
    return status, reasons


def is_loss_of_function(reasons, is_lof=None):
    if is_lof is not None:
        return is_lof
    lof_reasons = ['frame_shift', 'stop_gained', 'start_lost', 'splice_site']
    return any(reason in lof_reasons for reason in reasons)


def check_by_mut_type(cdna_chg, status, reasons, region, types_by_region):
    if region in types_by_region:
        for type_ in types_by_region[region]:
            if type_ in cdna_chg:
                tier = types_by_region[type_]
                status, reasons = update_status(status, reasons, statuses[tier], 'manually_curated')
    return status, reasons


def parse_genes_list(fpath):
    genes = []
    if fpath and verify_file(fpath):
        genes = [line.strip() for line in open(fpath)]
    return genes


def print_mutation(out_f, lines_written, status, reasons, fields, fm_data=None, is_output_fm=False):
    lines_written += 1
    if fm_data and is_output_fm:
        sample, platform, gene, pos, aa_chg_col, cdna_chg_col, chr_col, depth, allele_freq = fm_data
        out_f.write('\t'.join([sample, platform, 'short-variant', gene, status, fields[aa_chg_col], fields[cdna_chg_col],
                               fields[chr_col] + ':' + pos, str(depth), str(allele_freq * 100),
                               '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-'])  + '\n')
    else:
        out_f.write('\t'.join(fields + [status]) + ('\t' + ', '.join(reasons) + '\n'))
        # if status != fields[-2] or ','.join(reasons) != fields[-1]:
        #     out_f.write('\t'.join(fields[1:15] + fields[-3:] + [status]) + ('\t' + ','.join(reasons) + '\n'))
    return lines_written


def print_mutations_for_one_gene(out_f, lines_written, cur_gene_mutations, gene, genes_with_dependent_mutations,
                                 sensitizations, is_output_fm=False):
    for line in cur_gene_mutations:
        fields, status, reasons, gene_aachg, fm_data = line
        for sens_mut in genes_with_dependent_mutations[gene]:
            aa_chg = sensitization_aa_changes.keys()[sensitization_aa_changes.values().index(sens_mut)]
            if sens_mut not in sensitizations and status == 'known':
                status, reasons = update_status(status, reasons, 'likely', reasons, force=True)
            elif sens_mut not in sensitizations and status == 'likely':
                status, reasons = update_status(status, reasons, 'unlikely', aa_chg + '_required', force=True)
        lines_written = print_mutation(out_f, lines_written, status, reasons, fields, fm_data, is_output_fm)
    return lines_written


def set_filtering_params(cnf, bcbio_structure=None):
    cnf.variant_filtering.min_freq_vardict2mut = cnf.min_freq or cnf.variant_filtering.min_freq_vardict2mut or cnf.variant_filtering.min_freq
    if cnf.variant_filtering.min_freq_vardict2mut == 'bcbio' or cnf.variant_filtering.min_freq_vardict2mut is None:
        sample_min_freq = bcbio_structure.samples[0].min_af if bcbio_structure else None
        cnf.variant_filtering.min_freq_vardict2mut = sample_min_freq or defaults.default_min_freq
    cnf.variant_filtering.min_hotspot_freq = cnf.min_hotspot_freq or cnf.variant_filtering.min_hotspot_freq
    if cnf.variant_filtering.min_hotspot_freq is None or cnf.variant_filtering.min_hotspot_freq == 'default':
        cnf.variant_filtering.min_hotspot_freq = min(0.01, cnf.variant_filtering.min_freq_vardict2mut / 2)

if __name__ == '__main__':
    main()
