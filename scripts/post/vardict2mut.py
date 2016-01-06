#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc

from collections import defaultdict
from optparse import OptionParser
from os.path import join
from os.path import exists
import sys
import re

from source import info, verify_file
from source.config import Config
from source import logger
from source.file_utils import adjust_path
from source.logger import critical
from source.prepare_args_and_cnf import determine_run_cnf, check_genome_resources, \
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug
from source.prepare_args_and_cnf import determine_sys_cnf

aa_chg_pattern = re.compile('^([A-Z]\d+)[A-Z]$')
stop_gain_pattern = re.compile('^[A-Z]+\d+\*')
fs_pattern = re.compile('^[A-Z]+(\d+)fs')

statuses = ['unknown', 'likely', 'known']


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
        critical()
    logger.is_debug = opts.debug

    vcf2txt_res_fpath = verify_file(args[0])

    run_cnf = determine_run_cnf(opts)
    cnf = Config(opts.__dict__, determine_sys_cnf(opts), run_cnf)
    if not cnf.genome:
        critical('Please, specify the --genome option (e.g. --genome hg19)')

    check_genome_resources(cnf)

    if not cnf.output_file:
        critical('Please, specify the output fpath with -o')

    cnf.variant_filtering.min_freq = cnf.min_freq or cnf.variant_filtering.min_freq
    cnf.variant_filtering.min_hotspot_freq = cnf.min_hotspot_freq or cnf.variant_filtering.min_hotspot_freq
    if cnf.variant_filtering.min_hotspot_freq is None or cnf.variant_filtering.min_hotspot_freq == 'default':
        cnf.variant_filtering.min_hotspot_freq = min(0.01, cnf.variant_filtering.min_freq / 2)

    info()

    return cnf, vcf2txt_res_fpath, adjust_path(cnf.output_file)


def main():
    cnf, vcf2txt_res_fpath, cnf.out_fpath = get_args()

    info('-' * 70)
    info('Writing to ' + cnf.out_fpath)

    do_filtering(cnf, vcf2txt_res_fpath, cnf.out_fpath)

    info()
    info('Done, saved to ' + cnf.out_fpath)


def do_filtering(cnf, vcf2txt_res_fpath, out_fpath):
    suppressors = parse_genes_list(cnf.suppressor or [])
    oncogenes = parse_genes_list(cnf.oncogenes or [])

    tp53_hg19_positions = []
    tp53_groups = dict()
    if cnf.genome.ruledir:
        tp53_groups = {'Group 1': parse_mut_tp53(join(cnf.genome.ruledir, 'Rules', 'DNE.txt')),
                       'Group 2': parse_mut_tp53(join(cnf.genome.ruledir, 'Rules', 'TA0-25.txt')),
                       'Group 3': parse_mut_tp53(join(cnf.genome.ruledir, 'Rules', 'TA25-50_SOM_10x.txt'))}

        tp53_hg19_positions = (
            '7574034 7574035 7576851 7576852 7576927 7576928 7577017 7577018 7577156 7577157 7577497 '
            '7577498 7577609 7577610 7578175 7578176 7578290 7578291 7578369 7578370 7578555 7578556 '
            '7579310 7579311 7579591 7579592 7579698 7579699 7579722 7579723 7579837 7579838 7574033 '
            '7576926 7577155 7577608 7578289 7578554 7579590 7579721 7576853 7577019 7577499 7578177 '
            '7578371 7579312 7579700 7579839').split()

    # Set up common SNP filter
    filter_snp = set()
    if cnf.genome.filter_common_snp:
        with open(cnf.genome.filter_common_snp) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
                line = l.split('\t')
                filter_snp.add('-'.join(line[1:5]))

    snpeff_snp = set()
    snpeff_snp_ids = set()
    if cnf.snpeffect_export_polymorphic:
        with open(cnf.snpeffect_export_polymorphic) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
                line = l.split('\t')
                if len(line) > 11 and line[11]:
                    snpeff_snp.add('-'.join([line[11], line[2]]))
                elif line[5] != '-':
                    snpeff_snp_ids.add(line[5])

    filter_artifacts = set()
    if cnf.genome.filter_common_artifacts:
        with open(cnf.genome.filter_common_artifacts) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
                line = l.split('\t')
                filter_artifacts.add('-'.join(line[1:5]))

    actionable_hotspots = defaultdict(set)
    if cnf.actionable_hotspot:
        with open(cnf.actionable_hotspot) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
                line = l.split('\t')
                actionable_hotspots[line[0]].add(line[1])

    act_somatic = set()
    act_germline = set()

    rules = {}
    inframe_del = 'inframe-del'
    inframe_ins = 'inframe-ins'
    rules[inframe_del] = {}
    rules[inframe_ins] = {}
    if cnf.genome.actionable:
        with open(cnf.genome.actionable) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
                line = l.split('\t')
                if line[7] == 'germline':
                    key = '-'.join(line[1:5])
                    act_germline.add(key)
                elif line[7] == 'somatic':
                    if line[6] == 'rule':
                        if line[4] == '*' and len(line[3]) == 1:
                            key = '-'.join(line[1:4])
                            act_somatic.add(key)
                        elif line[5] == inframe_del:
                            rules[inframe_del].setdefault(line[0], []).append(line[1:5])
                        elif line[5] == inframe_ins:
                            rules[inframe_ins].setdefault(line[0], []).append(line[1:5])
                    else:
                        key = '-'.join(line[1:5])
                        act_somatic.add(key)
    hotspot_nucleotides = set()
    hotspot_proteins = set()
    if cnf.genome.compendia_ms7_hotspot:
        with open(cnf.genome.compendia_ms7_hotspot) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
                line = l.split('\t')
                hotspot_nucleotides.add('-'.join(line[1:5]))
                hotspot_proteins.add('-'.join([line[0], line[6]]))

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
    pcnt_sample_col = None

    header = ''
    with open(vcf2txt_res_fpath) as f, open(out_fpath, 'w') as out_f:
        if cnf.is_output_fm:
            out_f.write('SAMPLE ID\tANALYSIS FILE LOCATION\tVARIANT-TYPE\tGENE\tSOMATIC STATUS/FUNCTIONAL IMPACT\tSV-PROTEIN-CHANGE\tSV-CDS-CHANGE\tSV-GENOME-POSITION\tSV-COVERAGE\tSV-PERCENT-READS\tCNA-COPY-NUMBER\tCNA-EXONS\tCNA-RATIO\tCNA-TYPE\tREARR-GENE1\tREARR-GENE2\tREARR-DESCRIPTION\tREARR-IN-FRAME?\tREARR-POS1\tREARR-POS2\tREARR-NUMBER-OF-READS\n')
        for i, l in enumerate(f):
            l = l.strip()
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
                func_col = header.index('Functional_Class')
                gene_code_col = header.index('Gene_Coding')
                allele_freq_col = header.index('AlleleFreq')
                gene_col = header.index('Gene')
                depth_col = header.index('Depth')
                vd_col = header.index('VD')
                aa_chg_col = header.index('Amino_Acid_Change')
                pcnt_sample_col = header.index('Pcnt_sample')
                if not cnf.is_output_fm:
                    out_f.write(l + '\tStatus\tReason')
                continue
            line = l.split('\t')
            if line[pass_col] != 'TRUE':
                continue
            sample, chr, pos, ref, alt, aa_chg, gene, depth = \
                line[sample_col], line[chr_col], line[pos_col], line[ref_col], \
                line[alt_col], line[aa_chg_col], line[gene_col], float(line[depth_col])
            aa_chg = line[aa_chg_col]
            gene = line[gene_col]
            if 'chr' not in chr:
                chr = 'chr' + chr
            key = '-'.join([chr, pos, ref, alt])
            allele_freq = float(line[allele_freq_col])
            is_act = False
            if all([rules, act_somatic, act_germline, actionable_hotspots, tp53_hg19_positions, tp53_groups]):
                is_act = is_actionable(chr, pos, ref, alt, gene, aa_chg, rules, act_somatic,
                                       act_germline, actionable_hotspots, tp53_hg19_positions, tp53_groups)
            if not is_act and (key in filter_snp or ('-'.join([gene, aa_chg]) in snpeff_snp) or
                                   (key in filter_artifacts and allele_freq < 0.5)):
                continue
            if depth < cnf.variant_filtering.filt_depth or line[vd_col] < cnf.variant_filtering.min_vd:
                continue
            snps = re.findall(r'rs\d+', line[3])
            if any(snp in snpeff_snp_ids for snp in snps):
                continue

            sample_pattern = re.compile('-\d\d[_-]([^_\d]+?)$')
            platform = re.findall(sample_pattern, sample)[0] if sample_pattern.match(sample) else ''
            if cnf.reg_exp_sample:
                if re.compile(cnf.reg_exp_sample).match(sample):
                    sample = re.findall(cnf.reg_exp_sample, sample)[0]
                else:
                    continue
            status = 'unknown'
            reasons = []
            var_class, var_type, fclass, gene_coding = \
                line[class_col], line[type_col], line[func_col], line[gene_code_col]
            var_type = var_type.upper()
            status, reasons = check_by_var_class(var_class, status, reasons, line, header)
            status, reasons = check_by_type(var_type, status, reasons, aa_chg)

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
            if is_act:
                if allele_freq < cnf.variant_filtering.min_hotspot_freq or (allele_freq < 0.2 and key in act_germline):
                    continue
            else:
                if allele_freq < cnf.variant_filtering.min_freq or var_type.startswith('INTRON') or var_type.startswith('SYNONYMOUS') or fclass.upper() == 'SILENT' \
                        or (var_type == 'SPLICE_REGION_VARIANT' and not aa_chg):
                        continue
            if status != 'known' and (var_type.startswith('UPSTREAM') or var_type.startswith('DOWNSTREAM') or var_type.startswith('INTERGENIC') or
                                          var_type.startswith('INTRAGENIC') or ('UTR_' in var_type and 'CODON' not in var_type) or
                                              'NON_CODING' in gene_coding.upper() or fclass.upper().startswith('NON_CODING')):
                continue
            if status != 'known' and var_class == 'dbSNP':
                continue
            if status != 'known' and float(line[pcnt_sample_col]) > cnf.variant_filtering.max_ratio:
                continue
            if cnf.is_output_fm:
                out_f.write('\t'.join([sample, platform, 'short-variant', gene, status, line[aa_chg_col], line[header.index('cDNA_Change')], 'chr:' + line[chr_col],
                                 str(depth), str(allele_freq * 100), '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']))
            else:
                out_f.write('\t'.join(line + [status]) + ('\t' + ','.join(reasons)))


def is_actionable(chr, pos, ref, alt, gene, aa_chg, rules, act_som, act_germ, act_hotspots, tp53_hg19_pos, tp53_groups):
    key = '-'.join([chr, pos, ref, alt])
    if key in act_som:
        return act_som[key]
    if key in act_germ:
        return 'germline'
    if len(ref) == 1 and len(ref) == len(alt):
        key = '-'.join([chr, pos, ref])
        if key in act_som:
            return act_som[key]
    if aa_chg_pattern.match(aa_chg) and gene in act_hotspots:
        act_hotspots = act_hotspots[gene]
        if aa_chg in act_hotspots:
            return True
        aa_chg_trim = re.findall(aa_chg_pattern, aa_chg)[0]
        if aa_chg_trim in act_hotspots:
            return True
    if gene == 'TP53':
        tp53_group = classify_tp53(aa_chg, pos, ref, alt, tp53_hg19_pos, tp53_groups)
        if tp53_group != 'NA':
            return 'somatic'
    if gene in rules['inframe-del'] and len(ref) > len(alt) and (len(ref) - len(alt)) % 3 == 0:
        for r in rules['inframe-del'][gene]:
            if r[0] == chr and r[1] <= pos <= r[2] and len(ref) - len(alt) >= r[3]:
                return 'somatic'
    elif gene in rules['inframe-ins'] and len(ref) < len(alt) and (len(alt) - len(ref)) % 3 == 0:
        for r in rules['inframe-ins'][gene]:
            if r[0] == chr and r[1] <= pos <= r[2] and len(alt) - len(ref) >= r[3]:
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


def update_status(cur_status, reasons, new_status, new_reason):
    if statuses.index(new_status) < statuses.index(cur_status):
        return cur_status, reasons
    if cur_status != new_status:
        reasons = [new_reason]
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
            if int(line[header.index('COSMIC_Cnt')]) >= 5:
                status, reasons = update_status(status, reasons, 'likely', 'COSMIC_5+')
        else:
            status, reasons = update_status(status, reasons, 'likely', 'COSMIC')
    return status, reasons


def check_by_type(var_type, status, reasons, aa_chg):
    if 'FRAME_SHIFT' in var_type or 'FRAMESHIFT' in var_type:
        status, reasons = update_status(status, reasons, 'likely', 'frame_shift')
    elif stop_gain_pattern.match(aa_chg):
        status, reasons = update_status(status, reasons, 'likely', 'stop_gained')
    elif not aa_chg and 'SPLICE' in var_type and 'REGION_VARIANT' not in var_type:
        status, reasons = update_status(status, reasons, 'likely', 'splice_site')
    elif 'SPLICE_DONOR' in var_type or 'SPLICE_ACCEPTOR' in var_type:
        status, reasons = update_status(status, reasons, 'likely', 'splice_site')
    return status, reasons


def classify_tp53(aa_chg, pos, ref, alt, tp53_hg19_pos, tp53_groups):
    ref = ref.replace(' ', '')
    alt = alt.replace(' ', '')
    pos = pos.replace(' ', '')
    aa_chg = aa_chg.replace(' ', '')
    if pos in tp53_hg19_pos and len(ref) == 1 and len(alt) == 1:
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


def parse_genes_list(fpath):
    genes = [line.strip() for line in open(fpath)]
    return genes


if __name__ == '__main__':
    main()
