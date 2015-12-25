import getopt
from os.path import join
from os.path import exists
import sys
import re
import vardict2mut_config as config


long_options = "help depth= reads= min-freq= min-freq-hs= max-rate= output= " \
               "ruledir= filter_common_snp= snpeffect_export_polymorphic= filter_common_artifacts= " \
               "actionable_hotspot= actionable= compendia_ms7_hotspot= " \
               "annotation-dir= suppressors= oncogenes= report_reason".split()
short_options = "hD:V:f:F:R:o:n:"

aa_chg_pattern = re.compile('^([A-Z]\d+)[A-Z]$')
stop_gain_pattern = re.compile('^[A-Z]+\d+\*')
fs_pattern = re.compile('^[A-Z]+(\d+)fs')

statuses = ['unknown', 'likely', 'known']


def usage():
    print >> sys.stderr, "Usage: 0 [-n reg_name] input"
    print >> sys.stderr, "The program will filter the VarDict output after vcf2txt.pl " \
                         "to candidate interpretable mutations, somatic or germline."

    print >> sys.stderr, "Options:"
    print >> sys.stderr, "-h  --help                Print this usage"
    print >> sys.stderr, "-M                        Output in FM's format"
    print >> sys.stderr, "-D  --depth <int>         The minimum total depth"
    print >> sys.stderr, "-V  --reads <int>         The minimum reads supporting variant"
    print >> sys.stderr, "-n  <regexp>              The regular expression to extract sample name."
    print >> sys.stderr, "                          Default: Use as it is. For TCGA, '(TCGA-..-....)-01' is preferred for tumor sample."
    print >> sys.stderr, "-f --min-freq <double>    The minimum allele frequency for regular variants. Default: 0.05"
    print >> sys.stderr, "-R --max-rate <double>    If a variant is present in > [double] fraction of samples, it's deemed not a mutation."
    print >> sys.stderr, "                          [default: 1.0, or no filtering]."
    print >> sys.stderr, "                          Use with caution.  It'll filter even if it's in COSMIC, unless if actionable. "
    print >> sys.stderr, "                          Don't use it if the sample is homogeneous. Use only in heterogeneous samples."

    print >> sys.stderr, "-F --min-freq-hs <double> The minimum allele frequency hotspot somatic mutations, typically lower then -f."
    print >> sys.stderr, "                          Default: 0.01 or half -f, whichever is less"
    sys.exit(0)


def main(args):
    if not args:
        usage()
        sys.exit(0)

    try:
        options, vcf2txt_res_fpath = getopt.gnu_getopt(args, short_options, long_options)
    except getopt.GetoptError:
        _, exc_value, _ = sys.exc_info()
        print >> sys.stderr, exc_value
        print >> sys.stderr
        usage()
        sys.exit(2)

    vcf2txt_res_fpath = vcf2txt_res_fpath[0]
    if not exists(vcf2txt_res_fpath):
        sys.stderr('Error! Please specify VCF-file for filtering')
        sys.exit(1)

    suppressors = []
    oncogenes = []

    for opt, arg in options:
        if opt in ('-h', '--help'):
            usage()
            sys.exit(0)
        if opt == '-M':
            config.is_output_fm = True
        elif opt in ('-D', '--depth'):
            config.min_depth = float(arg)
        elif opt in ('-V', '--reads'):
            config.min_num_reads = int(arg)
        elif opt in ('-f', '--min-freq'):
            config.min_allele_freq = float(arg)
        elif opt in ('-F', '--min-freq-hs'):
            config.min_allele_freq_hotspot = float(arg)
        elif opt == '-R':
            config.max_rate = float(arg)
        elif opt == '-n':
            config.reg_exp_sample = arg
        elif opt in ('-o', '--output'):
            sys.stdout = open(arg, 'w')
        elif opt == '--report_reason':
            config.report_reason = True
        elif opt == '--ruledir':
            config.rule_dirpath = arg
        elif opt == '--filter_common_snp':
            config.filter_common_snp_fpath = arg
        elif opt == '--snpeffect_export_polymorphic':
            config.snpeffect_export_polymorphic_fpath = arg
        elif opt == '--filter_common_artifacts':
            config.filter_common_artifacts_fpath = arg
        elif opt == '--filter_common_snp':
            config.filter_common_snp_fpath = arg
        elif opt == '--actionable_hotspot':
            config.actionable_hotspot_txt_fpath = arg
        elif opt == '--actionable':
            config.actionable_txt_fpath = arg
        elif opt == '--compendia_ms7_hotspot':
            config.compendia_ms7_hotspot_fpath = arg
        elif opt == '--annotation-dir':
            config.annotation_dirpath_fpath = arg
        elif opt == '--suppressors':
            config.suppressors_fpath = arg
        elif opt == '--oncogenes':
            config.oncogenes_fpath = arg

    if not config.min_allele_freq_hotspot:
        config.min_allele_freq_hotspot = min(0.01, config.min_allele_freq / 2)

    if config.suppressors_fpath:
        suppressors = parse_genes_list(config.suppressors_fpath)
    if config.oncogenes_fpath:
        oncogenes = parse_genes_list(config.oncogenes_fpath)

    tp53_groups = {'Group 1': parse_mut_tp53(join(config.rule_dirpath, 'Rules', 'DNE.txt')),
                   'Group 2': parse_mut_tp53(join(config.rule_dirpath, 'Rules', 'TA0-25.txt')),
                   'Group 3': parse_mut_tp53(join(config.rule_dirpath, 'Rules', 'TA25-50_SOM_10x.txt'))}

    tp53_hg19_pos = (
        '7574034 7574035 7576851 7576852 7576927 7576928 7577017 7577018 7577156 7577157 7577497 '
        '7577498 7577609 7577610 7578175 7578176 7578290 7578291 7578369 7578370 7578555 7578556 '
        '7579310 7579311 7579591 7579592 7579698 7579699 7579722 7579723 7579837 7579838 7574033 '
        '7576926 7577155 7577608 7578289 7578554 7579590 7579721 7576853 7577019 7577499 7578177 '
        '7578371 7579312 7579700 7579839').split()

    # Set up common SNP filter
    filter_snp = set()
    with open(config.filter_common_snp_fpath) as f:
        for l in f:
            l = l.strip()
            if not l:
                continue
            line = l.split('\t')
            filter_snp.add('-'.join(line[1:5]))

    snpeff_snp = set()
    snpeff_snp_ids = set()
    with open(config.snpeffect_export_polymorphic_fpath) as f:
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
    with open(config.filter_common_artifacts_fpath) as f:
        for l in f:
            l = l.strip()
            if not l:
                continue
            line = l.split('\t')
            filter_artifacts.add('-'.join(line[1:5]))

    actionable_hotspots = {}
    with open(config.actionable_hotspot_txt_fpath) as f:
        for l in f:
            l = l.strip()
            if not l:
                continue
            line = l.split('\t')
            actionable_hotspots.setdefault(line[0], set()).add(line[1])

    act_somatic = set()
    act_germline = set()

    rules = {}
    inframe_del = 'inframe-del'
    inframe_ins = 'inframe-ins'
    rules[inframe_del] = {}
    rules[inframe_ins] = {}
    with open(config.actionable_txt_fpath) as f:
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
    with open(config.compendia_ms7_hotspot_fpath) as f:
        for l in f:
            l = l.strip()
            if not l:
                continue
            line = l.split('\t')
            hotspot_nucleotides.add('-'.join(line[1:5]))
            hotspot_proteins.add('-'.join([line[0], line[6]]))

    if config.is_output_fm:
        print 'SAMPLE ID\tANALYSIS FILE LOCATION\tVARIANT-TYPE\tGENE\tSOMATIC STATUS/FUNCTIONAL IMPACT\tSV-PROTEIN-CHANGE\tSV-CDS-CHANGE\tSV-GENOME-POSITION\tSV-COVERAGE\tSV-PERCENT-READS\tCNA-COPY-NUMBER\tCNA-EXONS\tCNA-RATIO\tCNA-TYPE\tREARR-GENE1\tREARR-GENE2\tREARR-DESCRIPTION\tREARR-IN-FRAME?\tREARR-POS1\tREARR-POS2\tREARR-NUMBER-OF-READS\n'
    header = ''
    with open(vcf2txt_res_fpath) as f:
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
                if not config.is_output_fm:
                    print l + '\tStatus' + ('\tReason' if config.report_reason else '')
                continue
            line = l.split('\t')
            if line[pass_col] != 'TRUE':
                continue
            sample, chr, pos, ref, alt, aa_chg, gene, depth = line[sample_col], line[chr_col], line[pos_col], line[ref_col], \
                                                                line[alt_col], line[aa_chg_col], line[gene_col], float(line[depth_col])
            aa_chg = line[aa_chg_col]
            gene = line[gene_col]
            if 'chr' not in chr:
                chr = 'chr' + chr
            key = '-'.join([chr, pos, ref, alt])
            allele_freq = float(line[allele_freq_col])
            is_act = is_actionable(chr, pos, ref, alt, gene, aa_chg, rules, act_somatic, act_germline, actionable_hotspots, tp53_hg19_pos, tp53_groups)
            if not is_act and (key in filter_snp or ('-'.join([gene, aa_chg]) in snpeff_snp) or (key in filter_artifacts and allele_freq < 0.5)):
                continue
            if depth < config.min_depth or line[vd_col] < config.min_num_reads:
                continue
            snps = re.findall(r'rs\d+', line[3])
            if any(snp in snpeff_snp_ids for snp in snps):
                continue

            sample_pattern = re.compile('-\d\d[_-]([^_\d]+?)$')
            platform = re.findall(sample_pattern, sample)[0] if sample_pattern.match(sample) else ''
            if config.reg_exp_sample:
                if re.compile(config.reg_exp_sample).match(sample):
                    sample = re.findall(config.reg_exp_sample, sample)[0]
                else:
                    continue
            status = 'unknown'
            reasons = []
            var_class, var_type, fclass, gene_coding = line[class_col], line[type_col], line[func_col], line[gene_code_col]
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
                if allele_freq < config.min_allele_freq_hotspot or (allele_freq < 0.2 and key in act_germline):
                    continue
            else:
                if allele_freq < config.min_allele_freq or var_type.startswith('INTRON') or var_type.startswith('SYNONYMOUS') or fclass.upper() == 'SILENT' \
                        or (var_type == 'SPLICE_REGION_VARIANT' and not aa_chg):
                        continue
            if status != 'known' and (var_type.startswith('UPSTREAM') or var_type.startswith('DOWNSTREAM') or var_type.startswith('INTERGENIC') or
                                          var_type.startswith('INTRAGENIC') or ('UTR_' in var_type and 'CODON' not in var_type) or
                                              'NON_CODING' in gene_coding.upper() or fclass.upper().startswith('NON_CODING')):
                continue
            if status != 'known' and var_class == 'dbSNP':
                continue
            if status != 'known' and float(line[pcnt_sample_col]) > config.max_rate:
                continue
            if config.is_output_fm:
                print '\t'.join([sample, platform, 'short-variant', gene, status, line[aa_chg_col], line[header.index('cDNA_Change')], 'chr:' + line[chr_col],
                                 str(depth), str(allele_freq * 100), '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-'])
            else:
                print '\t'.join(line + [status]) + ('\t' + ','.join(reasons) if config.report_reason else '')


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
    main(sys.argv[1:])
