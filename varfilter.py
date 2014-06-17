#!/usr/bin/env python
from collections import OrderedDict
import re

import sys
from source.config import Defaults
from source.logger import err

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from source.main import read_opts_and_cnfs, check_system_resources
from source.runner import run_one
from source.utils import get_java_tool_cmdline, intermediate_fname, iterate_file, call


def main(args):
    defaults = Defaults.variant_filtering

    cnf = read_opts_and_cnfs(
        description=
        'The program will filter an annotated VCF file by SnpEff using dbSNP and COSMIC, '
        'adding PASS or REJECT into the FILTER column.\n'
        '\n'
        'A novel variant (non-dbSNP, non-COSMIC) is considered false positive '
        'if all three conditions (-r -f -n) are met. Any variant meeting the -p '
        'or -q conditions are also considered likely false positive. '
        'False positive variants are annotated REJECT in column FILTER, PASS otherwise.',

        extra_opts=[
            (['--expression'], dict(
                dest='expression',
                help='Filtering line for SnpSift. Default is ' + defaults['expression']
            )),

            (['-u'], dict(
                dest='count_undetermined',
                action='store_false',
                help='Undeteremined won\'t be counted for the sample count.'
            )),

            (['-b'], dict(
                dest='bias',
                action='store_true',
                help='Novel or dbSNP variants with strand bias "2;1" or "2;0" '
                     'and AF < 0.3 will be considered as false positive.'
            )),

            (['-r'], dict(
                dest='fraction',
                type='double',
                help='When a novel variant is present in more than [fraction] '
                     'of samples and mean allele frequency is less than [freq], '
                     'it\'s considered as likely false positive. Default %f. '
                     'Used with -f and -n' % defaults['fraction'],
            )),
            (['-f'], dict(
                dest='freq',
                type='double',
                help='When the average allele frequency is also below the [freq], '
                     'the variant is considered likely false positive. '
                     'Default %f. Used with -r and -n' % defaults['freq'],
            )),
            (['-n'], dict(
                dest='sample_cnt',
                type='int',
                help='When the variant is detected in greater or equal [sample_cnt] '
                     'samples, the variant is considered likely false positive. '
                     'Default %d. Used with -r and -f' % defaults['sample_cnt'],
            )),

            (['-R'], dict(
                dest='max_ratio',
                type='double',
                help='When a variant is present in more than [fraction] of samples, '
                     'and AF < 0.3, it\'s considered as likely false positive, '
                     'even if it\'s in COSMIC. Default %f.' % defaults['max_ratio'],
            )),

            (['-F'], dict(
                dest='min_freq',
                type='double',
                help='When individual allele frequency < feq for variants, '
                     'it was considered likely false poitives. '
                     'Default %f' % defaults['min_freq'],
            )),

            (['-p'], dict(
                dest='min_p_mean',
                type='int',
                help='The minimum mean position in reads for variants.'
                     'Default %d bp' % defaults['min_p_mean'],
            )),
            (['-q'], dict(
                dest='min_q_mean',
                type='double',
                help='The minimum mean base quality phred score for variant.'
                     'Default %d' % defaults['min_q_mean'],
            )),
            (['-P'], dict(
                dest='filt_p_mean',
                type='int',
                help='The filtering mean position in reads for variants. '
                     'The raw variant will be filtered on first place if the mean '
                     'posititon is less then [filt_p_mean]. '
                     'Default %d bp' % defaults['filt_p_mean'],
            )),
            (['-Q'], dict(
                dest='filt_q_mean',
                type='double',
                help='The filtering mean base quality phred score for variants. '
                     'The raw variant will be filtered on first place  '
                     'if the mean quality is less then [filt_q_mean]. '
                     'Default %f' % defaults['filt_q_mean'],
            )),

            (['-M'], dict(
                dest='mean_mq',
                type='double',
                help='The filtering mean mapping quality score for variants. '
                     'The raw variant will be filtered if the mean mapping quality '
                     'score is less then specified. Default %d' % defaults['mean_mq'],
            )),
            (['-D'], dict(
                dest='filt_depth',
                type='int',
                help='The filtering total depth. The raw variant will be filtered '
                     'on first place if the total depth is less then [filt_depth]. '
                     'Default %d' % defaults['filt_depth'],
            )),
            (['-V'], dict(
                dest='mean_vd',
                type='int',
                help='The filtering variant depth. Variants with depth < [mean_vd] will '
                     'be considered false positive. Default is %d (meaning at least %d reads '
                     'are needed for a variant)' % (defaults['mean_vd'], defaults['mean_vd'])
            )),

            (['-m'], dict(
                dest='maf',
                type='double',
                help='If there is MAF with frequency, it will be considered dbSNP '
                     'regardless of COSMIC. Default MAF is %f' % defaults['maf'],
            )),
            (['-s'], dict(
                dest='signal_noise',
                type='int',
                help='Signal/noise value. Default %d' % defaults['signal']
            )),
            (['-c'], dict(
                dest='control',
                help='The control sample name. Any novel or COSMIC variants passing all '
                     'above filters but also detected in Control sample will be deemed '
                     'considered false positive. Use only when there\'s control sample.'
            )),
        ],
        required_keys=['vcf'],
        file_keys=['vcf'],
        key_for_sample_name='vcf')

    check_system_resources(cnf, required=['java', 'snpsift'], optional=[])

    for key in cnf.keys():
        if key in cnf['variant_filtering'].keys():
            cnf['variant_filtering'][key] = cnf[key]
            del cnf[key]

    run_one(cnf, process_one, finalize_one)


def process_one(cnf):
    anno_vcf_fpath = filter_variants(cnf, cnf['vcf'])
    return [anno_vcf_fpath]


def finalize_one(cnf, filtered_vcf_fpath):
    if filtered_vcf_fpath:
        err('Saved filtered VCF to ' + filtered_vcf_fpath)


def filter_variants(cnf, vcf_fpath):
    err('')
    err('*' * 70)

    filt_cnf = cnf['variant_filtering']

    vcf_fpath = remove_prev_pass(cnf, vcf_fpath)

    vcf_fpath = main_filtering(cnf, filt_cnf, vcf_fpath)

    vcf_fpath = run_snpsift(cnf, filt_cnf, vcf_fpath)

    return vcf_fpath


def _parse_info(line):
    return OrderedDict(f.split('=') if '=' in f else (f, True) for f in line.split(';'))


def _build_info(d):
    return ';'.join([k + ('=' + str(v)) if v else '' for k, v in d.items()])


def _reject(tokens, val='REJECT'):
    if tokens[6] == 'PASS':
        tokens[6] = val
    return '\t'.join(tokens)


def _pass(tokens):
    return '\t'.join(tokens)


def main_filtering(cnf, filt_cnf, vcf_fpath):
    control_dict = dict()
    sample_dict = dict()
    var_dict = dict()

    def __do(l):
        if l.startswith('#'):
            return l

        tokens = l.split('\t')
        vark = '\t'.join(tokens[:4])
        d = _parse_info(tokens[7])

        def less(a, b):
            return a in d and b in filt_cnf and int(d[a]) < filt_cnf[b]

        # FILTER FIRST
        for p in ('DP', 'filt_depth'), \
                 ('QUAL', 'filt_q_mean'), \
                 ('PMEAN', 'filt_p_mean'):
            if less(*p):
                return _reject(tokens)

        # FILTER NEXT
        if (filt_cnf['control'] and 'SAMPLE' in d and
            filt_cnf['control'] == d['SAMPLE']):
            reject = False

            for p in ('QUAL', 'min_q_mean'), \
                     ('PMEAN', 'min_q_mean'), \
                     ('AF', 'min_freq'),\
                     ('MQ', 'min_mq'),\
                     ('SN', 'signal_noise'),\
                     ('VD', 'mean_vd'):
                if less(*p):
                    reject = True

            id_ = tokens[2]
            cls = 'Novel'
            if 'COSM' in id_:
                cls = 'COSMIC'
            if 'rs' in id_:
                if check_clnsig(d['CLNSIG']):
                    cls = 'ClnSNP'
                else:
                    cls = 'dbSNP'

            # so that any novel variants showed up in control won't be filtered:
            if not reject and cls == 'Novel':
                control_dict[vark] = 1
            else:
                return _reject(tokens)

        return _pass(tokens)

        # Undetermined won't count toward samples
        if ('SAMPLE' in d and
            (filt_cnf['count_undetermined'] or 'Undetermined' not in d['SAMPLE'])):

            sample_dict[d['SAMPLE']] = 1
        #     d['AF']
        #     $sample{ $d{ SAMPLE } } = 1;
        # push( @{ $var{ $vark } }, $d{ AF } );

    return iterate_file(cnf, vcf_fpath, __do)


def check_clnsig(clnsig):
    if not clnsig:
        return 0

    for c in re.split('|,', clnsig):
        if 3 < c < 7 or c == 255:
            return 1

    return -1


def remove_prev_pass(cnf, vcf_fpath):
    def __do(l):
        if l.startswith('#'):
            return l

        tokens = l.split('\t')
        if tokens[6] == 'PASS':
            tokens[6] = '.'
        return '\t'.join(tokens)

    return iterate_file(cnf, vcf_fpath, __do)


def run_snpsift(cnf, vcf_cnf, vcf_fpath):
    expression = vcf_cnf.get('expression')
    if not expression:
        return vcf_fpath

    executable = get_java_tool_cmdline(cnf, 'snpsift')
    cmdline = '{executable} filter -i PASS -f {vcf_fpath} "{expression}"'.format(**locals())
    filtered_fpath = intermediate_fname(cnf, vcf_fpath, 'filtered')
    call(cnf, cmdline, filtered_fpath)
    return filtered_fpath


if __name__ == '__main__':
    main(sys.argv[1:])





































