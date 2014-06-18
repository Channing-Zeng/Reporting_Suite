#!/usr/bin/env python
from collections import OrderedDict
from genericpath import isfile
import os
from os.path import basename, join
import re
import shutil

import sys
import operator
from source.config import Defaults
from source.logger import err, step_greetings, info
from source.utils_from_bcbio import add_suffix

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from source.main import read_opts_and_cnfs, check_system_resources
from source.runner import run_one
from source.utils import get_java_tool_cmdline, intermediate_fname, iterate_file, call, mean


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
            (['--vcf', '--var'], dict(
                dest='vcf',
                help='Annotated variants to filter')
             ),

            (['-e', '--expression'], dict(
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
                type='float',
                help='When a novel variant is present in more than [fraction] '
                     'of samples and mean allele frequency is less than [freq], '
                     'it\'s considered as likely false positive. Default %f. '
                     'Used with -f and -n' % defaults['fraction'],
            )),
            (['-f'], dict(
                dest='freq',
                type='float',
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
                type='float',
                help='When a variant is present in more than [fraction] of samples, '
                     'and AF < 0.3, it\'s considered as likely false positive, '
                     'even if it\'s in COSMIC. Default %f.' % defaults['max_ratio'],
            )),

            (['-F'], dict(
                dest='min_freq',
                type='float',
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
                type='float',
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
                type='float',
                help='The filtering mean base quality phred score for variants. '
                     'The raw variant will be filtered on first place  '
                     'if the mean quality is less then [filt_q_mean]. '
                     'Default %f' % defaults['filt_q_mean'],
            )),

            (['-M'], dict(
                dest='mean_mq',
                type='float',
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
                type='float',
                help='If there is MAF with frequency, it will be considered dbSNP '
                     'regardless of COSMIC. Default MAF is %f' % defaults['maf'],
            )),
            (['--sn'], dict(
                dest='signal_noise',
                type='int',
                help='Signal/noise value. Default %d' % defaults['signal_noise']
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

    #check_system_resources(cnf, required=['java', 'snpsift'], optional=[])
    check_system_resources(cnf, required=['java'], optional=[])

    for key in cnf.keys():
        if key in cnf['variant_filtering'].keys():
            cnf['variant_filtering'][key] = cnf[key]
            del cnf[key]

    run_one(cnf, process_one, finalize_one)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


def process_one(cnf):
    filtered_vcf_fpath = filter_variants(cnf, cnf['vcf'])
    return [filtered_vcf_fpath]


def finalize_one(cnf, filtered_vcf_fpath):
    if filtered_vcf_fpath:
        info('Saved filtered VCF to ' + filtered_vcf_fpath)


def filter_variants(cnf, vcf_fpath):
    filt_cnf = cnf['variant_filtering']

    vcf_fpath = remove_prev_pass(cnf, vcf_fpath)

    vcf_fpath = main_filtering(cnf, filt_cnf, vcf_fpath)

    #vcf_fpath = run_snpsift(cnf, filt_cnf, vcf_fpath)

    final_vcf_fname = add_suffix(basename(cnf['vcf']), 'filt')
    final_vcf_fpath = join(cnf['output_dir'], final_vcf_fname)
    if isfile(final_vcf_fpath):
        os.remove(final_vcf_fpath)
    shutil.copyfile(vcf_fpath, final_vcf_fpath)

    return final_vcf_fpath


def _parse_info(line):
    return OrderedDict(f.split('=') if '=' in f else (f, True) for f in line.split(';'))


def _build_info(d):
    return ';'.join([k + ('=' + str(v)) if v else '' for k, v in d.items()])


def _reject(tokens, val='REJECT'):
    if tokens[6] == 'PASS' or tokens[6] == '.':
        tokens[6] = val
    return '\t'.join(tokens)


def _pass(tokens):
    return '\t'.join(tokens)


def main_filtering(cnf, filt_cnf, vcf_fpath):
    control_dict = dict()
    sample_dict = dict()
    var_dict = dict()
    data_dict = dict()

    def __comp(a, b, d, op=operator.lt):
        return a in d and b in filt_cnf and op(float(d[a]), filt_cnf[b])

    def __proc_line(l):
        if l.startswith('#'):
            return l

        tokens = l.split('\t')
        vark = ':'.join(tokens[0:2] + tokens[3:5])
        d = _parse_info(tokens[7])
        less = lambda x, y: __comp(x, y, d=d)

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
                     ('PMEAN', 'min_p_mean'), \
                     ('AF', 'min_freq'),\
                     ('MQ', 'min_mq'),\
                     ('SN', 'signal_noise'),\
                     ('VD', 'mean_vd'):
                if less(*p):
                    reject = True

            cls = get_class(d, tokens[2])

            # so that any novel variants showed up in control won't be filtered:
            if not reject or cls == 'Novel':
                control_dict[vark] = 1
            else:
                return _reject(tokens)

        # Undetermined won't count toward samples
        if not (filt_cnf['count_undetermined'] and ('SAMPLE' in d and 'undetermined' in d['SAMPLE'].lower())):
            sample_name = '' if 'SAMPLE' not in d else d['SAMPLE']
            sample_dict[sample_name] = 1
            if vark not in var_dict:
                var_dict[vark] = []
            var_dict[vark].append(0 if 'AF' not in d else float(d['AF']))

        # TODO: what should we do with lines without effects?
        if 'EFF' not in d:
            return _reject(tokens, val='NO_EFF')

        data_dict[vark] = d
        return _pass(tokens)

    def __post_proc_line(l):
        if l.startswith('#'):
            return l

        samples_n = len(sample_dict.keys())
        tokens = l.split('\t')
        vark = ':'.join(tokens[0:2] + tokens[3:5])
        d = _parse_info(tokens[7])
        less = lambda x, y: __comp(x, y, d=d)
        greater = lambda x, y: __comp(x, y, d=d, op=operator.gt)

        if vark not in data_dict or vark not in var_dict:
            return _pass(tokens)
        var_n = len(var_dict[vark])
        average_af = mean(var_dict[vark])
        fraction = float(var_n) / samples_n

        reject_val = 'PASS'
        if fraction > filt_cnf['fraction'] and var_n >= filt_cnf['sample_cnt'] \
            and average_af < filt_cnf['freq'] and tokens[2] == '.':
            reject_val = 'MULTI'
        if 'PSTD' in d and d['PSTD'] == 0 and \
           'BIAS' in d and not (d['BIAS'].endswith('0') or d['BIAS'].endswith('1')):
            reject_val = 'DUP'
        if fraction >= filt_cnf['max_ratio'] and 'AF' in d and float(d['AF']) < 0.3:
            reject_val = 'MAXRATE'

        for p in ('QUAL', 'min_q_mean'),\
                 ('PMEAN', 'min_p_mean'),\
                 ('MQ', 'min_mq'),\
                 ('SN', 'signal_noise'),\
                 ('AF', 'min_freq'),\
                 ('VD', 'mean_vd'):
            if less(*p):
                reject_val = p[0]

        if 'control' in filt_cnf and vark in control_dict:
            reject_val = 'CNTL'

        cls = get_class(d, tokens[2])
        if greater('GMAF', 'maf'): # if there's MAF with frequency, it'll be considered dbSNP regardless of COSMIC
            cls = 'dbSNP'
        ## Not needed in our python version of vcf2txt.pl:
        # Rescue deleterious dbSNP, such as rs80357372 (BRCA1 Q139* that is in dbSNP, but not in ClnSNP or COSMIC
        # if ( ($d->[6] eq "STOP_GAINED" || $d->[6] eq "FRAME_SHIFT") && $class eq "dbSNP" ) {
        #     my $pos = $1 if ( $d->[10] =~ /(\d+)/ );
        #     $class = "dbSNP_del" if ( $pos/$d->[11] < 0.95 );
        # }
        if 'bias' in filt_cnf and filt_cnf['bias'] and (cls == 'Novel' or cls == 'dbSNP') and \
           'BIAS' in d and (d['BIAS'] == "2;1" or d['BIAS'] == "2;0") and 'AF' in d and float(d['AF']) < 0.3:
            reject_val = 'BIAS'
        if check_clnsig(d) == -1 and cls != 'COSMIC':
            reject_val = 'NonClnSNP'
        return _reject(tokens, reject_val)

    step_greetings('Filtering based on Zhongwu\'s vcf2txt.pl.')

    vcf_fpath = iterate_file(cnf, vcf_fpath, __proc_line, suffix='zh1')
    return iterate_file(cnf, vcf_fpath, __post_proc_line, suffix='zh2')


def get_class(d, id_):
    cls = 'Novel'
    if 'COSM' in id_:
        cls = 'COSMIC'
    elif id_.startswith('rs'):
        if check_clnsig(d):
            cls = 'ClnSNP'
        else:
            cls = 'dbSNP'
    return cls


def check_clnsig(d):
    clnsig = None if 'CLNSIG' not in d else d['CLNSIG']
    if not clnsig:
        return 0

    for c in re.split('\||,', clnsig):
        if 3 < c < 7 or c == 255:
            return 1

    return -1


def remove_prev_pass(cnf, vcf_fpath):
    def __proc_line(l):
        if l.startswith('#'):
            return l

        tokens = l.split('\t')
        return _reject(tokens, val='.')

    step_greetings('Removing previous "PASS" values.')

    return iterate_file(cnf, vcf_fpath, __proc_line, suffix='rpp')


def run_snpsift(cnf, vcf_cnf, vcf_fpath):
    expression = vcf_cnf.get('expression')
    if not expression:
        return vcf_fpath

    step_greetings('Running SnpSift filter.')

    executable = get_java_tool_cmdline(cnf, 'snpsift')
    cmdline = '{executable} filter -i PASS -f {vcf_fpath} "{expression}"'.format(**locals())
    filtered_fpath = intermediate_fname(cnf, vcf_fpath, 'snpsift')
    call(cnf, cmdline, filtered_fpath)
    return filtered_fpath


if __name__ == '__main__':
    main(sys.argv[1:])





































