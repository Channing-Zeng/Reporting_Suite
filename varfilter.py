#!/usr/bin/env python
from collections import OrderedDict, defaultdict
from genericpath import isfile
import os
from os.path import basename, join
import re
import shutil

import sys
import operator
from source.config import Defaults
from source.logger import err, step_greetings, info, critical
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

            # (['-e', '--expression'], dict(
            #     dest='expression',
            #     help='Filtering line for SnpSift. Default is ' + defaults['expression']
            # )),
            (['-i', '--impact'], dict(
                dest='impact',
                help='Effect impact. Default: ' + defaults['impact']
            )),
            (['-e', '--effect-type'], dict(
                dest='effect_type',
                help='Effect type. Default: ' + defaults['effect_type']
            )),
            (['-b', '--bias'], dict(
                dest='bias',
                action='store_true',
                help='Novel or dbSNP variants with strand bias "2;1" or "2;0" '
                     'and AF < 0.3 will be considered as false positive.'
            )),
            (['-M', '--mean-mq'], dict(
                dest='mean_mq',
                type='float',
                help='The filtering mean mapping quality score for variants. '
                     'The raw variant will be filtered if the mean mapping quality '
                     'score is less then specified. Default %d' % defaults['mean_mq'],
            )),
            (['-D', '--filt-depth'], dict(
                dest='filt_depth',
                type='int',
                help='The filtering total depth. The raw variant will be filtered '
                     'on first place if the total depth is less then [filt_depth]. '
                     'Default %d' % defaults['filt_depth'],
            )),
            (['-V', '--mean-vd'], dict(
                dest='mean_vd',
                type='int',
                help='The filtering variant depth. Variants with depth < [mean_vd] will '
                     'be considered false positive. Default is %d (meaning at least %d reads '
                     'are needed for a variant)' % (defaults['mean_vd'], defaults['mean_vd'])
            )),

            (['-m', '--maf'], dict(
                dest='maf',
                type='float',
                help='If there is MAF with frequency, it will be considered dbSNP '
                     'regardless of COSMIC. Default MAF is %f' % defaults['maf'],
            )),

            (['-r', '--fraction'], dict(
                dest='fraction',
                type='float',
                help='When a novel variant is present in more than [fraction] '
                     'of samples and mean allele frequency is less than [freq], '
                     'it\'s considered as likely false positive. Default %f. '
                     'Used with -f and -n' % defaults['fraction'],
            )),
            (['-f', '--freq'], dict(
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

            (['-R', '--max-ratio'], dict(
                dest='max_ratio',
                type='float',
                help='When a variant is present in more than [fraction] of samples, '
                     'and AF < 0.3, it\'s considered as likely false positive, '
                     'even if it\'s in COSMIC. Default %f.' % defaults['max_ratio'],
            )),

            (['-F', '--min-freq'], dict(
                dest='min_freq',
                type='float',
                help='When individual allele frequency < feq for variants, '
                     'it was considered likely false poitives. '
                     'Default %f' % defaults['min_freq'],
            )),

            # (['-p'], dict(
            #     dest='min_p_mean',
            #     type='int',
            #     help='The minimum mean position in reads for variants.'
            #          'Default %d bp' % defaults['min_p_mean'],
            # )),
            # (['-q'], dict(
            #     dest='min_q_mean',
            #     type='float',
            #     help='The minimum mean base quality phred score for variant.'
            #          'Default %d' % defaults['min_q_mean'],
            # )),
            # (['-P'], dict(
            #     dest='filt_p_mean',
            #     type='int',
            #     help='The filtering mean position in reads for variants. '
            #          'The raw variant will be filtered on first place if the mean '
            #          'posititon is less then [filt_p_mean]. '
            #          'Default %d bp' % defaults['filt_p_mean'],
            # )),
            # (['-Q'], dict(
            #     dest='filt_q_mean',
            #     type='float',
            #     help='The filtering mean base quality phred score for variants. '
            #          'The raw variant will be filtered on first place  '
            #          'if the mean quality is less then [filt_q_mean]. '
            #          'Default %f' % defaults['filt_q_mean'],
            # )),

            (['--sn'], dict(
                dest='signal_noise',
                type='int',
                help='Signal/noise value. Default %d' % defaults['signal_noise']
            )),

            (['-u'], dict(
                dest='count_undetermined',
                action='store_false',
                help='Undeteremined won\'t be counted for the sample count.'
            )),

            (['-c', '--control'], dict(
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
    check_system_resources(cnf, required=['java'], optional=[])

    for key in cnf.keys():
        if key in cnf['variant_filtering'].keys():
            cnf['variant_filtering'][key] = cnf[key]
            del cnf[key]

    run_one(cnf, process_one, finalize_one)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


def finalize_one(cnf, filtered_vcf_fpath):
    if filtered_vcf_fpath:
        info('Saved filtered VCF to ' + filtered_vcf_fpath)


def process_one(cnf):
    vcf_fpath = cnf['vcf']
    filt_cnf = cnf['variant_filtering']

    vcf_fpath = remove_previous_filter_rejects(cnf, vcf_fpath)

    # vcf_fpath = run_snpsift(cnf, filt_cnf, vcf_fpath)

    vcf_fpath = main_filtering(cnf, filt_cnf, vcf_fpath)

    final_vcf_fname = add_suffix(basename(cnf['vcf']), 'filt')
    final_vcf_fpath = join(cnf['output_dir'], final_vcf_fname)
    if isfile(final_vcf_fpath):
        os.remove(final_vcf_fpath)
    shutil.copyfile(vcf_fpath, final_vcf_fpath)

    return [final_vcf_fpath]


def _parse_fields(tokens):
    d = OrderedDict(f.split('=') if '=' in f else (f, True) for f in tokens[7].split(';'))
    d['QUAL'] = tokens[5]
    return d


def _build_info(d):
    return ';'.join([k + ('=' + str(v)) if v else '' for k, v in d.items()])


def _add_reject(tokens, val='REJECT'):
    if tokens[6] in ['PASS', '.']:
        tokens[6] = val

    elif val not in tokens[6]:
        tokens[6] += ';' + val

    return tokens


def _make_var_line(tokens):
    return '\t'.join(tokens)


def _filter_effects(filt_cnf, d, i, tokens):
    if 'EFF' not in d:
        critical('Warning: in line ' + str(i + 1) + ', EFF field missing in INFO column')

    reject_values = []

    if filt_cnf['impact']:
        if filt_cnf['impact'] not in d['EFF']:
            reject_values.append('IMPACT')

    if filt_cnf['effect_type']:
        if filt_cnf['effect_type'] not in d['EFF']:
            reject_values.append('EFF_TYPE')

    for val in reject_values:
        tokens = _add_reject(tokens, val)
    return tokens


def main_filtering(cnf, filt_cnf, vcf_fpath):
    control_dict = dict()
    sample_dict = dict()
    var_dict = defaultdict(list)

    def __comp(real_key, test_key, d, line_num, op=operator.lt):
        assert test_key in filt_cnf

        if real_key not in d:
            critical('Warning: in line ' + str(line_num + 1) + ', value ' +
                     real_key + ' missing -- requied to test ' + test_key)

        return op(float(d[real_key]), filt_cnf[test_key])


    def __proc_line(l, i):
        if l.startswith('#'):
            return l

        tokens = l.split('\t')
        vark = ':'.join(tokens[0:2] + tokens[3:5])
        d = _parse_fields(tokens)
        less = lambda x, y: __comp(x, y, d, i)

        reject_values = []

        # FILTER FIRST
        for real_key, test_key in [
               ('DP', 'filt_depth'),
               ('QUAL', 'filt_q_mean'),
               # ('PMEAN', 'filt_p_mean')\
            ]:
            if less(real_key, test_key):
                reject_values.append(test_key.upper())

        for val in reject_values:
            tokens = _add_reject(tokens, val)
        if reject_values:
            return _make_var_line(tokens)

        # FILTER NEXT
        if 'SAMPLE' in d and d['SAMPLE']:
            if filt_cnf['control'] and filt_cnf['control'] == d['SAMPLE']:
                reject_values = []

                for real_key, test_key in [
                       ('QUAL', 'min_q_mean'),
                       # ('PMEAN', 'min_p_mean'),
                       ('AF', 'min_freq'),
                       ('MQ', 'min_mq'),
                       # ('SN', 'signal_noise'),
                       # ('VD', 'mean_vd')
                    ]:
                    if less(real_key, test_key):
                        reject_values.append(test_key.upper())

                cls = get_class(d, tokens[2])

                # So that any novel variants showed up in control won't be filtered:
                if reject_values == [] or cls == 'Novel':
                    control_dict[vark] = 1
                else:
                    for val in reject_values:
                        tokens = _add_reject(tokens, val)

            # Undetermined won't count toward samples
            if 'undetermined' not in d['SAMPLE'].lower() or filt_cnf['count_undetermined']:
                sample_dict[d['SAMPLE']] = 1
                var_dict[vark].append(0.0 if 'AF' not in d else float(d['AF']))

        tokens = _filter_effects(filt_cnf, d, i, tokens)

        return _make_var_line(tokens)

    def __post_proc_line(l, i):
        if l.startswith('#'):
            return l

        samples_n = len(sample_dict.keys())
        tokens = l.split('\t')
        vark = ':'.join(tokens[0:2] + tokens[3:5])
        d = _parse_fields(tokens[7])
        less = lambda x, y: __comp(x, y, d, i, op=operator.lt)
        greater = lambda x, y: __comp(x, y, d, i, op=operator.gt)

        if vark not in var_dict:  # Likely just in Undetermined
            if 'SAMPLE' in d and d['SAMPLE']:
                _add_reject(tokens, 'UNDET_SAMPLE')
            return _make_var_line(tokens)

        var_n = len(var_dict[vark])
        average_af = mean(var_dict[vark])
        fraction = float(var_n) / samples_n

        reject_values = []
        if fraction > filt_cnf['fraction'] and var_n >= filt_cnf['sample_cnt'] \
            and average_af < filt_cnf['freq'] and tokens[2] == '.':
            reject_values.append('MULTI')

        if 'PSTD' in d and d['PSTD'] == 0 and \
           'BIAS' in d and not (d['BIAS'].endswith('0') or d['BIAS'].endswith('1')):
            reject_values.append('DUP')

        if fraction >= filt_cnf['max_ratio'] and 'AF' in d and float(d['AF']) < 0.3:
            reject_values.append('MAXRATE')

        for real_key, test_key in [
                 ('QUAL', 'min_q_mean'),
                 # ('PMEAN', 'min_p_mean'),
                 ('MQ', 'min_mq'),
                 # ('SN', 'signal_noise'),
                 ('AF', 'min_freq'),
                 # ('VD', 'mean_vd')
                 ]:
            if less(real_key, test_key):
                reject_values.append(test_key.upper())

        if 'control' in filt_cnf and vark in control_dict:
            reject_values.append('CNTL')

        cls = get_class(d, tokens[2])
        if greater('GMAF', 'maf'):  # if there's MAF with frequency, it'll be considered
                                    # dbSNP regardless of COSMIC
            cls = 'dbSNP'
        ## Not needed in our python version of vcf2txt.pl:
        # Rescue deleterious dbSNP, such as rs80357372 (BRCA1 Q139* that is in dbSNP, but not in ClnSNP or COSMIC
        # if ( ($d->[6] eq "STOP_GAINED" || $d->[6] eq "FRAME_SHIFT") && $class eq "dbSNP" ) {
        #     my $pos = $1 if ( $d->[10] =~ /(\d+)/ );
        #     $class = "dbSNP_del" if ( $pos/$d->[11] < 0.95 );
        # }
        if 'bias' in filt_cnf and filt_cnf['bias'] and (cls == 'Novel' or cls == 'dbSNP') and \
           'BIAS' in d and (d['BIAS'] == "2;1" or d['BIAS'] == "2;0") and 'AF' in d and float(d['AF']) < 0.3:
            reject_values.append('BIAS')

        if check_clnsig(d) == -1 and cls != 'COSMIC':
            reject_values.append('NonClnSNP')

        for val in reject_values:
            tokens = _add_reject(tokens, val)

        return _make_var_line(tokens)

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


def remove_previous_filter_rejects(cnf, vcf_fpath):
    def __proc_line(l, i):
        if l.startswith('#'):
            return l

        tokens = l.split('\t')
        tokens[6] = 'PASS'
        return _make_var_line(tokens)

    step_greetings('Removing previous "PASS" values.')

    out_fpath = iterate_file(cnf, vcf_fpath, __proc_line, suffix='rpp')

    info('Done.')

    return out_fpath


def run_snpsift(cnf, vcf_cnf, vcf_fpath):
    expression = vcf_cnf.get('expression')
    if not expression:
        return vcf_fpath

    step_greetings('Running SnpSift filter.')

    executable = get_java_tool_cmdline(cnf, 'snpsift')
    cmdline = '{executable} filter -a EXPR -n -p -f {vcf_fpath} "{expression}"'.format(**locals())
    filtered_fpath = intermediate_fname(cnf, vcf_fpath, 'snpsift')
    call(cnf, cmdline, filtered_fpath)

    info('Done.')

    return filtered_fpath


if __name__ == '__main__':
    main(sys.argv[1:])





































