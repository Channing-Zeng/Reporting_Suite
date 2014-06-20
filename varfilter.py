#!/usr/bin/env python
import sys
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

import os
from os.path import basename, join, isfile
import re
import shutil
import operator
from collections import defaultdict

from source.main import read_opts_and_cnfs, check_system_resources
from source.config import Defaults
from source.logger import err, step_greetings, info, critical
from source.utils_from_bcbio import add_suffix
from source.runner import run_one
from source.utils import get_java_tool_cmdline, intermediate_fname, iterate_file, call, mean, convert_file

from ext_modules import vcf

def main():
    defaults = Defaults.variant_filtering

    cnf = read_opts_and_cnfs(
        description=
        'The program will filter an annotated VCF file by SnpEff using dbSNP and COSMIC, '
        'setting the value of the FILTER column.\n'
        '\n'
        'A novel variant (non-dbSNP, non-COSMIC) is considered false positive '
        'if all three conditions (-r -f -n) are met. False positive variants are '
        'annotated PASS in column FILTER if the conditions are satisfied, or with '
        'other value otherwise, where the value is ;-separated list of failed criteria.',

        extra_opts=[
            (['--vcf', '--var'], dict(
                dest='vcf',
                help='Annotated variants to filter')
             ),
            (['-i', '--impact'], dict(
                dest='impact',
                help='Effect impact. Default: ' + defaults['impact']
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

    #check_system_resources(cnf, required=['java', 'snpsift'], optional=[])
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

    info('Saved filtered VCF to ' + final_vcf_fpath)
    return [final_vcf_fpath]


def _add_reject(rec, val='REJECT'):
    if rec.FILTER in ['PASS', '.']:
        rec.FILTER = None

    elif val not in rec.FILTER:
        rec.add_filter(val)

    return rec


def _filter_effects(filt_cnf, d, line_num, rec):
    if 'EFF' not in d:
        critical('Error: in variant line ' + str(line_num + 1) + ', EFF field missing in INFO column')

    reject_values = []

    for f in ['impact', 'effect_type']:
        if filt_cnf[f]:
            values = filt_cnf[f].split('|')
            for v in values:
                if v.lower() not in ''.join(d['EFF']).lower():
                    reject_values.append('NO_' + v.upper())

    for val in reject_values:
        rec = _add_reject(rec, val)

    return rec


def _proc_vcf(inp_f, out_f, proc_line_fun):
    reader = vcf.Reader(inp_f)
    writer = vcf.Writer(out_f, reader)

    for i, rec in enumerate(reader):
        vark = ':'.join(map(str, [rec.CHROM, str(rec.POS), rec.REF, rec.ALT]))
        d = rec.INFO.copy()
        d['QUAL'] = rec.QUAL

        rec = proc_line_fun(rec, i, vark, d)
        if rec:
            writer.write_record(rec)


def main_filtering(cnf, filt_cnf, vcf_fpath):
    control_dict = dict()
    sample_dict = dict()
    var_dict = defaultdict(list)

    def __comp(real_key, test_key, d, line_num, op=operator.lt):
        assert test_key in filt_cnf

        if real_key not in d:
            critical('Error: in variant line ' + str(line_num + 1) + ', value ' +
                     real_key + ' missing -- requied to test ' + test_key)

        return op(float(d[real_key]), filt_cnf[test_key])

    def proc_line_1st_round(rec, i, vark, d):
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
            rec = _add_reject(rec, val)
        if reject_values:
            return rec

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

                cls = get_class(d, rec.ID)

                # So that any novel variants showed up in control won't be filtered:
                if reject_values == [] or cls == 'Novel':
                    control_dict[vark] = 1
                else:
                    for val in reject_values:
                        rec = _add_reject(rec, val)

            # Undetermined won't count toward samples
            if 'undetermined' not in d['SAMPLE'].lower() or filt_cnf['count_undetermined']:
                sample_dict[d['SAMPLE']] = 1
                var_dict[vark].append(0.0 if 'AF' not in d else float(d['AF']))

        rec = _filter_effects(filt_cnf, d, i, rec)

        return rec

    def proc_line_2nd_round(rec, i, vark, d):
        samples_n = len(sample_dict.keys())
        less = lambda x, y: __comp(x, y, d, i, op=operator.lt)
        greater = lambda x, y: __comp(x, y, d, i, op=operator.gt)

        if vark not in var_dict:  # Likely just in Undetermined
            if 'SAMPLE' in d and d['SAMPLE']:
                rec = _add_reject(rec, 'UNDET_SAMPLE')
            return rec

        var_n = len(var_dict[vark])
        average_af = mean(var_dict[vark])
        fraction = float(var_n) / samples_n

        reject_values = []
        if fraction > filt_cnf['fraction'] and var_n >= filt_cnf['sample_cnt'] \
            and average_af < filt_cnf['freq'] and rec.ID == '.':
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

        cls = get_class(d, rec.ID)
        if greater('GMAF', 'maf'):
            # if there's MAF with frequency, it'll be considered
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
            rec = _add_reject(rec, val)

        return rec

    proc_vcf_1st_round = lambda inp_f, out_f: _proc_vcf(inp_f, out_f, proc_line_1st_round)

    proc_vcf_2nd_round = lambda inp_f, out_f: _proc_vcf(inp_f, out_f, proc_line_2nd_round)

    step_greetings('Filtering')

    vcf_fpath = convert_file(cnf, vcf_fpath, proc_vcf_1st_round, suffix='zh1')

    vcf_fpath = convert_file(cnf, vcf_fpath, proc_vcf_2nd_round, suffix='zh2')

    return vcf_fpath


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
    step_greetings('Removing previous "PASS" values.')

    def proc_line(rec, i, vark, d):
        rec.FILTER = 'PASS'
        return rec

    proc_vcf = lambda inp_f, out_f: _proc_vcf(inp_f, out_f, proc_line)

    vcf_fpath = convert_file(cnf, vcf_fpath, proc_vcf, suffix='rpp')

    info('Done.')

    return vcf_fpath


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
    main()





































