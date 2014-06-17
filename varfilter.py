#!/usr/bin/env python

import sys
from source.config import Defaults
from source.logger import err

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from source.main import read_opts_and_cnfs, check_system_resources
from source.runner import run_one
from source.utils import get_java_tool_cmdline, intermediate_fname, iterate_file


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
                dest='undetermined',
                action='store_true',
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
                     'of samples and mean allele frequency is less than -f, '
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
                help='When indivisual allele frequency < feq for variants, '
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
            (['-s'], dict(
                dest='signal',
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

    check_system_resources(cnf,
                           required=['java', 'snpsift'],
                           optional=[])

    if 'expression' in cnf:
        cnf['variant_filtering'] = {'expression': cnf['expression']}

    run_one(cnf, process_one, finalize_one)


def process_one(cnf):
    anno_vcf_fpath = filter_variants(cnf, cnf['vcf'])
    return [anno_vcf_fpath]


def finalize_one(cnf, filtered_vcf_fpath):
    if filtered_vcf_fpath:
        err('Saved filtered VCF to ' + filtered_vcf_fpath)


def filter_variants(cnf, vcf_fpath):
    executable = get_java_tool_cmdline(cnf, 'snpsift')
    expression = cnf['variant_filtering']['expression'] or ''

    err('')
    err('*' * 70)
    vcf_fpath = iterate_file(cnf, vcf_fpath, lambda l: l.replace('\tPASS\t', '\t\t'))

    # awk_sort = subprocess.Popen( ["-c", "awk -f script.awk | sort > outfile.txt" ],
    # stdin= subprocess.PIPE, shell=True )
    # awk_sort.communicate( "input data\n" )
    # awk_sort.wait()

    cmdline = '{executable} filter -i PASS -f {vcf_fpath} "{expression}"'.format(**locals())
    filtered_fpath = intermediate_fname(cnf, vcf_fpath, 'filtered')
    # subprocess.call(cmdline, stdin=sys.stdin, stdout=)
    return filtered_fpath


if __name__ == '__main__':
    main(sys.argv[1:])