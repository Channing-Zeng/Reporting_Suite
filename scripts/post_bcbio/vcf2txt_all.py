#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc


import sys
import os
from traceback import format_exc
from os.path import join, basename, islink, isdir, isfile

from ext_modules.joblib import Parallel, delayed

from source.bcbio.bcbio_filtering import filter_bcbio_structure
from source.config import defaults
from source.logger import info, err, critical, warn
from source.bcbio.bcbio_structure import BCBioStructure, bcbio_summary_script_proc_params


def main():
    info(' '.join(sys.argv))
    info()

    description = '''
        The program filters an annotated VCF file by SnpEff using dbSNP and COSMIC,
        setting the value of the FILTER column.

        A novel variant (non-dbSNP, non-COSMIC) is considered false positive
        if all three conditions (-r -f -n) are met. False positive variants are
        annotated PASS in column FILTER if the conditions are satisfied, or with
        other value otherwise, where the value is ;-separated list of failed criteria.
        '''

    dfts = defaults['variant_filtering']
    extra_opts = [
        (['--caller'], dict(
            dest='caller',
            help='Variant caller name to process. If not set, processed all variant callers'
        )),

        (['--wgs'], dict(
            dest='is_wgs',
            action='store_true',
            default=False,
            help='Splits vcf2txt runs by samples, thus turns off cohort filtering'
        )),

        (['-b', '--bias'], dict(
            dest='bias',
            action='store_true',
            help='Novel or dbSNP variants with strand bias "2;1" or "2;0" '
                 'and AF < 0.3 will be considered as false positive.'
        )),

        (['-M', '--min-mq'], dict(
            dest='min_mq',
            type='float',
            help='The filtering mean mapping quality score for variants. '
                 'The raw variant will be filtered if the mean mapping quality '
                 'score is less then specified. Default %d' % dfts['min_mq'],
        )),

        (['-D', '--filt-depth'], dict(
            dest='filt_depth',
            type='int',
            help='The filtering total depth. The raw variant will be filtered '
                 'on first place if the total depth is less then [filt_depth]. '
                 'Default %s' % str(dfts['filt_depth']),
        )),

        (['-V', '--min-vd'], dict(
            dest='min_vd',
            type='int',
            help='The filtering variant depth. Variants with depth < [min_vd] will '
                 'be considered false positive. Default is %d (meaning at least %d reads '
                 'are needed for a variant)' % (dfts['min_vd'], dfts['min_vd'])
        )),

        (['-m', '--maf'], dict(
            dest='maf',
            type='float',
            help='If there is MAF with frequency, it will be considered dbSNP '
                 'regardless of COSMIC. Default MAF is %f' % dfts['maf'],
        )),

        (['-r', '--fraction'], dict(
            dest='fraction',
            type='float',
            help='When a novel variant is present in more than [fraction] '
                 'of samples and mean allele frequency is less than [freq], '
                 'it\'s considered as likely false positive. Default %f. '
                 'Used with -f and -n' % dfts['fraction'],
        )),

        (['-F', '--ave-freq'], dict(
            dest='ave_freq',
            type='float',
            help='When the average allele frequency is also below the [freq], '
                 'the variant is considered likely false positive. '
                 'Default %f. Used with -r and -n' % dfts['ave_freq'],
        )),

        (['--min-hotspot-freq'], dict(
            dest='min_hotspot_freq',
            type='float',
            help='The minimum allele frequency hotspot somatic mutations, typically lower then -f.'
                 'Default: 0.01 or half _min_freq_, whichever is less',
        )),

        (['-n'], dict(
            dest='sample_cnt',
            type='int',
            help='When the variant is detected in greater or equal [sample_cnt] '
                 'samples, the variant is considered likely false positive. '
                 'Default %d. Used with -r and -f' % dfts['sample_cnt'],
        )),

        (['-R', '--max-ratio'], dict(
            dest='max_ratio',
            type='float',
            help='When a variant is present in more than [fraction] of samples, '
                 'and AF < 0.3, it\'s considered as likely false positive, '
                 'even if it\'s in COSMIC. Default %f.' % dfts['max_ratio'],
        )),

        # This option moved to add_post_bcbio_args()
        # (['-f', '--min-freq'], dict(
        #     dest='min_freq',
        #     type='float',
        #     help='When individual allele frequency < freq for variants, '
        #          'it was considered likely false poitives. '
        #          'Default %f' % defaults['default_min_freq'],
        # )),

        (['-p'], dict(
            dest='min_p_mean',
            type='int',
            help='The minimum mean position in reads for variants.'
                 'Default %d bp' % dfts['min_p_mean'],
        )),

        (['-q'], dict(
            dest='min_q_mean',
            type='float',
            help='The minimum mean base quality phred score for variant.'
                 'Default %d' % dfts['min_q_mean'],
        )),

        (['-P'], dict(
            dest='filt_p_mean',
            type='int',
            help='The filtering mean position in reads for variants. '
                 'The raw variant will be filtered on first place if the mean '
                 'posititon is less then [filt_p_mean]. '
                 'Default %s bp' % str(dfts['filt_p_mean']),
        )),

        (['-Q'], dict(
            dest='filt_q_mean',
            type='float',
            help='The filtering mean base quality phred score for variants. '
                 'The raw variant will be filtered on first place  '
                 'if the mean quality is less then [filt_q_mean]. '
                 'Default %s' % str(dfts['filt_q_mean']),
        )),

        (['--sn'], dict(
            dest='signal_noise',
            type='int',
            help='Minimal signal/noise value. Default %d' % dfts['signal_noise']
        )),

        (['-u'], dict(
            dest='count_undetermined',
            action='store_false',
            default=True,
            help='Undeteremined won\'t be counted for the sample count.'
        )),

        (['-c', '--control'], dict(
            dest='control',
            help='The control sample name. Any novel or COSMIC variants passing all '
                 'above filters but also detected in Control sample will be deemed '
                 'considered false positive. Use only when there\'s control sample.'
        )),

        (['--datahub-path'], dict(
            dest='datahub_path',
            help='DataHub directory path to upload final MAFs and CNV (can be remote).',
        )),
    ]

    cnf, bcbio_structure = bcbio_summary_script_proc_params(
        BCBioStructure.varfilter_name,
        description=description,
        extra_opts=extra_opts)

    info('*' * 70)
    info()

    vcf2txt_all(cnf, bcbio_structure)


if __name__ == '__main__':
    main()

