#!/usr/bin/env python

from __future__ import print_function
import sys
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from os.path import join, pardir, basename, splitext, isfile, dirname, abspath, realpath
from site import addsitedir
source_dir = abspath(dirname(realpath(__file__)))
addsitedir(join(source_dir, 'ext_modules'))

from optparse import OptionParser
from joblib import Parallel, delayed
import os
import shutil

from source.utils_from_bcbio import add_suffix, safe_mkdir
from source.variants.filtering import Filtering
from source.variants.tsv import make_tsv
from source.variants.vcf_processing import vcf_one_per_line, remove_rejected, convert_to_maf
from source.reporting import read_sample_names, get_sample_report_fpaths_for_bcbio_final_dir
from source.variants.summarize_qc import make_summary_reports
from source.config import Defaults, Config
from source.main import check_keys, check_inputs, set_up_dirs
from source.logger import info


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
    parser = OptionParser(description=description)
    parser.add_option('-d', dest='bcbio_final_dir', help='Path to bcbio-nextgen final directory (default is pwd)')
    parser.add_option('-s', dest='samples', help='List of samples (default is samples.txt in bcbio final directory)')
    parser.add_option('--vcf-suf', dest='vcf_suf', help='Suffix to choose VCF files (mutect, ensembl, freebayes, etc). Multiple comma-separated values allowed.')
    parser.add_option('-o', '--output_dir', dest='output_dir', metavar='DIR', help='output directory (or directory name in case of bcbio final dir)')

    parser.add_option('-v', dest='verbose', action='store_true', help='Verbose')
    parser.add_option('-t', dest='threads', type='int', help='Number of threads for each process')
    parser.add_option('-w', dest='overwrite', action='store_true', help='Overwrite existing results')

    parser.add_option('--runner', dest='qsub_runner', help='Bash script that takes command line as the 1st argument. This script will be submitted to GRID. Default: ' + Defaults.qsub_runner)
    parser.add_option('--work-dir', dest='work_dir', metavar='DIR')
    parser.add_option('--sys-cnf', '--sys-info', '--sys-cfg', dest='sys_cnf', default=Defaults.sys_cnf,
                      help='System configuration yaml with paths to external tools and genome resources (see default one %s)' % Defaults.sys_cnf)
    parser.add_option('--run-cnf', '--run-info', '--run-cfg', dest='run_cnf', default=Defaults.run_cnf,
                      help='Run configuration yaml (see default one %s)' % Defaults.run_cnf
    )

    defaults = Defaults.variant_filtering
    extra_opts = [
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
                 'Default %s' % str(defaults['filt_depth']),
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
            help='When individual allele frequency < freq for variants, '
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
                 'Default %s bp' % str(defaults['filt_p_mean']),
        )),

        (['-Q'], dict(
            dest='filt_q_mean',
            type='float',
            help='The filtering mean base quality phred score for variants. '
                 'The raw variant will be filtered on first place  '
                 'if the mean quality is less then [filt_q_mean]. '
                 'Default %s' % str(defaults['filt_q_mean']),
        )),

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
        ))]
    for args, kwargs in extra_opts:
        parser.add_option(*args, **kwargs)

    (opts, args) = parser.parse_args()
    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)

    cnf.name = cnf['name'] or 'varfilter_all'
    set_up_dirs(cnf)

    if not cnf.samples:
        cnf.samples = join(cnf.bcbio_final_dir, 'samples.txt')

    info('BCBio "final" dir: ' + cnf.bcbio_final_dir + ' (set with -d)')
    info('Samples: ' + cnf.samples + ' (set with -s)')

    if not check_keys(cnf, ['bcbio_final_dir', 'samples']):
        parser.print_help()
        sys.exit(1)

    if 'qsub_runner' in cnf:
        cnf.qsub_runner = join(cnf.sys_cnf, pardir, cnf.qsub_runner)

    if not check_inputs(cnf, file_keys=['samples', 'qsub_runner'], dir_keys=['bcbio_final_dir']):
        sys.exit(1)

    info()
    info('*' * 70)

    sample_names = read_sample_names(cnf['samples'])

    filter_all(cnf, sample_names)


class VariantCaller:
    def __init__(self, suf):
        self.name = suf
        self.suf = suf
        self.anno_vcf_fpaths = []
        self.anno_filt_vcf_fpaths = []


cnfs_for_samples = dict()


def filter_all(cnf, sample_names):
    varannotate_dir = 'varannotate'

    vcf_sufs = cnf['vcf_suf'].split(',')
    callers = [VariantCaller(suf) for suf in vcf_sufs]

    filt_cnf = cnf['variant_filtering']

    for caller in callers:
        info('*' * 70)
        info('Running for ' + caller.name)
        info('*' * 70)

        anno_vcf_fpaths, sample_names = get_sample_report_fpaths_for_bcbio_final_dir(
            cnf['bcbio_final_dir'], sample_names, varannotate_dir,
            '-' + caller.suf + '.anno.vcf')

        filtering = Filtering(cnf, filt_cnf, anno_vcf_fpaths)
        filt_anno_vcf_fpaths = filtering.run_filtering()

        global cnfs_for_samples
        for sname in sample_names:
            cnf_copy = cnf.copy()
            cnf_copy['name'] = sname
            cnfs_for_samples[sname] = cnf_copy

        results = Parallel(n_jobs=len(anno_vcf_fpaths))(delayed(postprocess)(sname, anno_vcf_fpath, work_filt_vcf_fpath)
            for sname, anno_vcf_fpath, work_filt_vcf_fpath in zip(sample_names, anno_vcf_fpaths, filt_anno_vcf_fpaths))

        for res in results:
            finalize_one(cnf, *res)


def postprocess(sname, anno_vcf_fpath, work_filt_vcf_fpath):
    cnf = cnfs_for_samples[sname]

    final_vcf_fpath = add_suffix(anno_vcf_fpath, 'filt').replace('varannotate', 'varfilter_3')

    safe_mkdir(dirname(final_vcf_fpath))

    file_basepath = splitext(final_vcf_fpath)[0]
    final_vcf_fpath = file_basepath + '.vcf'
    final_tsv_fpath = file_basepath + '.tsv'
    final_maf_fpath = file_basepath + '.maf'

    # Moving final VCF
    if isfile(final_vcf_fpath):
        os.remove(final_vcf_fpath)
    shutil.move(work_filt_vcf_fpath, final_vcf_fpath)
    os.symlink(final_vcf_fpath, work_filt_vcf_fpath)

    # Converting to TSV
    if work_filt_vcf_fpath and 'tsv_fields' in cnf:
        tsv_fpath = make_tsv(cnf, work_filt_vcf_fpath)

        if isfile(final_tsv_fpath):
            os.remove(final_tsv_fpath)
        shutil.move(tsv_fpath, final_tsv_fpath)
    else:
        final_tsv_fpath = None

    # Converting to MAF
    if cnf.make_maf:
        clean_filtered_vcf_fpath = remove_rejected(cnf, final_vcf_fpath)
        if clean_filtered_vcf_fpath is None:
            info('All variants are rejected.')
            final_maf_fpath = None
        else:
            maf_fpath = convert_to_maf(cnf, clean_filtered_vcf_fpath)
            if isfile(final_maf_fpath):
                os.remove(final_maf_fpath)
            shutil.move(maf_fpath, final_maf_fpath)
    else:
        final_maf_fpath = None

    return [final_vcf_fpath, final_tsv_fpath, final_maf_fpath]


def finalize_one(cnf, filtered_vcf_fpath, tsv_fpath, maf_fpath):
    if filtered_vcf_fpath:
        info('Saved VCF to ' + filtered_vcf_fpath)
    if tsv_fpath:
        info('Saved TSV to ' + tsv_fpath)
    if maf_fpath:
        info('Saved MAF (only passed) to ' + maf_fpath)


if __name__ == '__main__':
    main()

















