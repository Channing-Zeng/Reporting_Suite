#!/usr/bin/env python

from __future__ import print_function
import sys
from source.bcbio_structure import BCBioStructure, load_bcbio_cnf, VariantCaller
from source.summary import _check_args
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from os.path import join, pardir, basename, splitext, isfile, dirname, abspath, realpath, islink
from site import addsitedir
source_dir = abspath(dirname(realpath(__file__)))
addsitedir(join(source_dir, 'ext_modules'))

from optparse import OptionParser
from joblib import Parallel, delayed
import os
import shutil

from source.file_utils import safe_mkdir, add_suffix
from source.variants.filtering import Filtering
from source.variants.tsv import make_tsv
from source.variants.vcf_processing import remove_rejected, convert_to_maf, vcf_is_empty, \
    igvtools_index
from source.config import Defaults, Config
from source.main import set_up_dirs
from source.logger import info, critical


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
    parser.add_option('-v', dest='verbose', action='store_true', help='Verbose')
    parser.add_option('-t', dest='threads', type='int', help='Number of threads for each process')
    parser.add_option('-w', dest='overwrite', action='store_true', help='Overwrite existing results')
    parser.add_option('--reuse', dest='overwrite', help='Reuse intermediate files from previous run', action='store_false')

    parser.add_option('--runner', dest='qsub_runner', help='Bash script that takes command line as the 1st argument. This script will be submitted to GRID. Default: ' + Defaults.qsub_runner)
    parser.add_option('--sys-cnf', '--sys-info', '--sys-cfg', dest='sys_cnf', default=Defaults.sys_cnf, help='System configuration yaml with paths to external tools and genome resources (see default one %s)' % Defaults.sys_cnf)
    parser.add_option('--run-cnf', '--run-info', '--run-cfg', dest='run_cnf', default=Defaults.run_cnf, help='Run configuration yaml (see default one %s)' % Defaults.run_cnf
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
            help='Minimal signal/noise value. Default %d' % defaults['signal_noise']
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

        (['--vcf-dir'], dict(
            dest='vcf_dir',
        ))]

    for args, kwargs in extra_opts:
        parser.add_option(*args, **kwargs)

    (opts, args) = parser.parse_args()
    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)
    if not opts.bcbio_final_dir and len(args) > 0:
        cnf.bcbio_final_dir = args[0]
    else:
        critical('Usage: ./post_for_bcbio.py <final_dir>')

    _check_args(parser, cnf)

    load_bcbio_cnf(cnf)

    bcbio_structure = BCBioStructure(cnf, cnf.bcbio_final_dir, cnf.bcbio_cnf)
    cnf.work_dir = bcbio_structure.work_dir
    cnf.name = BCBioStructure.varfilter_name
    set_up_dirs(cnf)

    info()
    info('*' * 70)

    filter_all(cnf, bcbio_structure)


cnfs_for_samples = dict()


def filter_all(cnf, bcbio_structure):
    callers = bcbio_structure.variant_callers

    filt_cnf = cnf['variant_filtering']

    for caller_name, caller in callers.items():
        info('*' * 70)
        info('Running for ' + caller.name)
        info('*' * 70)

        anno_vcf_by_sample = caller.get_anno_vcf_by_samples()

        filtering = Filtering(cnf, filt_cnf, bcbio_structure, caller)
        filt_anno_vcf_fpaths = filtering.run_filtering(anno_vcf_by_sample.values())

        global cnfs_for_samples
        for sample in anno_vcf_by_sample.keys():
            cnf_copy = cnf.copy()
            cnf_copy['name'] = sample.name
            cnfs_for_samples[sample.name] = cnf_copy

        results = Parallel(n_jobs=len(caller.samples)) \
            (delayed(postprocess)
             (sample.name, anno_vcf_by_sample[sample], work_filt_vcf_fpath)
              for sample, work_filt_vcf_fpath in
              zip(caller.samples, filt_anno_vcf_fpaths
            ))

        for res in results:
            finalize_one(cnf, *res)


def postprocess(sname, anno_vcf_fpath, work_filt_vcf_fpath):
    cnf = cnfs_for_samples[sname]

    final_vcf_fpath = add_suffix(anno_vcf_fpath, 'filt').replace('varAnnotate', 'varFilter')

    safe_mkdir(dirname(final_vcf_fpath))

    file_basepath = splitext(final_vcf_fpath)[0]
    final_vcf_fpath = file_basepath + '.vcf'
    final_clean_vcf_fpath = file_basepath + '.passed.vcf'
    final_tsv_fpath = file_basepath + '.tsv'
    final_clean_tsv_fpath = file_basepath + '.passed.tsv'
    final_clean_maf_fpath = file_basepath + '.passed.maf'

    # Moving final VCF
    if isfile(final_vcf_fpath):
        os.remove(final_vcf_fpath)
    shutil.move(work_filt_vcf_fpath, final_vcf_fpath)
    os.symlink(final_vcf_fpath, work_filt_vcf_fpath)
    igvtools_index(cnf, final_vcf_fpath)

    # Cleaning rejected variants
    clean_filtered_vcf_fpath = remove_rejected(cnf, work_filt_vcf_fpath)
    if vcf_is_empty(cnf, clean_filtered_vcf_fpath):
        info('All variants are rejected.')
    if isfile(final_clean_vcf_fpath):
        os.remove(final_clean_vcf_fpath)
    shutil.move(clean_filtered_vcf_fpath, final_clean_vcf_fpath)
    os.symlink(final_clean_vcf_fpath, clean_filtered_vcf_fpath)
    igvtools_index(cnf, final_clean_vcf_fpath)

    # Converting to TSV
    if work_filt_vcf_fpath and 'tsv_fields' in cnf:
        tsv_fpath = make_tsv(cnf, work_filt_vcf_fpath)

        if isfile(final_tsv_fpath):
            os.remove(final_tsv_fpath)
        shutil.move(tsv_fpath, final_tsv_fpath)
    else:
        final_tsv_fpath = None

    # Converting clean VCF to TSV
    if clean_filtered_vcf_fpath and 'tsv_fields' in cnf:
        clean_tsv_fpath = make_tsv(cnf, clean_filtered_vcf_fpath)

        if isfile(final_clean_tsv_fpath):
            os.remove(final_clean_tsv_fpath)
        shutil.move(clean_tsv_fpath, final_clean_tsv_fpath)
    else:
        final_clean_tsv_fpath = None

    # Converting to MAF
    if clean_filtered_vcf_fpath and cnf.make_maf:
        maf_fpath = convert_to_maf(cnf, clean_filtered_vcf_fpath)
        if isfile(final_clean_maf_fpath):
            os.remove(final_clean_maf_fpath)
        shutil.move(maf_fpath, final_clean_maf_fpath)
    else:
        final_clean_maf_fpath = None

    return [final_vcf_fpath, final_clean_vcf_fpath, final_tsv_fpath, final_clean_tsv_fpath, final_clean_maf_fpath]


def finalize_one(cnf, vcf_fpath, clean_vcf_fpath, tsv_fpath, clean_tsv_fpath, maf_fpath):
    if vcf_fpath:
        info('Saved VCF to ' + vcf_fpath)
    if clean_vcf_fpath:
        info('Saved VCF (only passed) to ' + clean_vcf_fpath)
    if tsv_fpath:
        info('Saved TSV to ' + clean_tsv_fpath)
    if tsv_fpath:
        info('Saved TSV (only passed) to ' + tsv_fpath)
    if maf_fpath:
        info('Saved MAF (only passed) to ' + maf_fpath)

    if cnf.make_soft_links:
        for fpath in [vcf_fpath, tsv_fpath, maf_fpath]:
            sl_path = join(dirname(fpath), pardir, basename(fpath))
            if islink(sl_path):
                os.unlink(sl_path)
            os.symlink(fpath, sl_path)

if __name__ == '__main__':
    main()

















