#!/usr/bin/env python
import sys
from source.variants.vcf_processing import get_trasncripts_fpath
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from os.path import join, pardir, basename, dirname, abspath, realpath, islink, isdir, relpath
from site import addsitedir
source_dir = abspath(dirname(realpath(__file__)))
addsitedir(join(source_dir, 'ext_modules'))

import os
#from memory_profiler import profile

from source.main import load_genome_resources
from source.variants.filtering import filter_for_variant_caller
from source.config import Defaults
from source.logger import info, send_email, err
from source.bcbio_structure import BCBioStructure
from source.prepare_args_and_cnf import summary_script_proc_params
from source.file_utils import safe_mkdir, symlink_plus, file_exists, verify_file


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
                 'Default %f' % Defaults.default_min_freq,
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
            help='The control sample name. Any novel or COSMIC varia    nts passing all '
                 'above filters but also detected in Control sample will be deemed '
                 'considered false positive. Use only when there\'s control sample.'
        )),

        (['--datahub-path'], dict(
            dest='datahub_path',
            help='DataHub directory path to upload final MAFs and CNV (can be remote).',
        )),

        (['--transcripts'], dict(
            dest='transcripts_fpath')
         ),
    ]

    cnf, bcbio_structure = summary_script_proc_params(
        BCBioStructure.varfilter_name,
        description=description,
        extra_opts=extra_opts)

    info('*' * 70)
    info()

    get_trasncripts_fpath(cnf)
    filter_all(cnf, bcbio_structure)


#@profile
def filter_all(cnf, bcbio_structure):
    info('Starting variant filtering.')
    info('-' * 70)

    for _, caller in bcbio_structure.variant_callers.items():
        filter_for_variant_caller(caller, cnf, bcbio_structure)

    msg = ['Filtering finished.']
    info('Results:')
    for sample in bcbio_structure.samples:
        finalize_one(cnf, bcbio_structure, sample, msg)

    if any(c.combined_filt_maf_fpath for c in bcbio_structure.variant_callers.values()):
        info()
        info('Combined MAF files:')
        msg.append('')
        msg.append('Combined MAF files:')
        for caller in bcbio_structure.variant_callers.values():
            if caller.combined_filt_maf_fpath:
                info('  ' + caller.name)
                info('     ' + caller.combined_filt_maf_fpath)
                # info('     ' + caller.combined_filt_pass_maf_fpath)
                msg.append('  ' + caller.name)
                msg.append('     ' + caller.combined_filt_maf_fpath)
                # msg.append('     ' + caller.combined_filt_pass_maf_fpath)

                if cnf.datahub_path:
                    copy_to_datahub(cnf, caller, cnf.datahub_path)

        # send_email('\n'.join(msg))


def copy_to_datahub(cnf, caller, datahub_dirpath):
    info('Copying to DataHub...')
    cmdl1 = 'ssh klpf990@ukapdlnx115.ukapd.astrazeneca.net \'bash -c ""\' '
    cmdl2 = 'scp {fpath} klpf990@ukapdlnx115.ukapd.astrazeneca.net:' + datahub_dirpath
    # caller.combined_filt_maf_fpath


def symlink_to_dir(fpath, dirpath):
    if not isdir(dirpath):
        safe_mkdir(dirpath)

    dst_path = join(dirpath, basename(fpath))
    if islink(dst_path):
        os.unlink(dst_path)

    symlink_plus(fpath, dst_path)


def finalize_one(cnf, bcbio_structure, sample, msg):
    info(sample.name + ':')
    msg.append(sample.name + ':')

    for caller_name, caller in bcbio_structure.variant_callers.items():
        msg.append('  ' + caller_name)
        for fpath in [
            sample.get_filt_vcf_fpath_by_callername(caller_name),
            sample.get_filt_tsv_fpath_by_callername(caller_name),
            sample.get_filt_maf_fpath_by_callername(caller_name)]:

            if not file_exists(fpath):
                if sample.phenotype == 'normal':
                    continue
                else:
                    err('Phenotype is ' + sample.phenotype + ', and ' + fpath +
                        ' for ' + sample.name + ', ' + caller_name + ' does not exist.')
                    continue

            info('  ' + caller_name + ': ' + fpath)
            msg.append('    ' + fpath)
            symlink_to_dir(fpath, join(dirname(fpath), pardir))
            symlink_to_dir(fpath, join(bcbio_structure.date_dirpath, BCBioStructure.var_dir))


if __name__ == '__main__':
    main()















