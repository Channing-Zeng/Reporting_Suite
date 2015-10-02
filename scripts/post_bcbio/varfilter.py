#!/usr/bin/env python
import __check_python_version

import sys
import os
from traceback import format_exc
from joblib import Parallel, delayed
from os.path import join, pardir, basename, dirname, islink, isdir, isfile
from source.calling_process import call
from source.tools_from_cnf import get_system_path, get_java_tool_cmdline
from source.variants.filtering import filter_for_variant_caller
from source.config import defaults
from source.logger import info, err, send_email, critical, warn
from source.bcbio_structure import BCBioStructure, summary_script_proc_params
from source.file_utils import safe_mkdir, symlink_plus, file_exists, num_lines, verify_file
from source.variants.vcf_processing import igvtools_index
from source.variants.vcf_processing import bgzip_and_tabix


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

    cnf, bcbio_structure = summary_script_proc_params(
        BCBioStructure.varfilter_name,
        description=description,
        extra_opts=extra_opts)

    info('*' * 70)
    info()

    filter_all(cnf, bcbio_structure)


glob_cnf = None


def filter_all(cnf, bcbio_structure):
    info('Starting variant filtering.')
    info('-' * 70)

    callers = bcbio_structure.variant_callers.values()
    if cnf.caller:
        try:
            callers = [next(c for c in callers if c.name == cnf.caller)]
        except StopIteration:
            critical('No variant caller ' + str(cnf.caller) + ' found')
        info('Running only for ' + callers[0].name)

    for caller in callers:
        filter_for_variant_caller(caller, cnf, bcbio_structure)
    info('Done filtering for all variant callers.')

    global glob_cnf
    glob_cnf = cnf

    threads_num = min(len(bcbio_structure.samples) * len(callers), cnf.threads)
    # write_vcfs(cnf, bcbio_structure.samples, vcf_fpaths, caller_name, vcf2txt_out_fpath, res, threads_num)

    info('Indexing final VCFs')
    Parallel(n_jobs=threads_num) \
        (delayed(_index_vcf)(sample, caller.name)
            for caller in callers for sample in caller.samples)

    email_msg = ['Variant filtering finished.']
    # info('Results:')

    errory = _symlink_vcfs(callers, bcbio_structure.var_dirpath)

    _combine_vcfs(cnf, callers, bcbio_structure.var_dirpath)

    if any(c.single_mut_res_fpath or c.paired_mut_res_fpath for c in callers):
        info()
        info('Final variants:')
        email_msg.append('')
        email_msg.append('Combined variants for each variant caller:')
        for caller in callers:
            info('  ' + caller.name)
            email_msg.append('  ' + caller.name)

            if caller.single_vcf2txt_res_fpath:
                msg = '     Single: ' + caller.single_vcf2txt_res_fpath + ', ' + str(num_lines(caller.single_vcf2txt_res_fpath) - 1) + ' variants'
                info(msg)
                email_msg.append(msg)
            if caller.paired_vcf2txt_res_fpath:
                msg = '     Paired: ' + caller.paired_vcf2txt_res_fpath + ', ' + str(num_lines(caller.paired_vcf2txt_res_fpath) - 1) + ' variants'
                info(msg)
                email_msg.append(msg)

            if caller.single_mut_res_fpath:
                msg = '     Single PASSed: ' + caller.single_mut_res_fpath + ', ' + str(num_lines(caller.single_mut_res_fpath) - 1) + ' variants'
                info(msg)
                email_msg.append(msg)
            if caller.paired_mut_res_fpath:
                msg = '     Paired PASSed: ' + caller.paired_mut_res_fpath + ', ' + str(num_lines(caller.paired_mut_res_fpath) - 1) + ' variants'
                info(msg)
                email_msg.append(msg)

                if cnf.datahub_path:
                    _copy_to_datahub(cnf, caller, cnf.datahub_path)

    if errory:
        err()
        err('For some samples and callers annotated VCFs could not be read:')
        for sample_name, fpath in errory:
            if not fpath:
                err('  For ' + str(sample_name) + ' VCF cannot be found')
            else:
                err('  For ' + str(sample_name) + ' VCF ' + str(fpath) + ' cannot be read')


def _symlink_vcfs(callers, datestamp_var_dirpath):
    info()
    info('Symlinking final VCFs:')
    errory = []
    for caller in callers:
        info(caller.name)
        for sample in caller.samples:
            info(sample.name)

            filt_vcf_fpath = sample.find_filt_vcf_by_callername(caller.name)
            if not verify_file(filt_vcf_fpath):
                errory.append([caller.name, filt_vcf_fpath])
            else:
                base_filt_fpath = filt_vcf_fpath[:-3] if filt_vcf_fpath.endswith('.gz') else filt_vcf_fpath
                for fpath in [base_filt_fpath + '.gz',
                              base_filt_fpath + '.idx',
                              base_filt_fpath + '.gz.tbi']:
                    if verify_file(fpath):
                        _symlink_to_dir(fpath, sample.dirpath)
                        # _symlink_to_dir(fpath, datestamp_var_dirpath)

            BCBioStructure.move_vcfs_to_var(sample)

    return errory


def _combine_vcfs(cnf, callers, datestamp_var_dirpath):
    info()
    info('Combining final VCFs:')
    for caller in callers:
        combined_vcf_fpath = join(datestamp_var_dirpath, caller.name + '.vcf')
        info(caller.name + ': writing to ' + combined_vcf_fpath)
        gatk = get_java_tool_cmdline(cnf, 'gatk')
        cmdl = '{gatk} -T CombineVariants -R {cnf.genome.seq}'.format(**locals())
        for sample in caller.samples:
            vcf_fpath = sample.find_filt_vcf_by_callername(caller.name)
            if vcf_fpath:
                cmdl += ' --variant:' + sample.name + ' ' + sample.find_filt_vcf_by_callername(caller.name)
        cmdl += ' -o ' + combined_vcf_fpath
        res = call(cnf, cmdl, output_fpath=combined_vcf_fpath, stdout_to_outputfile=False, exit_on_error=False)
        if res:
            info('Joined VCFs for caller ' + caller.name + ', saved into ' + combined_vcf_fpath)
            if isfile(combined_vcf_fpath + '.tx.idx'):
                try:
                    os.remove(combined_vcf_fpath + '.tx.idx')
                except OSError:
                    err(format_exc())
                    info()
            bgzip_and_tabix(cnf, combined_vcf_fpath)
        else:
            warn('Could not join vcfs for caller ' + caller.name)


def _symlink_to_dir(fpath, dirpath):
    if not isdir(dirpath):
        safe_mkdir(dirpath)

    dst_path = join(dirpath, basename(fpath))

    if islink(dst_path) or isfile(dst_path):
        try:
            os.remove(dst_path)
        except OSError:
            err('Cannot symlink ' + fpath + ' -> ' + dst_path + ': cannot remove ' + dst_path)
            return

    try:
        symlink_plus(fpath, dst_path)
    except OSError:
        err('Cannot symlink ' + fpath + ' -> ' + dst_path)


def _index_vcf(sample, caller_name):
    global glob_cnf
    cnf = glob_cnf

    info()
    info(sample.name + ', ' + caller_name + ': indexing')

    pass_fpath = sample.get_pass_filt_vcf_fpath_by_callername(caller_name, gz=False)
    filt_fpath = sample.get_filt_vcf_fpath_by_callername(caller_name, gz=False)

    for fpath in [pass_fpath, filt_fpath]:
        if not cnf.reuse_intermediate and not verify_file(fpath, silent=True):
            err(fpath + ' does not exist - cannot IGV index')
        else:
            if cnf.reuse_intermediate and verify_file(fpath + '.idx', silent=True):
                info('Reusing existing ' + fpath + '.idx')
            else:
                igvtools_index(cnf, fpath)

    if not cnf.reuse_intermediate and not verify_file(filt_fpath, silent=True):
        err(filt_fpath + ' does not exist - cannot gzip and tabix')
    else:
        if cnf.reuse_intermediate and verify_file(filt_fpath + '.gz', silent=True) and verify_file(filt_fpath + '.gz.tbi', silent=True):
            info(filt_fpath + '.gz and .gz.tbi exist; reusing')
        else:
            bgzip_and_tabix(cnf, filt_fpath)


def _copy_to_datahub(cnf, caller, datahub_dirpath):
    info('Copying to DataHub...')
    cmdl1 = 'ssh klpf990@ukapdlnx115.ukapd.astrazeneca.net \'bash -c ""\' '
    cmdl2 = 'scp {fpath} klpf990@ukapdlnx115.ukapd.astrazeneca.net:' + datahub_dirpath
    # caller.combined_filt_maf_fpath


if __name__ == '__main__':
    main()

