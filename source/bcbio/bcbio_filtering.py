from collections import OrderedDict, defaultdict
import os
from os.path import basename, join, isfile, islink, splitext, isdir

import source
from source.bcbio.bcbio_structure import BCBioStructure
from source.logger import err, warn, send_email, critical
from source.variants.filtering import run_vcf2txt_vardict2mut_for_samples, combine_vcfs, index_vcf
from source.file_utils import safe_mkdir, add_suffix, verify_file, symlink_plus, num_lines, file_transaction
from source.logger import info


def combine_results(cnf, vcf2txt_fpath_by_sample, variants_fpath):
    info('Combining vcf2txt variants')
    not_existing_snames = []
    if cnf.reuse_intermediate and isfile(variants_fpath) and verify_file(variants_fpath):
        info('Combined filtered results ' + variants_fpath + ' exist, reusing.')
    else:
        for i, (sname, vcf2txt_fpath) in enumerate(vcf2txt_fpath_by_sample.items()):
            if not verify_file(vcf2txt_fpath, description='variants file'):
                not_existing_snames.append(sname)
        if not_existing_snames:
            critical('For some samples do not exist, variants file was not found: ' + ', '.join(not_existing_snames))
        with file_transaction(cnf.work_dir, variants_fpath) as tx:
            with open(tx, 'w') as out:
                for i, (sname, vcf2txt_fpath) in enumerate(vcf2txt_fpath_by_sample.items()):
                    with open(vcf2txt_fpath) as f:
                        for j, l in enumerate(f):
                            if j == 0 and i == 0:
                                out.write(l)
                            if j > 0:
                                out.write(l)
        verify_file(variants_fpath, is_critical=True, description='combined mutation calls')
        info('Saved vcf2txt variants to ' + variants_fpath)

    info()
    info('Combining PASSed mutations')
    pass_variants_fpath = add_suffix(variants_fpath, source.mut_pass_suffix)
    not_existing_pass_snames = []
    if cnf.reuse_intermediate and isfile(pass_variants_fpath) and verify_file(pass_variants_fpath):
        info('Combined PASSed filtered results ' + pass_variants_fpath + ' exist, reusing.')
    else:
        for i, (sname, vcf2txt_fpath) in enumerate(vcf2txt_fpath_by_sample.items()):
            if not verify_file(add_suffix(vcf2txt_fpath, source.mut_pass_suffix), description='PASS variants file'):
                not_existing_pass_snames.append(sname)
        if not_existing_pass_snames:
            critical('For some samples do not exist, PASS variants file was not found: ' + ', '.join(s.name for s in not_existing_pass_snames))
        # if cnf.variant_filtering.max_ratio_vardict2mut < 1.0:
        info('*' * 70)
        info('Max ratio set to ' + str(cnf.variant_filtering.max_ratio_vardict2mut) + ', counting freqs in cohort')
        info('Calculating frequences of varaints in the cohort')
        info('*' * 70)
        count_in_cohort_by_vark = defaultdict(int)
        total_varks = 0
        total_duplicated_count = 0
        total_records_count = 0
        for i, (sname, vcf2txt_fpath) in enumerate(vcf2txt_fpath_by_sample.items()):
            met_in_this_sample = set()
            with open(add_suffix(vcf2txt_fpath, source.mut_pass_suffix)) as f:
                for j, l in enumerate(f):
                    if j > 0:
                        fs = l.replace('\n', '').split()
                        vark = ':'.join([fs[1], fs[2], fs[4], fs[5]])
                        if vark in met_in_this_sample:
                            warn(vark + ' already met for sample ' + sname)
                            total_duplicated_count += 1
                        else:
                            met_in_this_sample.add(vark)
                            count_in_cohort_by_vark[vark] += 1
                            total_varks += 1
                        total_records_count += 1
        info('Counted ' + str(len(count_in_cohort_by_vark)) + ' different variants '
             'in ' + str(len(vcf2txt_fpath_by_sample)) + ' samples with total ' + str(total_varks) + ' records')
        info('Duplicated variants: ' + str(total_duplicated_count) + ' out of total ' + str(total_records_count) + ' records')
        if cnf.variant_filtering.max_ratio_vardict2mut < 1.0:
            info('Saving passing threshold if cohort freq < ' + str(cnf.variant_filtering.max_ratio_vardict2mut) +
                 ' to ' + pass_variants_fpath)

        freq_in_cohort_by_vark = dict()
        max_freq = 0
        max_freq_vark = 0
        for vark, count in count_in_cohort_by_vark.items():
            f = float(count) / len(vcf2txt_fpath_by_sample)
            freq_in_cohort_by_vark[vark] = f
            if f > max_freq:
                max_freq = f
                max_freq_vark = vark
        info('Maximum frequency in cohort is ' + str(max_freq) + ' of ' + max_freq_vark)
        info()

        known_variants_count = 0
        act_variants_count = 0
        good_freq_variants_count = 0
        skipped_variants_count = 0
        written_lines_count = 0
        status_col, reason_col, pcnt_sample_col = None, None, None
        with file_transaction(cnf.work_dir, pass_variants_fpath) as tx:
            with open(tx, 'w') as out:
                for i, (sname, vcf2txt_fpath) in enumerate(vcf2txt_fpath_by_sample.items()):
                    with open(add_suffix(vcf2txt_fpath, source.mut_pass_suffix)) as f:
                        for j, l in enumerate(f):
                            fs = l.replace('\n', '').split('\t')
                            if j == 0 and i == 0:
                                out.write(l)
                                status_col = fs.index('Significance')
                                reason_col = status_col + 1
                                pcnt_sample_col = fs.index('Pcnt_sample')
                            if j > 0:
                                if cnf.variant_filtering.max_ratio_vardict2mut < 1.0:
                                    fs = l.replace('\n', '').split('\t')
                                    vark = ':'.join([fs[1], fs[2], fs[4], fs[5]])
                                    if len(fs) < reason_col:
                                        print l
                                    freq = freq_in_cohort_by_vark[vark]

                                    if fs[status_col] == 'known':
                                        known_variants_count += 1
                                    elif 'act_' in fs[reason_col] or 'actionable' in fs[reason_col]:
                                        act_variants_count += 1
                                    elif freq <= cnf.variant_filtering.max_ratio_vardict2mut:
                                        good_freq_variants_count += 1
                                    else:
                                        skipped_variants_count += 1
                                        continue
                                    fs[pcnt_sample_col] = str(freq)
                                    l = '\t'.join(fs) + '\n'
                                out.write(l)
                                written_lines_count += 1
        info('Skipped variants with cohort freq >= ' + str(cnf.variant_filtering.max_ratio_vardict2mut) +
             ': ' + str(skipped_variants_count))
        info('Actionable records: ' + str(act_variants_count))
        info('Not actionable, but known records: ' + str(known_variants_count))
        info('Unknown and not actionable records with freq < ' +
             str(cnf.variant_filtering.max_ratio_vardict2mut) + ': ' + str(good_freq_variants_count))
        verify_file(pass_variants_fpath, 'PASS variants file', is_critical=True)
        info('Written ' + str(written_lines_count) + ' records to ' + pass_variants_fpath)

        variants_fpath = verify_file(variants_fpath, is_critical=True)
        pass_variants_fpath = verify_file(pass_variants_fpath, is_critical=True)

        if not_existing_snames or not_existing_pass_snames:
            return None, None

    return variants_fpath, pass_variants_fpath


def finish_filtering_for_bcbio(cnf, bcbio_structure, callers):
    email_msg = ['Variant filtering finished.']

    info('')

    info('Combining resulting mutations')
    for c in callers:
        if c.get_single_samples():
            samples = c.get_single_samples()
            vcf2txt_by_sample = {s.name: s.get_vcf2txt_by_callername(c.name) for s in samples}
            combine_results(cnf, vcf2txt_by_sample, c.single_vcf2txt_res_fpath)
        if c.get_paired_samples():
            samples = c.get_paired_samples()
            vcf2txt_by_sample = {s.name: s.get_vcf2txt_by_callername(c.name) for s in samples}
            combine_results(cnf, vcf2txt_by_sample, c.paired_vcf2txt_res_fpath)

    info()
    info('Symlinking final VCFs:')
    errory = _symlink_vcfs(callers, bcbio_structure.var_dirpath)
    info()
    info('Combining final VCFs:')
    _combine_vcfs(cnf, callers, bcbio_structure.var_dirpath)

    if any(c.single_mut_res_fpath or c.paired_mut_res_fpath for c in callers):
        info()
        info('Final variants:')
        email_msg.append('')
        email_msg.append('Combined variants for each variant caller:')
        for c in callers:
            info('  ' + c.name)
            email_msg.append('  ' + c.name)

            if c.single_vcf2txt_res_fpath:
                msg = '     Single: ' + c.single_vcf2txt_res_fpath + ', ' + str(num_lines(c.single_vcf2txt_res_fpath) - 1) + ' variants'
                info(msg)
                email_msg.append(msg)
            if c.paired_vcf2txt_res_fpath:
                msg = '     Paired: ' + c.paired_vcf2txt_res_fpath + ', ' + str(num_lines(c.paired_vcf2txt_res_fpath) - 1) + ' variants'
                info(msg)
                email_msg.append(msg)

            if c.single_mut_res_fpath:
                __symlink_mut_pass(bcbio_structure, c.single_mut_res_fpath)
                msg = '     Single PASSed: ' + c.single_mut_res_fpath + ', ' + str(num_lines(c.single_mut_res_fpath) - 1) + ' variants'
                info(msg)
                email_msg.append(msg)
            if c.paired_mut_res_fpath:
                __symlink_mut_pass(bcbio_structure, c.paired_mut_res_fpath)
                msg = '     Paired PASSed: ' + c.paired_mut_res_fpath + ', ' + str(num_lines(c.paired_mut_res_fpath) - 1) + ' variants'
                info(msg)
                email_msg.append(msg)

                if cnf.datahub_path:
                    _copy_to_datahub(cnf, c, cnf.datahub_path)

    if errory:
        err()
        err('For some samples and callers annotated VCFs could not be read:')
        for sample_name, fpath in errory:
            if not fpath:
                err('  For ' + str(sample_name) + ' VCF cannot be found')
            else:
                err('  For ' + str(sample_name) + ' VCF ' + str(fpath) + ' cannot be read')


# def vcf2txt_bcbio_structure(cnf, bcbio_structure):
#     info('Starting vcf2txt.')
#     info('-' * 70)
#
#     callers = bcbio_structure.variant_callers.values()
#     if cnf.caller:
#         try:
#             callers = [next(c for c in callers if c.name == cnf.caller)]
#         except StopIteration:
#             critical('No variant caller ' + str(cnf.caller) + ' found')
#         info('Running only for ' + callers[0].name)
#
#     for c in callers:
#         filter_for_variant_caller(cnf, c, bcbio_structure)
#
#     info('Done vcf2txt for all variant callers.')
#
#
# def filter_bcbio_structure(cnf, bcbio_structure):
#     info('Starting variant filtering.')
#     info('-' * 70)
#
#     callers = bcbio_structure.variant_callers.values()
#     if cnf.caller:
#         try:
#             callers = [next(c for c in callers if c.name == cnf.caller)]
#         except StopIteration:
#             critical('No variant caller ' + str(cnf.caller) + ' found')
#         info('Running only for ' + callers[0].name)
#
#     for c in callers:
#         filter_for_variant_caller(cnf, c, bcbio_structure)
#     info('Done filtering for all variant callers.')
#
#     global glob_cnf
#     glob_cnf = cnf
#
#     threads_num = min(len(bcbio_structure.samples) * len(callers), cnf.threads)
#     # write_vcfs(cnf, bcbio_structure.samples, vcf_fpaths, caller_name, vcf2txt_out_fpath, res, threads_num)
#
#     info('Indexing final VCFs')
#     Parallel(n_jobs=threads_num) \
#         (delayed(index_vcf)(
#                 None,
#                 sample.name,
#                 sample.get_pass_filt_vcf_fpath_by_callername(caller.name, gz=False),
#                 sample.get_filt_vcf_fpath_by_callername(caller.name, gz=False),
#                 caller.name)
#             for caller in callers for sample in caller.samples)
#
#     finish_filtering_for_bcbio(cnf, bcbio_structure, callers)


def _combine_vcfs(cnf, callers, datestamp_var_dirpath):
    for caller in callers:
        combined_vcf_fpath = join(datestamp_var_dirpath, caller.name + '.vcf')
        vcf_fpath_by_sname = {s.name: s.find_filt_vcf_by_callername(caller.name) for s in caller.samples}
        vcf_fpath_by_sname = {s_name: vcf_fpath for s_name, vcf_fpath in vcf_fpath_by_sname.items() if vcf_fpath}
        info(caller.name + ': writing to ' + combined_vcf_fpath)
        combine_vcfs(cnf, vcf_fpath_by_sname, combined_vcf_fpath)


# def filtering_cohorts(cnf, caller, bcbio_structure):
#     all_vcf_by_sample = caller.find_anno_vcf_by_sample()
#     if len(all_vcf_by_sample) == 0:
#         err('No vcfs for ' + caller.name + '. Skipping.')
#         return caller
#
#     def fill_in(batches):
#         vcf_by_sample = OrderedDict()
#         for b in batches:
#             info('Batch ' + b.name)
#             if b.normal:
#                 info('  normal sample: ' + b.normal.name)
#                 vcf_fpath = all_vcf_by_sample.get(b.normal.name)
#                 if vcf_fpath:
#                     vcf_by_sample[b.normal.name] = vcf_fpath
#                     info('  normal VCF: ' + vcf_fpath)
#             if len(b.tumor) > 1:
#                 err('  ERROR: ' + caller.name + ': ' + str(len(b.tumor)) + ' tumor samples (' + ', '.join(t.name for t in b.tumor) + ') for batch ' + b.name)
#             if len(b.tumor) > 0:
#                 info('  tumor sample: ' + b.tumor[0].name)
#                 vcf_fpath = all_vcf_by_sample.get(b.tumor[0].name)
#                 if vcf_fpath:
#                     vcf_by_sample[b.tumor[0].name] = vcf_fpath
#                     info('  tumor VCF: ' + vcf_fpath)
#         return vcf_by_sample
#
#     paired_batches = [b for b in bcbio_structure.batches.values() if b.paired]
#     info('Paired batches: ' + ', '.join(b.name for b in paired_batches))
#     paired_vcf_by_sample = fill_in(paired_batches)
#
#     single_batches = [b for b in bcbio_structure.batches.values() if not b.paired]
#     info('Single batches: ' + ', '.join(b.name for b in single_batches))
#     single_vcf_by_sample = fill_in(single_batches)
#
#     if single_vcf_by_sample:
#         info('*' * 70)
#         info('Single samples (total ' + str(len(single_vcf_by_sample)) + '):')
#
#         vcf2txt_fname = source.mut_fname_template.format(caller_name=caller.name)
#         if paired_vcf_by_sample:
#             vcf2txt_fname = add_suffix(vcf2txt_fname, source.mut_single_suffix)
#
#         vcf2txt_fpath, mut_fpath = __proc_caller_samples(cnf, bcbio_structure, caller, single_vcf_by_sample, vcf2txt_fname)
#
#         caller.single_vcf2txt_res_fpath = vcf2txt_fpath
#         caller.single_mut_res_fpath = mut_fpath
#         if mut_fpath:
#             info('Done filtering with vcf2txt/vardict2mut for single samples, result is ' + str(mut_fpath))
#
#     if paired_vcf_by_sample:
#         info('*' * 70)
#         info('Paired samples (total ' + str(len(paired_vcf_by_sample)) + '):')
#
#         vcf2txt_fname = source.mut_fname_template.format(caller_name=caller.name)
#         if single_vcf_by_sample:
#             vcf2txt_fname = add_suffix(vcf2txt_fname, source.mut_paired_suffix)
#
#         vcf2txt_fpath, mut_fpath = __proc_caller_samples(cnf, bcbio_structure, caller, paired_vcf_by_sample, vcf2txt_fname)
#
#         caller.paired_vcf2txt_res_fpath = vcf2txt_fpath
#         caller.paired_mut_res_fpath = mut_fpath
#         if mut_fpath:
#             info('Done filtering with vcf2txt/vardict2mut for paired samples, result is ' + str(mut_fpath))
#
#     info('-' * 70)
#     info()
#
#     return caller



# def vcf2txt_for_variant_caller(cnf, caller, bcbio_structure):
#     info('Running for ' + caller.name)
#
#     all_vcf_by_sample = caller.find_anno_vcf_by_sample()
#     if len(all_vcf_by_sample) == 0:
#         err('No vcfs for ' + caller.name + '. Skipping.')
#         return caller
#
#     def fill_in(batches):
#         vcf_by_sample = OrderedDict()
#         for b in batches:
#             info('Batch ' + b.name)
#             if b.normal:
#                 info('  normal sample: ' + b.normal.name)
#                 vcf_fpath = all_vcf_by_sample.get(b.normal.name)
#                 if vcf_fpath:
#                     vcf_by_sample[b.normal.name] = vcf_fpath
#                     info('  normal VCF: ' + vcf_fpath)
#             if len(b.tumor) > 1:
#                 err('  ERROR: ' + caller.name + ': ' + str(len(b.tumor)) + ' tumor samples (' + ', '.join(t.name for t in b.tumor) + ') for batch ' + b.name)
#             if len(b.tumor) > 0:
#                 info('  tumor sample: ' + b.tumor[0].name)
#                 vcf_fpath = all_vcf_by_sample.get(b.tumor[0].name)
#                 if vcf_fpath:
#                     vcf_by_sample[b.tumor[0].name] = vcf_fpath
#                     info('  tumor VCF: ' + vcf_fpath)
#         return vcf_by_sample
#
#     paired_batches = [b for b in bcbio_structure.batches.values() if b.paired]
#     info('Paired batches: ' + ', '.join(b.name for b in paired_batches))
#     paired_vcf_by_sample = fill_in(paired_batches)
#
#     single_batches = [b for b in bcbio_structure.batches.values() if not b.paired]
#     info('Single batches: ' + ', '.join(b.name for b in single_batches))
#     single_vcf_by_sample = fill_in(single_batches)
#
#     if single_vcf_by_sample:
#         info('*' * 70)
#         info('Single samples - total ' + str(len(single_vcf_by_sample)))
#
#         vcf2txt_fname = source.mut_fname_template.format(caller_name=caller.name)
#         if paired_vcf_by_sample:
#             vcf2txt_fname = add_suffix(vcf2txt_fname, source.mut_single_suffix)
#
#         vcf2txt_fpath, mut_fpath = __proc_caller_samples(cnf, bcbio_structure, caller, single_vcf_by_sample, vcf2txt_fname)
#
#         caller.single_vcf2txt_res_fpath = vcf2txt_fpath
#         caller.single_mut_res_fpath = mut_fpath
#         if mut_fpath:
#             info('Done filtering with vcf2txt/vardict2mut for single samples, result is ' + str(mut_fpath))
#
#     if paired_vcf_by_sample:
#         info('*' * 70)
#         info('Paired samples - total ' + str(len(paired_vcf_by_sample)))
#
#         vcf2txt_fname = source.mut_fname_template.format(caller_name=caller.name)
#         if single_vcf_by_sample:
#             vcf2txt_fname = add_suffix(vcf2txt_fname, source.mut_paired_suffix)
#
#         vcf2txt_fpath, mut_fpath = __proc_caller_samples(cnf, bcbio_structure, caller, paired_vcf_by_sample, vcf2txt_fname)
#
#         caller.paired_vcf2txt_res_fpath = vcf2txt_fpath
#         caller.paired_mut_res_fpath = mut_fpath
#         if mut_fpath:
#             info('Done filtering with vcf2txt/vardict2mut for paired samples, result is ' + str(mut_fpath))
#
#     info('-' * 70)
#     info()
#
#     return caller


# def filter_for_variant_caller(cnf, caller, bcbio_structure):
#     info('Running for ' + caller.name)
#
#     all_vcf_by_sample = caller.find_anno_vcf_by_sample()
#     if len(all_vcf_by_sample) == 0:
#         err('No vcfs for ' + caller.name + '. Skipping.')
#         return caller
#
#     def fill_in(batches):
#         vcf_by_sample = OrderedDict()
#         for b in batches:
#             info('Batch ' + b.name)
#             if b.normal:
#                 info('  normal sample: ' + b.normal.name)
#                 vcf_fpath = all_vcf_by_sample.get(b.normal.name)
#                 if vcf_fpath:
#                     vcf_by_sample[b.normal.name] = vcf_fpath
#                     info('  normal VCF: ' + vcf_fpath)
#             if len(b.tumor) > 1:
#                 err('  ERROR: ' + caller.name + ': ' + str(len(b.tumor)) +
#                     ' tumor samples (' + ', '.join(t.name for t in b.tumor) + ') for batch ' + b.name)
#             if len(b.tumor) > 0:
#                 info('  tumor sample: ' + b.tumor[0].name)
#                 vcf_fpath = all_vcf_by_sample.get(b.tumor[0].name)
#                 if vcf_fpath:
#                     vcf_by_sample[b.tumor[0].name] = vcf_fpath
#                     info('  tumor VCF: ' + vcf_fpath)
#         return vcf_by_sample
#
#     paired_batches = [b for b in bcbio_structure.batches.values() if b.paired]
#     info('Paired batches: ' + ', '.join(b.name for b in paired_batches))
#     paired_vcf_by_sample = fill_in(paired_batches)
#
#     single_batches = [b for b in bcbio_structure.batches.values() if not b.paired]
#     info('Single batches: ' + ', '.join(b.name for b in single_batches))
#     single_vcf_by_sample = fill_in(single_batches)
#
#     if single_vcf_by_sample:
#         info('*' * 70)
#         info('Single samples (total ' + str(len(single_vcf_by_sample)) + '):')
#
#         vcf2txt_fname = source.mut_fname_template.format(caller_name=caller.name)
#         if paired_vcf_by_sample:
#             vcf2txt_fname = add_suffix(vcf2txt_fname, source.mut_single_suffix)
#
#         vcf2txt_fpath, mut_fpath = __proc_caller_samples(cnf, bcbio_structure, caller, single_vcf_by_sample, vcf2txt_fname)
#
#         caller.single_vcf2txt_res_fpath = vcf2txt_fpath
#         caller.single_mut_res_fpath = mut_fpath
#         if mut_fpath:
#             info('Done filtering with vcf2txt/vardict2mut for single samples, result is ' + str(mut_fpath))
#
#     if paired_vcf_by_sample:
#         info('*' * 70)
#         info('Paired samples (total ' + str(len(paired_vcf_by_sample)) + '):')
#
#         vcf2txt_fname = source.mut_fname_template.format(caller_name=caller.name)
#         if single_vcf_by_sample:
#             vcf2txt_fname = add_suffix(vcf2txt_fname, source.mut_paired_suffix)
#
#         vcf2txt_fpath, mut_fpath = __proc_caller_samples(cnf, bcbio_structure, caller, paired_vcf_by_sample, vcf2txt_fname)
#
#         caller.paired_vcf2txt_res_fpath = vcf2txt_fpath
#         caller.paired_mut_res_fpath = mut_fpath
#         if mut_fpath:
#             info('Done filtering with vcf2txt/vardict2mut for paired samples, result is ' + str(mut_fpath))
#
#     info('-' * 70)
#     info()
#
#     return caller


def __proc_caller_samples(cnf, bcbio_structure, caller, vcf_by_sample, vcf2txt_fname):
    vcf2txt_fpath = join(bcbio_structure.var_dirpath, vcf2txt_fname)
    var_samples = []
    for s_name, vcf_fpath in vcf_by_sample.items():
        sample = next(s for s in caller.samples if s.name == s_name)
        var_s = source.VarSample(s_name, sample.dirpath)
        var_s.anno_vcf_fpath = vcf_fpath
        var_s.filt_vcf_fpath = sample.get_filt_vcf_fpath_by_callername(caller.name, gz=False)
        var_s.pass_filt_vcf_fpath = sample.get_pass_filt_vcf_fpath_by_callername(caller.name, gz=False)
        var_s.filt_tsv_fpath = sample.get_filt_tsv_fpath_by_callername(caller.name)
        var_s.varfilter_dirpath = join(sample.dirpath, BCBioStructure.varfilter_dir)
        var_samples.append(var_s)

    info()
    info('-' * 70)
    info('Filtering using vcf2txt...')

    mut_fpath = run_vcf2txt_vardict2mut_for_samples(
        cnf, var_samples, bcbio_structure.var_dirpath, vcf2txt_fpath,
        caller_name=caller.name, sample_min_freq=bcbio_structure.samples[0].min_af)

    if mut_fpath:
        __symlink_mut_pass(bcbio_structure, mut_fpath)
    return vcf2txt_fpath, mut_fpath


def __symlink_mut_pass(bcbio_structure, mut_fpath):
    # symlinking
    pass_txt_basefname = basename(mut_fpath)
    pass_txt_fpath_symlink = join(bcbio_structure.date_dirpath, pass_txt_basefname)

    if islink(pass_txt_fpath_symlink):
        try:
            os.unlink(pass_txt_fpath_symlink)
        except OSError:
            pass
    if isfile(pass_txt_fpath_symlink):
        try:
            os.remove(pass_txt_fpath_symlink)
        except OSError:
            pass
    try:
        symlink_plus(mut_fpath, pass_txt_fpath_symlink)
    except OSError:
        err('Cannot symlink ' + mut_fpath + ' -> ' + pass_txt_fpath_symlink)


def _copy_to_datahub(cnf, caller, datahub_dirpath):
    info('Copying to DataHub...')
    cmdl1 = 'ssh klpf990@ukapdlnx115.ukapd.astrazeneca.net \'bash -c ""\' '
    cmdl2 = 'scp {fpath} klpf990@ukapdlnx115.ukapd.astrazeneca.net:' + datahub_dirpath
    # caller.combined_filt_maf_fpath


def _symlink_vcfs(callers, datestamp_var_dirpath):
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
                    if verify_file(fpath, silent=True):
                        _symlink_to_dir(fpath, sample.dirpath)
                        # _symlink_to_dir(fpath, datestamp_var_dirpath)

            BCBioStructure.move_vcfs_to_var(sample)

    return errory


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

