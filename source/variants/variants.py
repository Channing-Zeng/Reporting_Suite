import os
from collections import defaultdict
from genericpath import isfile, isdir
from os.path import basename, join

import source
from source.logger import info, err, debug, critical
from source.qsub_utils import submit_job, wait_for_jobs
from source.reporting.reporting import FullReport
from source.tools_from_cnf import get_system_path, get_script_cmdline
from source.file_utils import verify_file, safe_mkdir, add_suffix, file_transaction
from source.variants.filtering import make_vcf2txt_cmdl_params, run_vcf2txt_with_retries, count_cohort_freqs
from source.variants.vcf_processing import verify_vcf


def run_variants(cnf, samples, variants_fpath=None):
    info('Annotating...')
    _annotate(cnf, samples)
    info()
    if cnf.qc:
        _summarize_varqc(cnf, output_dir=cnf.output_dir, samples=samples, caption=cnf.qc_caption)
    info('')

    variants_fname = (cnf.caller or 'variants') + '.txt'
    if variants_fpath:
        variants_fname = basename(variants_fpath)

    info('Filtering...')
    variants_fpath = variants_fpath or join(cnf.output_dir, variants_fname)
    variants_fpath, pass_variants_fpath = _filter(cnf, samples, variants_fpath, variants_fname)

    if variants_fpath and pass_variants_fpath:
        info()
        info('*' * 70)
        info('Saved results:')
        info('  ' + variants_fpath)
        info('  ' + pass_variants_fpath)


def _annotate(cnf, samples):
    varannotate_cmdl = (get_script_cmdline(cnf, 'python', join('scripts', 'post', 'varannotate.py')) +
        ' --sys-cnf ' + cnf.sys_cnf +
        ' --run-cnf ' + cnf.run_cnf +
        ' --project-name ' + cnf.project_name +
       (' --reuse ' if cnf.reuse_intermediate else '') +
        ' --log-dir -' +
        ' --genome ' + cnf.genome.name +
       (' --no-check ' if cnf.no_check else '') +
       (' --qc ' if cnf.qc else '') +
      ((' --caller ' + cnf.caller) if cnf.caller else '')
    )

    total_reused = 0
    total_processed = 0
    total_success = 0
    total_failed = 0

    not_submitted_samples = samples
    while not_submitted_samples:
        jobs_to_wait = []
        submitted_samples = []
        reused_samples = []
        for sample in not_submitted_samples:
            if not sample.varannotate_dirpath:
                sample.varannotate_dirpath = join(sample.dirpath, source.varannotate_name)
            if not sample.anno_vcf_fpath:
                sample.anno_vcf_fpath = join(sample.varannotate_dirpath, add_suffix(basename(sample.vcf), 'anno'))
            output_fpath = sample.anno_vcf_fpath
            if not output_fpath.endswith('.gz'):
                output_fpath += '.gz'
            debug('Checking ' + output_fpath)
            if cnf.reuse_intermediate and isfile(output_fpath) and verify_vcf(output_fpath):
                info('Annotated results ' + output_fpath + ' exist, reusing.')
                reused_samples.append(sample)
                info()
                continue

            work_dir = join(cnf.work_dir, source.varannotate_name + '_' + sample.name)
            j = submit_job(
                cnf,
                cmdline=varannotate_cmdl +
                    ' --vcf ' + sample.vcf +
                    ' -o ' + sample.varannotate_dirpath +
                    ' -s ' + sample.name +
                    ' --work-dir ' + work_dir +
                    ' --output-file ' + output_fpath,
                job_name='VA_' + cnf.project_name + '_' + sample.name,
                output_fpath=output_fpath,
                stdout_to_outputfile=False,
                work_dir=work_dir
            )
            if not j.is_done:
                jobs_to_wait.append(j)
            submitted_samples.append(sample)
            if len(jobs_to_wait) >= cnf.threads:
                not_submitted_samples = [s for s in not_submitted_samples if
                                         s not in submitted_samples and
                                         s not in reused_samples]

                if not_submitted_samples:
                    info('Submitted ' + str(len(jobs_to_wait)) + ' jobs, waiting them to finish before '
                             'submitting more ' + str(len(not_submitted_samples)))
                else:
                    info('Submitted ' + str(len(jobs_to_wait)) + ' last jobs.')
                info()
                break
            info()

        info()
        info('-' * 70)
        if jobs_to_wait:
            info('Submitted ' + str(len(jobs_to_wait)) + ' jobs, waiting...')
            jobs_to_wait = wait_for_jobs(cnf, jobs_to_wait)
        else:
            info('No annotation jobs to submit.')
        info('')
        info('-' * 70)
        info('Finihsed annotating ' + str(len(jobs_to_wait)) + ' jobs')
        for j in jobs_to_wait:
            if j.is_done and not j.is_failed and not verify_vcf(j.output_fpath):
                j.is_failed = True
            if j.is_done and not j.is_failed:
                if isdir(j.work_dir):
                    os.system('rm -rf ' + j.work_dir)
                else:
                    err('Job was done, but ' + j.work_dir + ' does not exist')

        processed = sum(1 for j in jobs_to_wait if j.is_done)
        failed = sum(1 for j in jobs_to_wait if j.is_failed)
        success = sum(1 for j in jobs_to_wait if j.is_done and not j.is_failed)
        total_failed += failed
        total_reused += len(reused_samples)
        total_processed += processed
        total_success += success
        info('Reused: ' + str(len(reused_samples)))
        info('Processed: ' + str(processed))
        info('Success: ' + str(success))
        info('Failed: ' + str(failed))
        info()

        not_submitted_samples = [s for s in not_submitted_samples if
                                 s not in submitted_samples and
                                 s not in reused_samples]

    info('-' * 70)
    info('Done with all ' + str(len(samples)) + ' samples.')
    info('Total reused: ' + str(total_reused))
    info('Total processed: ' + str(total_processed))
    info('Total success: ' + str(total_success))
    info('Total failed: ' + str(total_failed))
    info()


def _filter(cnf, samples, variants_fpath, variants_fname):
    # if cohort_mode:
    #     info('Running vcf2txt.pl in cohort mode')
    #     vcf2txt = get_script_cmdline(cnf, 'perl', 'vcf2txt', is_critical=True)
    #     vcf_fpath_by_sample = {s.name: s.anno_vcf_fpath for s in samples}
    #     cmdline = vcf2txt + ' ' + make_vcf2txt_cmdl_params(cnf, vcf_fpath_by_sample)
    #     res = run_vcf2txt_with_retries(cnf, cmdline, variants_fpath)
    #     if not res:
    #         critical('Error: vcf2txt.pl crashed')

    total_reused = 0
    total_processed = 0
    total_success = 0
    total_failed = 0

    cohort_freqs_fpath = None
    # if cnf.variant_filtering.max_ratio_vardict2mut < 1.0:
    #     cohort_freqs_fpath = join(cnf.work_dir, 'cohort_freqs.tsv')
    #     info('*' * 70)
    #     info('Max ratio set to ' + str(cnf.variant_filtering.max_ratio_vardict2mut) + ', counting freqs in cohort')
    #     # cnf.variant_filtering.max_ratio < 1.0 or \
    #     # cnf.fraction < 1.0
    #     cohort_freqs_fpath = count_cohort_freqs(cnf, samples, cohort_freqs_fpath, max_ratio=cnf.variant_filtering.max_ratio_vardict2mut)
    #     info('*' * 70)
    # info()

    not_submitted_samples = samples
    while not_submitted_samples:
        reused_samples = []
        jobs_to_wait = []
        submitted_samples = []
        for sample in not_submitted_samples:
            output_dirpath = sample.varfilter_dirpath = join(sample.dirpath, source.varfilter_name)
            output_fpath = sample.variants_fpath = join(sample.varfilter_dirpath, variants_fname)
            pass_output_fpath = add_suffix(sample.variants_fpath, source.mut_pass_suffix)

            if cnf.reuse_intermediate and isfile(output_fpath) and verify_file(output_fpath) \
                    and isfile(pass_output_fpath) and verify_file(pass_output_fpath):
                info('Filtered results ' + output_fpath + ' and ' + pass_output_fpath + ' exist, reusing.')
                reused_samples.append(sample)
                info()
                continue

            varfilter_py = get_script_cmdline(cnf, 'python', join('scripts', 'post', 'varfilter.py'))
            work_dir = join(cnf.work_dir, 'filt_' + sample.name)
            cmdl = ('{varfilter_py}' +
                    ' --sys-cnf ' + cnf.sys_cnf +
                    ' --run-cnf ' + cnf.run_cnf +
                    ' --log-dir -' +
                    ' --vcf {sample.anno_vcf_fpath}' +
                    ' --sample {sample.name}' +
                    ' -o {output_dirpath}' +
                    ' --output-file {sample.variants_fpath}' +
                    ' --project-name ' + cnf.project_name +
                    ' --genome {cnf.genome.name}' +
                    ' --work-dir {work_dir}' +
                   (' --cohort-freqs {cohort_freqs_fpath}' if cohort_freqs_fpath else '') +
                   (' --reuse ' if cnf.reuse_intermediate else '') +
                  ((' --caller ' + cnf.caller) if cnf.caller else '') +
                   (' --qc' if cnf.qc else '') +
                   (' --no-tsv' if not cnf.tsv else '')
                ).format(**locals())
            j = submit_job(cnf, cmdl,
                job_name='_filt_' + sample.name,
                output_fpath=pass_output_fpath,
                stdout_to_outputfile=False,
                work_dir=work_dir)
            if not j.is_done:
                jobs_to_wait.append(j)
            jobs_to_wait.append(j)
            submitted_samples.append(sample)
            if len(jobs_to_wait) >= cnf.threads:
                not_submitted_samples = [s for s in not_submitted_samples if
                                         s not in submitted_samples and
                                         s not in reused_samples]
                if not_submitted_samples:
                    info('Submitted ' + str(len(jobs_to_wait)) + ' jobs, waiting them to finish before '
                             'submitting more ' + str(len(not_submitted_samples)))
                else:
                    info('Submitted ' + str(len(jobs_to_wait)) + ' last jobs.')
                info()
                break
            info()

        info()
        info('-' * 70)
        if jobs_to_wait:
            info('Submitted ' + str(len(jobs_to_wait)) + ' jobs, waiting...')
            jobs_to_wait = wait_for_jobs(cnf, jobs_to_wait)
        else:
            info('No filtering jobs to submit.')
        info('')
        info('-' * 70)
        info('Finihsed filtering ' + str(len(jobs_to_wait)) + ' jobs')
        for j in jobs_to_wait:
            if j.is_done and not j.is_failed and not verify_file(j.output_fpath):
                j.is_failed = True
            if j.is_done and not j.is_failed and not cnf.debug:
                if isdir(j.work_dir):
                    os.system('rm -rf ' + j.work_dir)
                else:
                    err('Job was done, but ' + j.work_dir + ' does not exist')

        processed = sum(1 for j in jobs_to_wait if j.is_done)
        failed = sum(1 for j in jobs_to_wait if j.is_failed)
        success = sum(1 for j in jobs_to_wait if j.is_done and not j.is_failed)
        total_failed += failed
        total_reused += len(reused_samples)
        total_processed += processed
        total_success += success
        info('Reused: ' + str(len(reused_samples)))
        info('Processed: ' + str(processed))
        info('Success: ' + str(success))
        info('Failed: ' + str(failed))
        info()

        not_submitted_samples = [s for s in not_submitted_samples if
                                 s not in submitted_samples and
                                 s not in reused_samples]
    info('-' * 70)
    info('Done with all ' + str(len(samples)) + ' samples.')
    info('Total reused: ' + str(total_reused))
    info('Total processed: ' + str(total_processed))
    info('Total success: ' + str(total_success))
    info('Total failed: ' + str(total_failed))
    info()

    info('Combining results...')
    variants_fpath, pass_variants_fpath = _combine_results(cnf, samples, variants_fpath)
    return variants_fpath, pass_variants_fpath


def _combine_results(cnf, samples, variants_fpath, cohort_mode=False):
    not_existing = []
    for i, s in enumerate(samples):
        if not verify_file(s.variants_fpath, description='variants file'):
            not_existing.append(s)
    if not_existing:
        err('For some samples do not exist, variants file was not found: ' + ', '.join(s.name for s in not_existing))
    else:
        with file_transaction(cnf.work_dir, variants_fpath) as tx:
            with open(tx, 'w') as out:
                for i, s in enumerate(samples):
                    with open(s.variants_fpath) as f:
                        for j, l in enumerate(f):
                            if j == 0 and i == 0:
                                out.write(l)
                            if j > 0:
                                out.write(l)
        verify_file(variants_fpath, is_critical=True, description='combined mutation calls')

    pass_variants_fpath = add_suffix(variants_fpath, source.mut_pass_suffix)
    not_existing_pass = []
    for i, s in enumerate(samples):
        if not verify_file(add_suffix(s.variants_fpath, source.mut_pass_suffix), description='PASS variants file'):
            not_existing_pass.append(s)
    if not_existing_pass:
        err('For some samples do not exist, PASS variants file was not found: ' + ', '.join(s.name for s in not_existing_pass))
    else:
        # if cnf.variant_filtering.max_ratio_vardict2mut < 1.0:
        info('*' * 70)
        info('Max ratio set to ' + str(cnf.variant_filtering.max_ratio_vardict2mut) + ', counting freqs in cohort')
        info('Calculating frequences of varaints in the cohort')
        info('*' * 70)
        count_in_sample_by_vark = defaultdict(int)
        total_varks = 0
        for i, s in enumerate(samples):
            with open(add_suffix(s.variants_fpath, source.mut_pass_suffix)) as f:
                for j, l in enumerate(f):
                    if j > 0:
                        fs = l.replace('\n', '').split()
                        vark = fs[1], fs[2], fs[4], fs[5]
                        count_in_sample_by_vark[vark] += 1
                        total_varks += 1
        info('Counted ' + str(len(count_in_sample_by_vark)) + ' different variants '
             'in ' + str(len(samples)) + ' with total ' + str(total_varks) + ' records')
        if cnf.variant_filtering.max_ratio_vardict2mut < 1.0:
            info('Saving passing threshold if cohort freq > ' + str(cnf.variant_filtering.max_ratio_vardict2mut) +
                 ' to ' + pass_variants_fpath)

        freq_in_sample_by_vark = dict()
        for vark, count in count_in_sample_by_vark.items():
            freq_in_sample_by_vark[vark] = float(count) / len(samples)

        skipped_variants = 0
        written_lines = 0
        status_col, reason_col = None, None
        with file_transaction(cnf.work_dir, pass_variants_fpath) as tx:
            with open(tx, 'w') as out:
                for i, s in enumerate(samples):
                    with open(add_suffix(s.variants_fpath, source.mut_pass_suffix)) as f:
                        for j, l in enumerate(f):
                            fs = l.replace('\n', '').split('\t')
                            if j == 0 and i == 0:
                                out.write(l)
                                status_col = fs.index('Significance')
                                reason_col = status_col + 1
                            if j > 0:
                                if cnf.variant_filtering.max_ratio_vardict2mut < 1.0:
                                    fs = l.replace('\n', '').split('\t')
                                    vark = fs[1], fs[2], fs[4], fs[5]
                                    if len(fs) < reason_col:
                                        print l
                                    if fs[status_col] != 'known' and 'act_' not in fs[reason_col] and 'actionable' not in fs[reason_col] \
                                            and count_in_sample_by_vark[vark] >= cnf.variant_filtering.max_ratio_vardict2mut:
                                        skipped_variants += 1
                                        continue
                                out.write(l)
                                written_lines += 1
        if skipped_variants:
            info('Skipped variants with cohort freq > ' + str(cnf.variant_filtering.max_ratio_vardict2mut) +
                 ': ' + str(skipped_variants))
        info('Saved ' + str(written_lines) + ' mutations to ' + pass_variants_fpath)

    variants_fpath = verify_file(variants_fpath, is_critical=True)
    pass_variants_fpath = verify_file(pass_variants_fpath, is_critical=True)

    if not_existing or not_existing_pass:
        return None, None

    if cnf.qc:
        _summarize_varqc(cnf, cnf.output_dir, samples, cnf.project_name, post_filter=True)

    return variants_fpath, pass_variants_fpath


def _summarize_varqc(cnf, output_dir, samples, caption, post_filter=False):
    name = source.varqc_name
    if post_filter:
        name = source.varqc_after_name
    varqc_dir = join(output_dir, name)
    safe_mkdir(varqc_dir)

    info('VarQC ' + ('(post-filtering) ' if post_filter else '') + 'summary, saving to ' + output_dir)

    jsons_by_sample = dict()
    for s in samples:
        fpath = join((s.varannotate_dirpath if not post_filter else s.varfilter_dirpath), 'qc',
                     s.name + (('-' + cnf.caller) if cnf.caller else '') + '.' + name + '.json')
        if verify_file(fpath):
            jsons_by_sample[s.name] = fpath

    htmls_by_sample = dict()
    for s in samples:
        fpath = join((s.varannotate_dirpath if not post_filter else s.varfilter_dirpath), 'qc',
                     s.name + (('-' + cnf.caller) if cnf.caller else '') + '.' + name + '.html')
        if verify_file(fpath):
            htmls_by_sample[s.name] = fpath

    report = FullReport.construct_from_sample_report_jsons(
        samples, output_dir, jsons_by_sample=jsons_by_sample, htmls_by_sample=htmls_by_sample)
    full_summary_fpaths = report.save_into_files(cnf, join(varqc_dir, name),
        caption='Variant QC' + (' post-varfilter' if post_filter else '') +
                ((', ' + caption) if caption else ''))

    info()
    info('*' * 70)
    for fpath in full_summary_fpaths:
        if fpath:
            info(fpath)

    return full_summary_fpaths
