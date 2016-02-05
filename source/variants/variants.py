import os
from genericpath import isfile, isdir
from os.path import basename, join

import source
from source.logger import info, err, debug, critical
from source.qsub_utils import submit_job, wait_for_jobs
from source.reporting.reporting import FullReport
from source.tools_from_cnf import get_system_path, get_script_cmdline
from source.file_utils import verify_file, safe_mkdir, add_suffix, file_transaction
from source.variants.vcf_processing import verify_vcf


def run_variants(cnf, samples, main_script_name=None, mut_fpath=None):
    info('Annotating...')
    _annotate(cnf, samples)
    info()
    _summarize_varqc(cnf, output_dir=cnf.output_dir, samples=samples, caption=cnf.qc_caption)
    info('')

    mut_fname = (cnf.caller or 'variants') + '.txt'
    if mut_fpath:
        mut_fname = basename(mut_fpath)

    info('Filtering...')
    _filter(cnf, samples, mut_fname)
    info()
    info('Combining results...')
    mut_fpath = _combine_results(cnf, samples, mut_fpath or join(cnf.output_dir, mut_fname))

    info()
    info('*' * 70)
    info('Saved results to ' + mut_fpath)


def _annotate(cnf, samples):
    jobs_to_wait = []

    varannotate_cmdl = (get_script_cmdline(cnf, 'python', join('scripts', 'post', 'varannotate.py')) +
        ' --sys-cnf ' + cnf.sys_cnf +
        ' --run-cnf ' + cnf.run_cnf +
        ' --project-name ' + cnf.project_name +
       (' --reuse ' if cnf.reuse_intermediate else '') +
        ' --log-dir -' +
        ' --genome ' + cnf.genome.name +
       (' --no-check ' if cnf.no_check else '') +
        ' --qc ' +
      ((' --caller ' + cnf.caller) if cnf.caller else '')
    )

    reused = 0
    for sample in samples:
        if not sample.varannotate_dirpath:
            sample.varannotate_dirpath = join(sample.dirpath, source.varannotate_name)
        if not sample.anno_vcf_fpath:
            sample.anno_vcf_fpath = join(sample.varannotate_dirpath, add_suffix(basename(sample.vcf), 'anno'))
        output_fpath = sample.anno_vcf_fpath
        if not output_fpath.endswith('.gz'):
            output_fpath += '.gz'
        debug('Checking ' + output_fpath)
        if cnf.reuse_intermediate and isfile(output_fpath) and verify_vcf(output_fpath):
            if 'validation_match' in sample.__dict__:
                info('validation_match found for ' + sample.name + ', rerunning this sample.')
                os.system('rm -rf ' + sample.dirpath)
                safe_mkdir(sample.dirpath)
            else:
                info('Annotated results ' + output_fpath + ' exist, reusing.')
                reused += 1
                info()
                continue

        info('Annotating "' + basename(sample.vcf) + '"')
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
        jobs_to_wait.append(j)
        info()

    info()
    info('-' * 70)
    if jobs_to_wait:
        info('Submittion of annotation jobs is finished.')
        jobs_to_wait = wait_for_jobs(cnf, jobs_to_wait)
    else:
        info('No annotation jobs to submit.')
    info('Reused: ' + str(reused))
    info('')
    info('-' * 70)
    info('Done annotating')

    for j in jobs_to_wait:
        if j.is_done and not j.is_failed and not verify_vcf(j.output_fpath):
            j.is_failed = True
        if j.is_done and not j.is_failed:
            if isdir(j.work_dir):
                os.system('rm -rf ' + j.work_dir)
            else:
                err('Job was done, but ' + j.work_dir + ' does not exist')

    if any(j.is_failed for j in jobs_to_wait):
        critical('Error: ' + str(sum(1 for j in jobs_to_wait if j.is_failed)) +
                 ' annotation jobs out of ' + str(len(jobs_to_wait)) + ' are failed.')


def _filter(cnf, samples, mut_fname):
    reused = 0
    jobs_to_wait = []
    for sample in samples:
        output_dirpath = sample.varfilter_dirpath = join(sample.dirpath, source.varfilter_name)
        output_fpath = sample.mut_fpath = join(sample.varfilter_dirpath, mut_fname)
        # sample.filt_tsv_fpath = join(sample.varfilter_dirpath)

        if cnf.reuse_intermediate and isfile(output_fpath) and verify_file(output_fpath):
            info('Filtered results ' + output_fpath + ' exist, reusing.')
            reused += 1
            info()
        else:
            varfilter_py = get_script_cmdline(cnf, 'python', join('scripts', 'post', 'varfilter.py'))
            work_dir = join(cnf.work_dir, 'filt_' + sample.name)
            cmdl = ('{varfilter_py}' +
                    ' --sys-cnf ' + cnf.sys_cnf +
                    ' --run-cnf ' + cnf.run_cnf +
                    ' --log-dir -' +
                    ' --vcf {sample.anno_vcf_fpath}' +
                    ' --sample {sample.name}' +
                    ' -o {output_dirpath}' +
                    ' --output-file {sample.mut_fpath}' +
                    ' --project-name ' + cnf.project_name +
                    ' --genome {cnf.genome.name}' +
                    ' --work-dir {work_dir}' +
                   (' --reuse ' if cnf.reuse_intermediate else '') +
                  ((' --caller ' + cnf.caller) if cnf.caller else '') +
                    ' --qc' +
                   (' --no-tsv' if not cnf.tsv else '')
                ).format(**locals())
            j = submit_job(cnf, cmdl, job_name='_filt_' + sample.name,
                output_fpath=output_fpath, stdout_to_outputfile=False)
            jobs_to_wait.append(j)
            info()
    info()
    info('-' * 70)
    if jobs_to_wait:
        info('Submittion of filtering jobs is finished.')
        jobs_to_wait = wait_for_jobs(cnf, jobs_to_wait)
    else:
        info('No filtering jobs to submit.')
    info('Reused: ' + str(reused))

    info('')
    info('-' * 70)
    info('Done filtering')
    if any(j.is_failed for j in jobs_to_wait):
        critical('Error: ' + str(sum(1 for j in jobs_to_wait if j.is_failed)) + ' filtering jobs out of ' + str(len(jobs_to_wait)) + ' are failed.')


def _combine_results(cnf, samples, mut_fpath):
    if cnf.reuse_intermediate and isfile(mut_fpath) and verify_file(mut_fpath):
        info('Combined filtered results ' + mut_fpath + ' exist, reusing.')
    with file_transaction(cnf.work_dir, mut_fpath) as tx:
        with open(tx, 'w') as out:
            for var_s in samples:
                verify_file(var_s.mut_fpath, is_critical=True, description='mutations file')
                with open(var_s.mut_fpath) as f:
                    out.write(f.read())
    verify_file(mut_fpath, is_critical=True, description='final combined mutation calls')
    info('Saved all mutations to ' + mut_fpath)

    _summarize_varqc(cnf, cnf.output_dir, samples, cnf.project_name, post_filter=True)

    return mut_fpath


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
