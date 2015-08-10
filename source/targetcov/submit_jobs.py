import hashlib
import base64
import sys
import os
from os.path import splitext, abspath, basename, join, isfile, dirname, pardir
from os import listdir

import source
from source.logger import info, err, warn, critical
from source.bcbio_runner import Step, fix_bed_for_qualimap
from source.file_utils import safe_mkdir, verify_file, verify_dir
from source.qsub_utils import submit_job, wait_for_jobs
from source.utils import get_system_path
from source.calling_process import call
from source.targetcov.summarize_targetcov import summarize_targqc


def run_targqc(cnf, samples, main_script_name, target_bed, exons_bed, exons_no_genes_bed, genes_fpath):
    # if not target_bed:
    #     target_bed = exons_bed
    #     info('No target_bed, using exons_bed instead')

    max_threads = cnf.threads
    threads_per_sample = 1  # max(max_threads / len(samples), 1)
    summary_threads = min(len(samples), max_threads)
    info('Number of threads to run summary: ' + str(summary_threads))

    jobs_to_wait = []
    if not cnf.only_summary:
        targetcov_step, ngscat_step, qualimap_step = \
            _prep_steps(cnf, threads_per_sample, summary_threads,
                samples, target_bed, exons_bed, exons_no_genes_bed, main_script_name)

        summary_wait_for_steps = []

        for sample in samples:
            info('Processing ' + basename(sample.bam))

            info('TargetSeq for "' + basename(sample.bam) + '"')
            j = _submit_job(cnf, targetcov_step, sample.name, threads=threads_per_sample, bam=sample.bam, sample=sample.name)
            jobs_to_wait.append(j)
            summary_wait_for_steps.append(targetcov_step.job_name(sample.name))

            # if not cnf.reuse_intermediate or not sample.ngscat_done():
            #     info('NgsCat for "' + basename(sample.bam) + '"')
            #     _submit_job(cnf, ngscat_step, sample.name, threads=threads_per_sample, bam=sample.bam, sample=sample.name, is_critical=False)
            #     summary_wait_for_steps.append(ngscat_step.job_name(sample.name))

            if qualimap_step and (not cnf.reuse_intermediate or not sample.qualimap_done()):
                info('Qualimap "' + basename(sample.bam) + '"')
                j = _submit_job(cnf, qualimap_step, sample.name, threads=threads_per_sample, bam=sample.bam, sample=sample.name, is_critical=False)
                jobs_to_wait.append(j)
                summary_wait_for_steps.append(qualimap_step.job_name(sample.name))

            info('Done ' + basename(sample.bam))
            info()

    wait_for_jobs(jobs_to_wait)

    info('Making targqc summary')
    return summarize_targqc(cnf, summary_threads, cnf.output_dir, samples, target_bed, exons_bed, genes_fpath)


def _prep_steps(cnf, threads_per_sample, summary_threads, samples, bed_fpath, exons_fpath, exons_no_genes_bed, main_script_name):
    hasher = hashlib.sha1(cnf.output_dir)
    path_hash = base64.urlsafe_b64encode(hasher.digest()[0:4])[:-2]
    run_id = path_hash + '_' + cnf.project_name

    basic_params = \
        ' --sys-cnf ' + cnf.sys_cnf + \
        ' --run-cnf ' + cnf.run_cnf

    params_for_one_sample = basic_params + \
        ' -t ' + str(threads_per_sample) + \
       (' --reuse ' if cnf.reuse_intermediate else '') + \
        ' --genome ' + cnf.genome.name + \
        ' --project-name ' + cnf.project_name + ' ' \
        ' --log-dir ' + cnf.log_dir

    targetcov_params = params_for_one_sample + \
        ' -s {sample}' + \
        ' -o ' + join(cnf.output_dir, '{sample}_' + source.targetseq_name) + \
        ' --work-dir ' + join(cnf.work_dir, '{sample}_' + source.targetseq_name) + \
        ' --bam {bam}' + \
       (' --bed ' + bed_fpath if bed_fpath else '') + \
       (' --exons ' + exons_fpath if exons_fpath else '') + \
       (' --exons-no-genes ' + exons_no_genes_bed if exons_no_genes_bed else '') + \
       (' --reannotate ' if cnf.reannotate else '') + \
       (' --dedup ' if cnf.dedup else '') + \
        ' --no-prep-bed'

    targetcov_step = Step(cnf, run_id,
        name=source.targetseq_name, short_name='tc',
        interpreter='python',
        script=join('scripts', 'post', 'targetcov.py'),
        paramln=targetcov_params
    )

    ngscat_step = None
    if bed_fpath or exons_no_genes_bed:
        ngscat_params = params_for_one_sample + \
            ' -s {sample} ' + \
            ' -o ' + join(cnf.output_dir, '{sample}_' + source.ngscat_name) + \
            ' --work-dir ' + join(cnf.work_dir, '{sample}_' + source.ngscat_name) + \
            ' --bam {bam}' + \
            ' --bed ' + (bed_fpath or exons_no_genes_bed) + \
            ' --saturation y '

        ngscat_step = Step(cnf, run_id,
            name=source.ngscat_name, short_name='nc',
            interpreter='python',
            script=join('scripts', 'post', 'ngscat.py'),
            paramln=ngscat_params
        )

    qualimap_bed_fpath = join(cnf.work_dir, 'tmp_qualimap.bed')
    if bed_fpath:
        fix_bed_for_qualimap(bed_fpath, qualimap_bed_fpath)

    qualimap_params = \
        params_for_one_sample + \
        ' --bam {bam}' + \
       (' --bed ' + qualimap_bed_fpath if bed_fpath else '') + \
        ' -o ' + join(cnf.output_dir, '{sample}_' + source.qualimap_name)

    qualimap_step = None
    if cnf.qualimap:
        qualimap_step = Step(cnf, run_id,
            name=source.qualimap_name, short_name='qm',
            interpreter='python',
            script=join('scripts', 'post', 'qualimap.py'),
            paramln=qualimap_params,
        )

    return targetcov_step, ngscat_step, qualimap_step


def _submit_job(cnf, step, sample_name='', wait_for_steps=None, threads=1, is_critical=True, **kwargs):
    tool_cmdline = get_system_path(cnf, step.interpreter, step.script, is_critical=is_critical)
    if not tool_cmdline:
        return False

    cmdline = tool_cmdline + ' ' + step.param_line.format(**kwargs)

    info(step.name)

    job = submit_job(cnf, cmdline,
        job_name=step.job_name(sample_name),
        wait_for_steps=wait_for_steps,
        threads=threads)

    info()
    return job