import hashlib
import base64
from os.path import basename, join

import source
from source.logger import info
from source.bcbio.bcbio_runner import Step
from source.qsub_utils import submit_job, wait_for_jobs
from source.targetcov.bam_and_bed_utils import fix_bed_for_qualimap, prepare_beds, extract_gene_names_and_filter_exons
from source.targetcov.summarize_targetcov import summarize_targqc
from source.tools_from_cnf import get_system_path


def run_targqc(cnf, output_dir, samples, target_bed, features_bed, genes_fpath=None):
    max_threads = cnf.threads
    threads_per_sample = 1  # max(max_threads / len(samples), 1)
    summary_threads = min(len(samples), max_threads)
    info('Number of threads to run summary: ' + str(summary_threads))

    jobs_to_wait = []
    if not cnf.only_summary:
        original_target_bed = target_bed
        features_bed, features_no_genes_bed, target_bed, seq2c_bed = prepare_beds(cnf, features_bed, target_bed)
        gene_keys_set, gene_keys_list, target_bed, features_bed, features_no_genes_bed = \
            extract_gene_names_and_filter_exons(cnf, target_bed, features_bed, features_no_genes_bed)
        if not genes_fpath:
            genes_fpath = join(cnf.work_dir, 'genes.txt')
            with open(genes_fpath, 'w') as f:
                f.write('\n'.join(g + '\t' + c for g, c in gene_keys_list))

        info('*' * 70)
        info()

        step = _prep_steps(cnf, threads_per_sample, summary_threads, samples, target_bed, original_target_bed, features_bed, features_no_genes_bed, genes_fpath)

        summary_wait_for_steps = []

        for sample in samples:
            info('Processing ' + basename(sample.name))
            input_params = ''
            if sample.bam:
                input_params = ' --bam ' + sample.bam
            elif sample.l_fpath and sample.r_fpath:
                input_params = ' -1 ' + sample.l_fpath + ' -2 ' + sample.r_fpath

            j = _submit_job(cnf, step, sample.name, threads=threads_per_sample, input_params=input_params, targqc_dirpath=sample.targqc_dirpath)
            jobs_to_wait.append(j)
            summary_wait_for_steps.append(step.job_name(sample.name))

            info('Done ' + basename(sample.name))
            info()

    wait_for_jobs(cnf, jobs_to_wait)

    info('Making targqc summary')
    return summarize_targqc(cnf, summary_threads, output_dir, samples, bed_fpath=target_bed, features_fpath=features_bed)


def _prep_steps(cnf, threads_per_sample, summary_threads, samples, bed_fpath, original_bed_fpath, exons_fpath, exons_no_genes_bed, genes_fpath):
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
        ' --log-dir -'

    targetcov_params = params_for_one_sample + \
        ' -s {sample_name}' + \
        ' -o ' + join('{targqc_dirpath}') + \
        ' --work-dir ' + join(cnf.work_dir, '{sample_name}') + \
       (' --bed ' + bed_fpath if bed_fpath else '') + \
       (' --original-bed ' + original_bed_fpath if original_bed_fpath else '') + \
       (' --exons ' + exons_fpath if exons_fpath else '') + \
       (' --exons-no-genes ' + exons_no_genes_bed if exons_no_genes_bed else '') + \
       (' --genes ' + genes_fpath if genes_fpath else '') + \
       (' --reannotate ' if cnf.reannotate else '') + \
        ' --no-prep-bed' + \
        ' {input_params}'

    targetcov_step = Step(cnf, run_id,
        name=source.targetseq_name, short_name='tc',
        interpreter='python',
        script=join('scripts', 'post', 'targetcov.py'),
        paramln=targetcov_params
    )

    return targetcov_step


def _submit_job(cnf, step, sample_name='', wait_for_steps=None, threads=1, is_critical=True, **kwargs):
    tool_cmdline = get_system_path(cnf, step.interpreter, step.script, is_critical=is_critical)
    if not tool_cmdline:
        return False

    kwargs['sample_name'] = sample_name
    cmdline = tool_cmdline + ' ' + step.param_line.format(**kwargs)

    info(step.name)

    job = submit_job(cnf, cmdline,
         job_name=step.job_name(sample_name),
         wait_for_steps=wait_for_steps,
         threads=threads)

    info()
    return job