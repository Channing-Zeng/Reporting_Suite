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


def run_targqc(cnf, samples, main_script_name, target_bed, exons_bed, genes_fpath):
    # if not target_bed:
    #     target_bed = exons_bed
    #     info('No target_bed, using exons_bed instead')

    max_threads = cnf.threads
    threads_per_sample = 1  # max(max_threads / len(samples), 1)
    summary_threads = min(len(samples), max_threads)
    info('Number of threads to run summary: ' + str(summary_threads))

    jobs_to_wait = []
    if not cnf.only_summary:
        exons_bed, exons_no_genes_bed, target_bed, seq2c_bed = prepare_beds(cnf, exons_bed, target_bed)
        gene_names_set, gene_names_list, target_bed, exons_bed, exons_no_genes_bed = \
            extract_gene_names_and_filter_exons(cnf, target_bed, exons_bed, exons_no_genes_bed, genes_fpath)
        if not genes_fpath:
            genes_fpath = join(cnf.work_dir, 'genes.txt')
            with open(genes_fpath, 'w') as f:
                f.write('\n'.join(gene_names_list))

        info('*' * 70)
        info()

        targetcov_step, ngscat_step, qualimap_step = \
            _prep_steps(cnf, threads_per_sample, summary_threads,
                samples, target_bed, exons_bed, exons_no_genes_bed, genes_fpath, main_script_name)

        summary_wait_for_steps = []

        for sample in samples:
            info('Processing ' + basename(sample.bam))

            info('TargetSeq for "' + basename(sample.bam) + '"')
            j = _submit_job(cnf, targetcov_step, sample, threads=threads_per_sample, bam=sample.bam)
            jobs_to_wait.append(j)
            summary_wait_for_steps.append(targetcov_step.job_name(sample.name))

            # if not cnf.reuse_intermediate or not sample.ngscat_done():
            #     info('NgsCat for "' + basename(sample.bam) + '"')
            #     _submit_job(cnf, ngscat_step, sample, threads=threads_per_sample, bam=sample.bam, sample=sample.name, is_critical=False)
            #     summary_wait_for_steps.append(ngscat_step.job_name(sample.name))

            # if qualimap_step and (not cnf.reuse_intermediate or not sample.qualimap_done()):
            #     info('Qualimap "' + basename(sample.bam) + '"')
            #     j = _submit_job(cnf, qualimap_step, sample, threads=threads_per_sample, bam=sample.bam, sample=sample.name, is_critical=False)
            #     jobs_to_wait.append(j)
            #     summary_wait_for_steps.append(qualimap_step.job_name(sample.name))

            info('Done ' + basename(sample.bam))
            info()

    wait_for_jobs(cnf, jobs_to_wait)

    info('Making targqc summary')
    return summarize_targqc(cnf, summary_threads, cnf.output_dir, samples, target_bed, exons_bed, genes_fpath)


def _prep_steps(cnf, threads_per_sample, summary_threads, samples,
                bed_fpath, exons_fpath, exons_no_genes_bed, genes_fpath, main_script_name):
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
        ' -s {sample.name}' + \
        ' -o ' + join('{sample.targqc_dirpath}') + \
        ' --work-dir ' + join(cnf.work_dir, '{sample.name}') + \
        ' --bam {bam}' + \
       (' --bed ' + bed_fpath if bed_fpath else '') + \
       (' --exons ' + exons_fpath if exons_fpath else '') + \
       (' --exons-no-genes ' + exons_no_genes_bed if exons_no_genes_bed else '') + \
       (' --genes ' + genes_fpath if genes_fpath else '') + \
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
            ' -s {sample.name} ' + \
            ' -o {sample.ngscat_dirpath}' + \
            ' --work-dir ' + join(cnf.work_dir, source.ngscat_name + '_{sample.name}') + \
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
        ' -o ' + join(cnf.output_dir, '{sample.qualimap_dirpath}')

    qualimap_step = None
    if cnf.qualimap:
        qualimap_step = Step(cnf, run_id,
            name=source.qualimap_name, short_name='qm',
            interpreter='python',
            script=join('scripts', 'post', 'qualimap.py'),
            paramln=qualimap_params,
        )

    return targetcov_step, ngscat_step, qualimap_step


def _submit_job(cnf, step, sample='', wait_for_steps=None, threads=1, is_critical=True, **kwargs):
    tool_cmdline = get_system_path(cnf, step.interpreter, step.script, is_critical=is_critical)
    if not tool_cmdline:
        return False

    kwargs['sample'] = sample
    cmdline = tool_cmdline + ' ' + step.param_line.format(**kwargs)

    info(step.name)

    job = submit_job(cnf, cmdline,
        job_name=step.job_name(sample.name),
        wait_for_steps=wait_for_steps,
        threads=threads)

    info()
    return job