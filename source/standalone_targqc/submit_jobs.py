import hashlib
import base64
import sys
import os
from os.path import splitext, abspath, basename, join, isfile, dirname, pardir
from os import listdir

from source.logger import info, err, warn, critical
from source.bcbio_runner import Step, fix_bed_for_qualimap
from source.bcbio_structure import BCBioStructure, Sample
from source.file_utils import safe_mkdir, verify_file, verify_dir
from source.utils import get_system_path
from source.calling_process import call
from source.standalone_targqc.summarize import summarize_targqc


def run(cnf, bed_fpath, bam_fpaths, main_script_name):
    samples = [
        Sample(basename(splitext(bam_fpath)[0]), bam=bam_fpath, bed=bed_fpath, genome=cnf.genome.name)
            for bam_fpath in bam_fpaths]

    max_threads = cnf.threads or 40
    threads_per_sample = max(max_threads / len(samples), 1)

    if not cnf.only_summary:
        targetcov_step, ngscat_step, qualimap_step, targqc_summary_step = \
            _prep_steps(cnf, max_threads, threads_per_sample, bed_fpath, main_script_name)

        summary_wait_for_steps = []

        for sample in samples:
            info('Processing "' + basename(sample.bam) + '"')

            info('TargetSeq for "' + basename(sample.bam) + '"')
            _submit_job(cnf, targetcov_step, sample.name, threads=threads_per_sample,
                bam=sample.bam, sample=sample.name)

            info('NgsCat for "' + basename(sample.bam) + '"')
            _submit_job(cnf, ngscat_step, sample.name, threads=threads_per_sample,
                bam=sample.bam, sample=sample.name)

            info('Qualimap for "' + basename(sample.bam) + '"')
            _submit_job(cnf, qualimap_step, sample.name, threads=threads_per_sample,
                bam=sample.bam, sample=sample.name)

            summary_wait_for_steps.append(targetcov_step.job_name(sample.name))
            summary_wait_for_steps.append(ngscat_step.job_name(sample.name))
            summary_wait_for_steps.append(qualimap_step.job_name(sample.name))

        _submit_job(cnf, targqc_summary_step, wait_for_steps=summary_wait_for_steps)

    else:
        info('Making targqc summary')
        summarize_targqc(cnf, samples, bed_fpath)


def _prep_steps(cnf, max_threads, threads_per_sample, bed_fpath, main_script_name):
    hasher = hashlib.sha1(cnf.output_dir)
    path_hash = base64.urlsafe_b64encode(hasher.digest()[0:4])[:-1]
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
        ' -o ' + join(cnf.output_dir, '{sample}_' + BCBioStructure.targetseq_name) + \
        ' --work-dir ' + join(cnf.work_dir, '{sample}_' + BCBioStructure.targetseq_name) + \
        ' --bam {bam}' + \
        ' --bed ' + cnf.bed + \
       ('--exons ' + cnf.exons if cnf.exons else '')

    targetcov_step = Step(cnf, run_id,
        name=BCBioStructure.targetseq_name, short_name='tc',
        interpreter='python',
        script=join('sub_scripts', 'targetcov.py'),
        paramln=targetcov_params
    )

    ngscat_params = params_for_one_sample + \
        ' -s {sample} ' + \
        ' -o ' + join(cnf.output_dir, '{sample}_' + BCBioStructure.ngscat_name) + \
        ' --work-dir ' + join(cnf.work_dir, '{sample}_' + BCBioStructure.ngscat_name) + \
        ' --bam {bam}' + \
        ' --bed ' + cnf.bed + \
        ' --saturation y '

    ngscat_step = Step(cnf, run_id,
        name=BCBioStructure.ngscat_name, short_name='nc',
        interpreter='python',
        script=join('sub_scripts', 'ngscat.py'),
        paramln=ngscat_params
    )

    qualimap_bed_fpath = join(cnf.work_dir, 'tmp_qualimap.bed')
    fix_bed_for_qualimap(bed_fpath, qualimap_bed_fpath)

    qualimap_params = \
        ' bamqc' + \
        ' -nt ' + str(threads_per_sample) + \
        ' --java-mem-size=24G' + \
        ' -nr 5000 ' + \
        ' -bam {bam}' + \
        ' -outdir ' + join(cnf.output_dir, '{sample}_' + BCBioStructure.qualimap_name) + \
        ' -gff ' + qualimap_bed_fpath + \
        ' -c' + \
        ' -gd HUMAN'

    qualimap_step = Step(cnf, run_id,
        name=BCBioStructure.qualimap_name, short_name='qm',
        script='qualimap',
        paramln=qualimap_params
    )

    summary_cmdline_params = ' '.join(sys.argv[1:]) + ' --only-summary'

    targqc_summary_step = Step(
        cnf, run_id,
        name=BCBioStructure.targqc_name + '_summary', short_name='targqc',
        interpreter='python',
        script=main_script_name,
        paramln=summary_cmdline_params
    )

    return targetcov_step, ngscat_step, qualimap_step, targqc_summary_step


def _submit_job(cnf, step, sample_name='', wait_for_steps=None, threads=1, **kwargs):
    log_fpath = join(cnf.log_dir, (step.name + ('_' + sample_name if sample_name else '') + '.log'))

    if isfile(log_fpath):
        try:
            os.remove(log_fpath)
        except OSError:
            err('Warning: cannot remove log file ' + log_fpath + ', probably permission denied.')

    safe_mkdir(dirname(log_fpath))

    tool_cmdline = get_system_path(cnf, step.interpreter, step.script)
    if not tool_cmdline: sys.exit(1)
    cmdline = tool_cmdline + ' ' + step.param_line.format(**kwargs)

    hold_jid_line = '-hold_jid ' + ','.join(wait_for_steps or ['_'])
    job_name = step.job_name(sample_name)
    qsub = get_system_path(cnf, 'qsub')
    threads = str(threads)
    queue = cnf.queue
    runner_script = cnf.qsub_runner
    qsub_cmdline = (
        '{qsub} -pe smp {threads} -S /bin/bash -q {queue} '
        '-j n -o {log_fpath} -e {log_fpath} {hold_jid_line} '
        '-N {job_name} {runner_script} "{cmdline}"'.format(**locals()))

    info(step.name)
    info(qsub_cmdline)

    call(cnf, qsub_cmdline, silent=True)

    info()