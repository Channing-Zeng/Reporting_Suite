import hashlib
import base64
from os.path import basename, join

import source
from source.logger import info
from source.bcbio.bcbio_runner import Step
from source.qsub_utils import submit_job, wait_for_jobs
from source.reporting.reporting import FullReport
from source.tools_from_cnf import get_system_path, get_script_cmdline
from source.file_utils import verify_file


def run_variants(cnf, samples, main_script_name):
    max_threads = cnf.threads
    threads_per_sample = 1  # max(max_threads / len(samples), 1)
    summary_threads = min(len(samples), max_threads)
    info('Number of threads to run summary: ' + str(summary_threads))

    jobs_to_wait = []
    if not cnf.only_summary:
        varannotate_cmdl = (get_script_cmdline(cnf, 'python', '/gpfs/group/ngs/src/az.reporting/scripts/post/varannotate.py') +
            ' --sys-cnf ' + cnf.sys_cnf +
            ' --run-cnf ' + cnf.run_cnf +
            ' --project-name ' + cnf.project_name +
           (' --reuse ' if cnf.reuse_intermediate else '') +
            ' --log-dir -' +
            ' --genome ' + cnf.genome.name +
           (' --no-check ' if cnf.no_check else '') +
            ' --qc '
        )

        for sample in samples:
            info('Processing ' + basename(sample.bam))

            info('TargetSeq for "' + basename(sample.bam) + '"')
            j = submit_job(cnf, varannotate_cmdl + ' --vcf ' + sample.vcf,
                           job_name='VA_' + cnf.project_name + '_' + sample.name,
                           threads=threads_per_sample, bam=sample.bam)
            jobs_to_wait.append(j)

            info('Done submitting ' + basename(sample.vcf))
            info()

    wait_for_jobs(cnf, jobs_to_wait)

    summarize_varqc(cnf, cnf.output_dir, samples, cnf.project_name)


def summarize_varqc(cnf, output_dir, samples, caption):
    info('VarQC summary...')

    jsons_by_sample = dict()
    for s in samples:
        fpath = join(s.dirpath, 'varAnnotate', 'qc', s.name + '.varQC.json')
        if verify_file(fpath):
            jsons_by_sample[s.name] = fpath

    htmls_by_sample = dict()
    for s in samples:
        fpath = join(s.dirpath, 'varAnnotate', 'qc', s.name + '.varQC.html')
        if verify_file(fpath):
            jsons_by_sample[s.name] = fpath

    report = FullReport.construct_from_sample_report_jsons(samples, output_dir, jsons_by_sample, htmls_by_sample)
    full_summary_fpaths = report.save_into_files(cnf, join(output_dir, 'varQC'), caption='Variant QC, ' + caption)

    info()
    info('*' * 70)
    for fpath in full_summary_fpaths:
        if fpath:
            info(fpath)

    return full_summary_fpaths