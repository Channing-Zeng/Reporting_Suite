#!/usr/bin/env python

import os
from os.path import join, isfile, basename
from time import sleep

from source.calling_process import call
from source.tools_from_cnf import get_system_path
from source.logger import info
from source.file_utils import make_tmpfile, adjust_system_path


class JobRunning:
    def __init__(self, job_id, log_fpath, qsub_cmdline, done_marker,
                 output_fpath=None, sample=None):
        self.job_id = job_id
        self.log_fpath = log_fpath
        self.qsub_cmdline = qsub_cmdline
        self.done_marker = done_marker
        self.repr = job_id
        self.is_done = False
        self.output_fpath = output_fpath
        self.sample = sample


def submit_job(cnf, cmdline, job_name, wait_for_steps=None, threads=1, **kwargs):
    qsub = get_system_path(cnf, 'qsub', is_critical=True)
    bash = get_system_path(cnf, 'bash', is_critical=True)

    f, marker_fpath = make_tmpfile(cnf)
    if isfile(marker_fpath): os.remove(marker_fpath)
    job_id = job_name + '_' + basename(marker_fpath)
    log_fpath = join(cnf.log_dir, job_id + '.log')

    queue = cnf.queue
    runner_script = adjust_system_path(cnf.qsub_runner)
    hold_jid_line = '-hold_jid ' + ','.join(wait_for_steps or ['_'])
    qsub_cmdline = (
        '{qsub} -pe smp {threads} -S {bash} -q {queue} '
        '-j n -o {log_fpath} -e {log_fpath} {hold_jid_line} '
        '-N {job_id} {runner_script} {marker_fpath} "{cmdline}"'.format(**locals()))
    info(job_name)
    info(qsub_cmdline)
    job = JobRunning(job_id, log_fpath, qsub_cmdline, marker_fpath, **kwargs)
    call(cnf, qsub_cmdline, silent=True)
    return job


def wait_for_jobs(jobs):
    info()
    jobs = filter(None, jobs)
    waiting = False
    while True:
        # set flags for all done jobs
        for j in jobs:
            if not j.is_done and isfile(j.done_marker):
                j.is_done = True
                os.remove(j.done_marker)
                if waiting:
                    info('', print_date=False)
                info('Done ' + j.repr)
                waiting = False

        # check flags and wait if not all are done
        if not all(j.is_done for j in jobs):
            if not waiting:
                waiting = True
                info('Waiting for the jobs to be proccesed on a GRID (monitor with qstat). Jobs running: ' +
                     ', '.join(set([j.repr for j in jobs if not j.is_done])))
                info('', print_date=True, ending='')
            sleep(10)
            info('.', print_date=False, ending='')
        else:
            break