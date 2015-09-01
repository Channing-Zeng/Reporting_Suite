#!/usr/bin/env python

import os
from os.path import join, isfile, basename
from time import sleep

from source.calling_process import call
from source.tools_from_cnf import get_system_path
from source.logger import info, err
from source.file_utils import make_tmpfile, adjust_system_path, verify_file


class JobRunning:
    def __init__(self, job_id, log_fpath, qsub_cmdline, done_marker,
                 output_fpath=None, tx_output_fpath=None, sample=None, **kwargs):
        self.job_id = job_id
        self.log_fpath = log_fpath
        self.qsub_cmdline = qsub_cmdline
        self.done_marker = done_marker
        self.repr = job_id
        self.is_done = False
        self.output_fpath = output_fpath
        self.tx_output_fpath = tx_output_fpath
        self.sample = sample
        for k, v in kwargs.items():
            self.__dict__[k] = v


def submit_job(cnf, cmdline, job_name, wait_for_steps=None, threads=1,
               output_fpath=None, stdout_to_outputfile=True, **kwargs):
    out_fpath = None
    if output_fpath:
        if cnf.reuse_intermediate and verify_file(output_fpath, silent=True):
            info(output_fpath + ' exists, reusing')
            return None
        if stdout_to_outputfile:
            out_fpath = output_fpath + '.tx'
            if isfile(out_fpath):
                os.remove(out_fpath)
        else:
            if isfile(output_fpath):
                os.remove(output_fpath)
        cmdline += ' > ' + output_fpath

    qsub = get_system_path(cnf, 'qsub', is_critical=True)
    bash = get_system_path(cnf, 'bash', is_critical=True)

    f, marker_fpath = make_tmpfile(cnf, prefix=job_name + '_' + str(cnf.project_name) + '_', suffix='.done_marker')
    if isfile(marker_fpath):
        os.remove(marker_fpath)
    job_id = basename(marker_fpath.split('.')[0])
    log_fpath = join(cnf.log_dir, job_id + '.log')

    queue = cnf.queue
    runner_script = adjust_system_path(cnf.qsub_runner)
    hold_jid_line = '-hold_jid ' + ','.join(wait_for_steps or ['_'])
    mem = threads * 15
    qsub_cmdline = (
        '{qsub} -pe smp {threads} -S {bash} -q {queue} '
        '-j n -o {log_fpath} -e {log_fpath} {hold_jid_line} '
        '-N {job_id} {runner_script} {marker_fpath} "{cmdline}"'.format(**locals()))
    info('Submitting job ' + job_id)
    info(qsub_cmdline)
    job = JobRunning(job_id, log_fpath, qsub_cmdline, marker_fpath,
        output_fpath=output_fpath, tx_output_fpath=out_fpath, **kwargs)
    call(cnf, qsub_cmdline, silent=True)
    return job


def wait_for_jobs(cnf, jobs):
    try:
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
                    if j.output_fpath:
                        if not verify_file(j.tx_output_fpath, description='j.tx_output_fpath for ' + str(j.name)):
                            err('Job ' + j.repr + ' was unsucsessful: ' + j.tx_output_fpath + ' does not exist or empty')
                        else:
                            os.rename(j.tx_output_fpath, j.output_fpath)
                            info('Done ' + j.repr + ', saved to ' + j.output_fpath)
                    else:
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
    except KeyboardInterrupt:
        qdel = get_system_path(cnf, 'qdel', is_critical=False)
        for j in jobs:
            if not j.is_done:
                if qdel:
                    call(cnf, qdel + ' ' + j.job_name, exit_on_error=False)
        raise
    return jobs