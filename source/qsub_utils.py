#!/usr/bin/env python
from datetime import datetime

import os
from os.path import join, isfile, basename, splitext
from time import sleep

from source.calling_process import call
from source.tools_from_cnf import get_system_path
from source.logger import info, err, warn, timestamp
from source.file_utils import make_tmpfile, adjust_system_path, verify_file
from source.utils import is_us


class JobRunning:
    def __init__(self, job_id, log_fpath, qsub_cmdline, done_marker_fpath, error_marker_fpath,
                 output_fpath=None, tx_output_fpath=None, stdout_to_outputfile=True, sample=None, **kwargs):
        self.job_id = job_id
        self.log_fpath = log_fpath
        self.qsub_cmdline = qsub_cmdline
        self.done_marker = done_marker_fpath
        self.error_marker = error_marker_fpath
        self.repr = job_id
        self.is_done = False
        self.is_failed = False
        self.output_fpath = output_fpath
        self.tx_output_fpath = tx_output_fpath
        self.stdout_to_outputfile = stdout_to_outputfile
        self.sample = sample
        for k, v in kwargs.items():
            self.__dict__[k] = v


def submit_job(cnf, cmdline, job_name, wait_for_steps=None, threads=1,
               output_fpath=None, stdout_to_outputfile=True, run_on_chara=False, **kwargs):
    tx_output_fpath = None
    if output_fpath:
        if cnf.reuse_intermediate and verify_file(output_fpath, silent=True):
            info(output_fpath + ' exists, reusing')
            return None
        if stdout_to_outputfile:
            tx_output_fpath = output_fpath + '.tx'
            if isfile(tx_output_fpath):
                os.remove(tx_output_fpath)
            cmdline += ' > ' + tx_output_fpath
        else:
            if isfile(output_fpath):
                os.remove(output_fpath)

    qsub = get_system_path(cnf, 'qsub', is_critical=True)
    bash = get_system_path(cnf, 'bash', is_critical=True)

    prefix = str(cnf.project_name) + '_'
    if job_name:
        prefix += job_name + '_'
    prefix += datetime.now().strftime("%Y_%m_%d_%H_%M_%S") + '_'
    f, done_marker_fpath = make_tmpfile(cnf, prefix=prefix, suffix='.done')
    f, error_marker_fpath = make_tmpfile(cnf, prefix=prefix, suffix='.error')
    if isfile(done_marker_fpath): os.remove(done_marker_fpath)
    if isfile(error_marker_fpath): os.remove(error_marker_fpath)
    job_id = basename(splitext(done_marker_fpath)[0])
    if cnf.log_dir:
        err_fpath = log_fpath = join(cnf.log_dir, job_id + '.log')
    else:
        fd, fpath = make_tmpfile(cnf, suffix=job_id + '.log', text=True)
        err_fpath = log_fpath = fpath

    queue = cnf.queue
    runner_script = adjust_system_path(cnf.qsub_runner)
    verify_file(runner_script, is_critical=True, description='qsub_runner')
    hold_jid_line = '-hold_jid ' + ','.join(wait_for_steps or ['_'])
    mem = threads * 15
    priority = 0
    if cnf.qsub_priority:
        priority = cnf.qsub_priority
    extra_qsub_opts = ''
    if run_on_chara and is_us():
        extra_qsub_opts += '-l h="chara|rask"'
    cmdline = cmdline.replace('"', '\\"').replace('\\\\"', '\\"')
    qsub_cmdline = (
        '{qsub} -pe smp {threads} {extra_qsub_opts} -S {bash} -q {queue} -p {priority} '
        '-j n -o {log_fpath} -e {err_fpath} {hold_jid_line} '
        '-N {job_id} {runner_script} {done_marker_fpath} {error_marker_fpath} "{cmdline}"'.format(**locals()))
    info('Submitting job ' + job_id)
    info(qsub_cmdline)
    job = JobRunning(job_id, log_fpath, qsub_cmdline, done_marker_fpath, error_marker_fpath,
                     output_fpath=output_fpath, tx_output_fpath=tx_output_fpath,
                     stdout_to_outputfile=stdout_to_outputfile, **kwargs)
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
                        if j.stdout_to_outputfile:
                            if not verify_file(j.tx_output_fpath, description='j.tx_output_fpath for ' + str(j.repr)):
                                err('Job ' + j.repr + ' was unsuccessful: ' + str(j.tx_output_fpath) + ' does not exist or empty.' +
                                   ((' Log saved to ' + j.log_fpath) if j.log_fpath else ''))
                                j.is_failed = True
                            else:
                                os.rename(j.tx_output_fpath, j.output_fpath)
                                info('Done ' + j.repr + ', saved to ' + j.output_fpath +
                                    ((' Log saved to ' + j.log_fpath) if j.log_fpath else ''))
                        else:
                            if not verify_file(j.output_fpath, description='j.output_fpath for ' + str(j.repr)):
                                err('Job ' + j.repr + ' was unsuccessful: ' + str(j.output_fpath) + ' does not exist or empty.' +
                                   ((' Log saved to ' + j.log_fpath) if j.log_fpath else ''))
                                j.is_failed = True
                            else:
                                info('Done ' + j.repr + ', saved to ' + j.output_fpath +
                                    ((' Log saved to ' + j.log_fpath) if j.log_fpath else ''))

                    else:
                        info('Done ' + j.repr + ((' Log saved to ' + j.log_fpath) if j.log_fpath else ''))

                    waiting = False
                if not j.is_done and isfile(j.error_marker):
                    j.is_done = True
                    j.is_failed = True
                    err('Job ' + j.repr + ' returned non-0. ' +
                       ((' Log saved to ' + j.log_fpath) if j.log_fpath else ''))
                    if waiting:
                        info('', print_date=False)

            # check flags and wait if not all are done
            if not all(j.is_done for j in jobs):
                if not waiting:
                    waiting = True
                    info('Waiting for the jobs to be processed on a GRID (monitor with qstat). '
                         'Jobs running: ' + str(len(sorted([j.repr for j in jobs if not j.is_done]))))
                    info('', print_date=True, ending='')
                sleep(10)
                info('.', print_date=False, ending='')
            else:
                break
    except KeyboardInterrupt:
        not_done = sum(1 for j in jobs if not j.is_done)
        done = sum(1 for j in jobs if j.is_done)
        failed = sum(1 for j in jobs if j.is_failed)
        warn('Interrupted. Done: ' + str(done) + ', not done: ' + str(not_done) + ', failed: ' + str(failed))
    except SystemExit:
        not_done = sum(1 for j in jobs if not j.is_done)
        done = sum(1 for j in jobs if j.is_done)
        failed = sum(1 for j in jobs if j.is_failed)
        warn('System exit. Done: ' + str(done) + ', not done: ' + str(not_done) + ', failed: ' + str(failed))
    except:
        not_done = sum(1 for j in jobs if not j.is_done)
        done = sum(1 for j in jobs if j.is_done)
        failed = sum(1 for j in jobs if j.is_failed)
        warn('Exception. Done: ' + str(done) + ', not done: ' + str(not_done) + ', failed: ' + str(failed))
        raise
    finally:
        del_jobs(cnf, jobs)

    not_done = sum(1 for j in jobs if not j.is_done)
    done = sum(1 for j in jobs if j.is_done)
    failed = sum(1 for j in jobs if j.is_failed)
    (warn if failed else info)('Done: ' + str(done) + ', not done: ' + str(not_done) + ', failed: ' + str(failed))
    return jobs


def del_jobs(cnf, jobs_running):
    done_job_ids = [j.job_id for j in jobs_running if not j.is_done and not j.not_wait]
    if done_job_ids:
        qdel = get_system_path(cnf, 'qdel', is_critical=False)
        command = ' '.join(done_job_ids)
        if qdel:
            res = call(cnf, qdel + ' ' + command, exit_on_error=False)
            if res == 0:
                info('All running jobs for this project has been deleted from queue.')
            else:
                warn('Can\'t run qdel. Please kill the remaning jobs manually using the following command:')
                warn('  qdel ' + command)
        else:
            warn('Can\'t find qdel. Please kill the remaning jobs manually using the following command:')
            warn('  qdel ' + command)
        info()