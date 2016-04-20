import base64
import hashlib
import sys
import subprocess
import os
import shutil
from os.path import isfile, exists, join, islink
import time
import select

from source.logger import info, err, warn, critical, silent_err
from source.file_utils import file_exists, file_transaction, verify_file, verify_obj_by_path


def call_pipe(cnf, cmdline, *args, **kwargs):
    return call_subprocess(cnf, cmdline, return_proc=True, *args, **kwargs)


def call_check_output(cnf, cmdline, *args, **kwargs):
    return call_subprocess(cnf, cmdline, check_output=True, *args, **kwargs)


def call(cnf, cmdline, output_fpath=None, *args, **kwargs):
    return call_subprocess(cnf, cmdline, None, output_fpath, *args, **kwargs)


def call_subprocess(cnf, cmdline, input_fpath_to_remove=None, output_fpath=None,
         stdout_to_outputfile=True, to_remove=list(), output_is_dir=False,
         stdin_fpath=None, exit_on_error=True, silent=False,

         overwrite=False, check_output=False, verify_output_not_empty=True, return_proc=False, print_stderr=True,
         return_err_code=False, stderr_dump=None, max_number_of_tries=20,

         env_vars=None):
    """
    Required arguments:
    ------------------------------------------------------------
    cnf                             dict with the following _optional_ fields:
                                      - reuse_intermediate
                                      - keep_intermediate
                                      - verbose
                                      - log
                                      - work_dir
    cmdline                         called using subprocess.Popen
    ------------------------------------------------------------

    Optional arguments:
    ------------------------------------------------------------
    input_fpath_to_remove           removed if not keep_intermediate
    output_fpath                    overwritten if reuse_intermediate
    stdout_to_outputfile            stdout=open(output_fpath, 'w')
    to_remove                       list of files removed after the process finished
    output_is_dir                   output_fpath is a directory
    stdin_fpath                     stdin=open(stdin_fpath)
    exit_on_error                   is return code != 0, exit

    overwrite                       overwrite even if reuse_intermediate=True
    check_output                    subprocess.check_output; returns stdout
    return_proc                     proc = subprocess.Popen; returns proc
    return_err_code                 if return code !=0, return this code (only if exit_on_error=False)
    stderr_dump                     list. also store stderr here

    env_vars                        dictionary of environment variables to set only for this subprocess call
    ------------------------------------------------------------
    """

    if output_fpath is None:
        stdout_to_outputfile = False
    elif islink(output_fpath):
        os.unlink(output_fpath)

    # NEEDED TO REUSE?
    if output_fpath and cnf.reuse_intermediate and not overwrite:
        if verify_obj_by_path(output_fpath, silent=True):
            info(output_fpath + ' exists, reusing')
            return output_fpath
        if not output_fpath.endswith('.gz') and verify_obj_by_path(output_fpath + '.gz', silent=True):
            info(output_fpath + '.gz exists, reusing')
            return output_fpath
    if output_fpath and isfile(output_fpath):
        if output_is_dir:
            shutil.rmtree(output_fpath)
        else:
            os.remove(output_fpath)

    def clean():
        for fpath in to_remove:
            if fpath and isfile(fpath):
                os.remove(fpath)

        if not cnf.keep_intermediate and input_fpath_to_remove:
            os.remove(input_fpath_to_remove)

    env = os.environ.copy()
    if env_vars:

        for k, v in env_vars.items():
            if v is None:
                if k in env:
                    del env[k]
            else:
                env[k] = v

    # RUN AND PRINT OUTPUT
    def do(cmdl, out_fpath=None, stderr_dump=None):
        stdout = subprocess.PIPE    # goes to proc.stdout
        stderr = subprocess.STDOUT  # goes to proc.stderr

        if cnf.verbose or return_proc or check_output:
            if out_fpath:
                # STDOUT TO PIPE OR TO FILE
                if stdout_to_outputfile:
                    if not silent:
                        info(cmdl + (' < ' + stdin_fpath if stdin_fpath else '') + ' > ' + out_fpath)
                    stdout = open(out_fpath, 'w')
                    stderr = subprocess.PIPE
                else:
                    if output_fpath:
                        cmdl = cmdl.replace(output_fpath, out_fpath)
                    if not silent:
                        info(cmdl + (' < ' + stdin_fpath if stdin_fpath else ''))
                    stdout = subprocess.PIPE
                    stderr = subprocess.STDOUT  #TODO: fix interlacing err/out and change to PIPE here
            else:
                if not silent:
                    info(cmdl + (' < ' + stdin_fpath if stdin_fpath else ''))

            if check_output:
                res = subprocess.check_output(cmdl, shell=True, stderr=stderr,
                    stdin=open(stdin_fpath) if stdin_fpath else None, env=env,
                    executable='/bin/bash')
                clean()
                return res

            proc = subprocess.Popen(cmdl, shell=True, stdout=stdout, stderr=stderr,
                stdin=open(stdin_fpath) if stdin_fpath else None, env=env,
                executable='/bin/bash')

            if stderr_dump is None:
                stderr_dump = []

            if return_proc:
                # TODO: make this yield (as well as other returns),
                # TODO: ...move cleaning to finally, and use this function from with statement
                clean()
                return proc

            else:
                # streams = []
                # if proc.stderr and print_stderr:
                #     streams.append(proc.stderr.fileno())
                # if proc.stdout:
                #     streams.append(proc.stdout.fileno())
                #
                # ret = select.select(streams, [], [])
                # for fd in ret[0]:
                #     if proc.stdout and fd == proc.stdout.fileno():
                #         line = proc.stdout.readline().rstrip()
                #         info('  ' + line)
                #     if proc.stderr and fd == proc.stderr.fileno():
                #         line = proc.stderr.readline().rstrip()
                #         if stdout_to_outputfile:  # if output fpath is specified, then we usually write stdout to the file, and piping stderr (redirect to stdout)
                #             info('  ' + line)
                #         else:
                #             warn('  ' + line)
                #             stderr_dump.append(line)

                # PRINT STDOUT AND STDERR
                if proc.stdout:
                    for line in iter(proc.stdout.readline, ''):
                        if stderr == subprocess.PIPE:
                            silent_err('   ' + line.strip())
                        else:
                            info('   ' + line.strip())
                elif proc.stderr and print_stderr:
                    for line in iter(proc.stderr.readline, ''):
                        silent_err('   ' + line.strip())
                        stderr_dump.append(line)
                    info()

            # CHECK RES CODE
            ret_code = proc.wait()
            if ret_code != 0:
                for to_remove_fpath in to_remove:
                    if to_remove_fpath and isfile(to_remove_fpath):
                        os.remove(to_remove_fpath)
                warn()
                warn('Command returned status ' + str(ret_code) +
                    ('. Log saved to ' + cnf.log if cnf.log is not None else '.'))

                msg = 'Command returned status ' + str(ret_code)
                if cnf.log: msg += '. Log saved to ' + cnf.log
                msg += '\n'
                msg += '\n'
                if stderr_dump:
                    msg += 'Stderr:\n'
                    for l in stderr_dump:
                        msg += l
                # send_email(msg)

                if exit_on_error:
                    clean()
                    critical(msg)
                else:
                    if return_err_code:
                        return ret_code
                    else:
                        return None
            return output_fpath

        else:  # NOT VERBOSE, KEEP STDERR TO ERR FILE
            # ERR FILE TO STORE STDERR. IF SUBPROCESS FAIL, STDERR PRINTED
            hasher = hashlib.sha1(str(cmdline))
            cmdl_hash = base64.urlsafe_b64encode(hasher.digest()[0:4])[:-2]
            err_fpath_ = join(cnf.work_dir, cmdl_hash + '.subprocess_stderr.txt')
            to_remove.append(err_fpath_)
            if exists(err_fpath_):
                os.remove(err_fpath_)

            if out_fpath:
                # STDOUT TO TO FILE
                if stdout_to_outputfile:
                    if not silent:
                        info(cmdl + (' < ' + stdin_fpath if stdin_fpath else '') + ' > ' + out_fpath)
                    stdout = open(out_fpath, 'w')
                    if err_fpath_:
                        stderr = open(err_fpath_, 'a')
                        info('Writing err to ' + err_fpath_)
                    else:
                        info('Writing err to /dev/null')
                        stderr = open('/dev/null')
                else:
                # STDOUT TO PIPE
                    if output_fpath:
                        cmdl = cmdl.replace(output_fpath, out_fpath)
                    if not silent:
                        info(cmdl + (' < ' + stdin_fpath if stdin_fpath else ''))
                    stdout = subprocess.STDOUT
                    if err_fpath_:
                        stderr = open(err_fpath_, 'a')
                        info('Writing err to ' + err_fpath_)
                    else:
                        info('Writing err to /dev/null')
                        stderr = open('/dev/null')
            else:
                if not silent:
                    info(cmdl + (' < ' + stdin_fpath if stdin_fpath else ''))
                    stderr = None
                    stdout = None    # goes to proc.stdout

            ret_code = subprocess.call(
                cmdl, shell=True, stdout=stdout, stderr=stderr,
                stdin=open(stdin_fpath) if stdin_fpath else None, env=env,
                executable='/bin/bash')

            # PRINT STDOUT AND STDERR
            if ret_code != 0:
                if err_fpath_ and isfile(err_fpath_):
                    with open(err_fpath_) as err_f:
                        stderr_dump = err_f.read()

                    info('')
                    warn(stderr_dump)
                    info('')

                for to_remove_fpath in to_remove:
                    if to_remove_fpath and isfile(to_remove_fpath):
                        os.remove(to_remove_fpath)
                err()
                msg = ('Command returned status ' + str(ret_code) +
                       (('. Log saved to ' + cnf.log) if 'log' in cnf else '.'))
                if exit_on_error:
                    if cnf.log:
                        msg += '. Log saved to ' + cnf.log
                    msg += '\n'
                    msg += '\n'
                    if stderr_dump:
                        msg += 'Stderr:\n'
                        msg += stderr_dump
                if exit_on_error:
                    clean()
                    critical(msg)
                else:
                    err(msg)
                    return None
            else:
                if cnf.log and err_fpath_ and isfile(err_fpath_):
                    with open(err_fpath_) as err_f, \
                         open(cnf.log, 'a') as log_f:
                        log_f.write('')
                        log_f.write(err_f.read())
                        log_f.write('')
            return output_fpath

    def do_handle_oserror(cmdl, out_fpath=None, stderr_dump=None, max_number_of_tries=20):
        res_ = None
        counter = 0
        slept = 0
        timeout = 30
        limit = 60 * 10
        while True:
            try:
                res_ = do(cmdl, out_fpath, stderr_dump=stderr_dump)
                break
            except OSError, e:
                counter += 1
                if counter >= max_number_of_tries:
                    break
                err('OSError: ' + str(e))
                err()
                if 'Cannot allocate memory' not in str(e):
                    break
                else:
                    if slept >= limit:
                        return None
                    else:
                        err('Waiting ' + str(timeout) + ' seconds...')
                        time.sleep(timeout)
                        slept += timeout
                        err('Retrying...')
                        err()
        return res_

    res = None  # = proc or output_fpath
    if output_fpath and not output_is_dir:
        with file_transaction(cnf.work_dir, output_fpath) as tx_out_fpath:
            res = do_handle_oserror(cmdline, tx_out_fpath, stderr_dump=stderr_dump, max_number_of_tries=max_number_of_tries)
    else:
        res = do_handle_oserror(cmdline, stderr_dump=stderr_dump, max_number_of_tries=max_number_of_tries)
        if res is not None:
            clean()
            return res

    clean()

    if res:
        if output_fpath and not output_is_dir:
            info('Saved to ' + output_fpath)
            if not verify_file(output_fpath, verify_size=verify_output_not_empty):
                if exit_on_error:
                    clean()
                    critical('Error: the output file for the program is empty or does not exist; exiting.')
                else:
                    return None

        return output_fpath
    else:
        return res


# def process_out(p):
#     while True:
#         reads = [p.stdout.fileno(), p.stderr.fileno()]
#         ret = select.select(reads, [], [])
#
#         for fd in ret[0]:
#             if fd == p.stdout.fileno():
#                 read = p.stdout.readline()
#                 info(read)
#             if fd == p.stderr.fileno():
#                 read = p.stderr.readline()
#                 warn(read)
#
#         if p.poll() is not None:
#             break