import sys
import subprocess
import os
import shutil
from os.path import isfile, exists, join
from source.file_utils import make_tmpfile

from source.logger import info, err
from source.transaction import file_transaction
from source.utils_from_bcbio import file_exists


def call_pipe(cnf, cmdline, *args, **kwargs):
    return call_subprocess(cnf, cmdline, return_proc=True, *args, **kwargs)


def call_check_output(cnf, cmdline, *args, **kwargs):
    return call_subprocess(cnf, cmdline, check_output=True, *args, **kwargs)


def call(cnf, cmdline, output_fpath=None, *args, **kwargs):
    return call_subprocess(cnf, cmdline, None, output_fpath, *args, **kwargs)


def call_subprocess(cnf, cmdline, input_fpath_to_remove=None, output_fpath=None,
         stdout_to_outputfile=True, to_remove=list(), output_is_dir=False,
         stdin_fpath=None, exit_on_error=True, silent=False,

         overwrite=False, check_output=False, return_proc=False, print_stderr=True):
    """
    Required arguments:
    ------------------------------------------------------------
    cnf:                            dict with the following _optional_ fields:
                                      - reuse_intermediate
                                      - keep_intermediate
                                      - log
                                      - tmp_dir
    cmdline:                        called using subprocess.Popen
    ------------------------------------------------------------

    Optional arguments:
    ------------------------------------------------------------
    input_fpath_to_remove:          removed if not keep_intermediate
    output_fpath:                   overwritten if reuse_intermediate
    stdout_to_outputfile:           stdout=open(output_fpath, 'w')
    to_remove:                      list of files removed after the process finished
    output_is_dir                   output_fpath is a directory
    stdin_fpath:                    stdin=open(stdin_fpath)
    exit_on_error:                  is return code != 0, exit

    overwrite                       overwrite even if reuse_intermediate=True
    check_output                    subprocess.check_output; returns stdout
    return_proc                     proc = subprocess.Popen; returns proc

    ------------------------------------------------------------
    """

    if output_fpath is None:
        stdout_to_outputfile = False

    # NEEDED TO REUSE?
    if output_fpath and cnf.reuse_intermediate and not overwrite:
        if file_exists(output_fpath):
            info(output_fpath + ' exists, reusing')
            return output_fpath
    if output_fpath and file_exists(output_fpath):
        if output_is_dir:
            shutil.rmtree(output_fpath)
        else:
            os.remove(output_fpath)

    # ERR FILE TO STORE STDERR. IF SUBPROCESS FAIL, STDERR PRINTED
    err_fpath = join(cnf['work_dir'], '.subprocess_stderr.txt')
    if exists(err_fpath):
        os.remove(err_fpath)

    to_remove.append(err_fpath)

    def clean():
        for fpath in to_remove:
            if fpath and isfile(fpath):
                os.remove(fpath)

        if not cnf.keep_intermediate and input_fpath_to_remove:
            os.remove(input_fpath_to_remove)

    # RUN AND PRINT OUTPUT
    def do(cmdl, out_fpath=None):
        stdout = subprocess.PIPE
        stderr = subprocess.STDOUT

        if cnf['verbose'] or return_proc or check_output:
            if out_fpath:
                # STDOUT TO PIPE OR TO FILE
                if stdout_to_outputfile:
                    if not silent:
                        info(cmdl + ' > ' + out_fpath + (' < ' + stdin_fpath if stdin_fpath else ''))
                    stdout = open(out_fpath, 'w')
                    stderr = subprocess.PIPE
                else:
                    if output_fpath:
                        cmdl = cmdl.replace(output_fpath, out_fpath)
                    if not silent:
                        info(cmdl + (' < ' + stdin_fpath if stdin_fpath else ''))
                    stdout = subprocess.PIPE
                    stderr = subprocess.STDOUT
            else:
                if not silent:
                    info(cmdl + (' < ' + stdin_fpath if stdin_fpath else ''))

            if check_output:
                res = subprocess.check_output(
                    cmdl, shell=True, stderr=stderr, stdin=open(stdin_fpath) if stdin_fpath else None)
                clean()
                return res

            proc = subprocess.Popen(cmdl, shell=True, stdout=stdout, stderr=stderr,
                                    stdin=open(stdin_fpath) if stdin_fpath else None)
            if return_proc:
                # TODO: make this yield (as well as other returns),
                # TODO: ...move cleaning to finally, and use this function from with statement
                clean()
                return proc

            else:
                # PRINT STDOUT AND STDERR
                if proc.stdout:
                    for line in iter(proc.stdout.readline, ''):
                        info('   ' + line.strip())
                elif proc.stderr and print_stderr:
                    for line in iter(proc.stderr.readline, ''):
                        info('   ' + line.strip())

            # CHECK RES CODE
            ret_code = proc.wait()
            if ret_code != 0:
                for to_remove_fpath in to_remove:
                    if to_remove_fpath and isfile(to_remove_fpath):
                        os.remove(to_remove_fpath)
                err()
                err('Command returned status ' + str(ret_code) +
                    ('. Log saved to ' + cnf.log if cnf.log is not None else '.'))
                if exit_on_error:
                    clean()
                    sys.exit(1)
                else:
                    return None
            return output_fpath

        else:  # NOT VERBOSE, KEEP STDERR TO ERR FILE
            if out_fpath:
                # STDOUT TO PIPE OR TO FILE
                if stdout_to_outputfile:
                    if not silent:
                        info(cmdl + ' > ' + out_fpath + (' < ' + stdin_fpath if stdin_fpath else ''))
                    stdout = open(out_fpath, 'w')
                    stderr = open(err_fpath, 'a') if err_fpath else open('/dev/null')
                else:
                    if output_fpath:
                        cmdl = cmdl.replace(output_fpath, out_fpath)
                    if not silent:
                        info(cmdl + (' < ' + stdin_fpath if stdin_fpath else ''))
                    stdout = open(err_fpath, 'a') if err_fpath else open('/dev/null')
                    stderr = subprocess.STDOUT
            else:
                if not silent:
                    info(cmdl + (' < ' + stdin_fpath if stdin_fpath else ''))

            ret_code = subprocess.call(
                cmdl, shell=True, stdout=stdout, stderr=stderr,
                stdin=open(stdin_fpath) if stdin_fpath else None)

            # PRINT STDOUT AND STDERR
            if ret_code != 0:
                with open(err_fpath) as err_f:
                    info('')
                    info(err_f.read())
                    info('')
                for to_remove_fpath in to_remove:
                    if to_remove_fpath and isfile(to_remove_fpath):
                        os.remove(to_remove_fpath)
                err()
                err('Command returned status ' + str(ret_code) +
                    ('. Log saved to ' + cnf.log if 'log' in cnf else '.'))
                if exit_on_error:
                    clean()
                    sys.exit(1)
                else:
                    return None
            else:
                if cnf.log and err_fpath:
                    with open(err_fpath) as err_f, \
                         open(cnf.log, 'a') as log_f:
                        log_f.write('')
                        log_f.write(err_f.read())
                        log_f.write('')
            return output_fpath

    res = None
    if output_fpath and not output_is_dir:
        with file_transaction(cnf, output_fpath) as tx_out_fpath:
            res = do(cmdline, tx_out_fpath)
    else:
        res = do(cmdline)
        if res is not None:
            clean()
            return res

    clean()

    if res:
        if output_fpath and not output_is_dir:
            info('Saved to ' + output_fpath)

        return output_fpath
    else:
        return res
