import contextlib
import hashlib
import sys
import subprocess
import tempfile
import os
import shutil
import re
from os.path import join, basename, isfile, isdir, getsize, exists, expanduser, dirname, abspath
from distutils.version import LooseVersion

from source.logger import info, err, critical
from source.transaction import file_transaction
from source.utils_from_bcbio import add_suffix, file_exists, which, open_gzipsafe, safe_mkdir


def remove_quotes(s):
    if s and s[0] == '"':
        s = s[1:]
    if s and s[-1] == '"':
        s = s[:-1]
    return s


def _tryint(s):
    try:
        return int(s)
    except ValueError:
        return s


def _alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [_tryint(c) for c in re.split('([0-9]+)', s)]


def human_sorted(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=_alphanum_key)
    return l


def verify_module(name):
    try:
        __import__(name)
        return True
    except ImportError:
        return False


def verify_file(fpath, description=''):
    if not fpath:
        err((description + ': f' if description else 'F') + 'ile name is empty.')
        return None

    fpath = expanduser(fpath)
    if not exists(fpath):
        err((description + ': ' if description else '') + fpath + ' does not exist.')
        return None

    if not isfile(fpath):
        err((description + ': ' if description else '') + fpath + ' is not a file.')
        return None

    if getsize(fpath) <= 0:
        err((description + ': ' if description else '') + fpath + ' is empty.')
        return None

    return fpath


def verify_dir(fpath, description=''):
    if not fpath:
        err((description + ': d' if description else 'D') + 'ir name is empty.')
        return None

    fpath = expanduser(fpath)
    if not exists(fpath):
        err((description + ': ' if description else '') + fpath + ' does not exist.')
        return None

    if not isdir(fpath):
        err((description + ': ' if description else '') + fpath + ' is not a directory.')
        return None

    return fpath


def make_tmpdir(cnf, prefix='ngs_reporting_tmp', *args, **kwargs):
    base_dir = cnf.tmp_base_dir or cnf.work_dir or os.getcwd()
    if not verify_dir(base_dir, 'Base directory for temporary files'):
        sys.exit(1)

    return tempfile.mkdtemp(dir=base_dir, prefix=prefix)


@contextlib.contextmanager
def tmpdir(cnf, *args, **kwargs):
    prev_tmp_dir = cnf.tmp_dir

    cnf.tmp_dir = make_tmpdir(cnf, *args, **kwargs)
    try:
        yield cnf.tmp_dir
    finally:
        try:
            shutil.rmtree(cnf.tmp_dir)
        except OSError:
            pass
        cnf.tmp_dir = prev_tmp_dir


def make_tmpfile(cnf, *args, **kwargs):
    base_dir = cnf.tmp_base_dir or cnf.work_dir or os.getcwd()
    if not verify_dir(base_dir, 'Base directory for temporary files'):
        sys.exit(1)

    return tempfile.mkstemp(dir=base_dir, *args, **kwargs)


@contextlib.contextmanager
def tmpfile(cnf, *args, **kwargs):
    tmp_file, fpath = make_tmpfile(cnf, *args, **kwargs)
    try:
        yield fpath
    finally:
        try:
            os.remove(fpath)
        except OSError:
            pass


def iterate_file(cnf, input_fpath, proc_line_fun, suffix=None,
                 keep_original_if_not_keep_intermediate=False,
                 reuse_intermediate=True):
    output_fpath = intermediate_fname(cnf, input_fpath, suf=suffix or 'tmp')

    if suffix and cnf.reuse_intermediate and reuse_intermediate:
        if file_exists(output_fpath):
            info(output_fpath + ' exists, reusing')
            return output_fpath

    with file_transaction(cnf.tmp_dir, output_fpath) as tx_fpath:
        with open(input_fpath) as vcf, open(tx_fpath, 'w') as out:
            for i, line in enumerate(vcf):
                clean_line = line.strip()
                if clean_line:
                    new_l = proc_line_fun(clean_line)
                    if new_l is not None:
                        out.write(new_l + '\n')
                else:
                    out.write(line)

    if not suffix:
        os.rename(output_fpath, input_fpath)
        output_fpath = input_fpath
    else:
        if (not cnf.keep_intermediate and
            not keep_original_if_not_keep_intermediate and
                input_fpath):
            os.remove(input_fpath)
    return output_fpath


def index_bam(cnf, bam_fpath):
    samtools = get_tool_cmdline(cnf, 'samtools')
    if not samtools:
        sys.exit(1)

    cmdline = '{samtools} index {bam_fpath}'.format(**locals())
    call_subprocess(cnf, cmdline, None, None)


def bgzip_and_tabix_vcf(cnf, vcf_fpath):
    bgzip = get_tool_cmdline(cnf, 'bgzip', suppress_warn=True)
    tabix = get_tool_cmdline(cnf, 'tabix', suppress_warn=True)

    gzipped_fpath = join(cnf['work_dir'], basename(vcf_fpath) + '.gz')
    tbi_fpath = gzipped_fpath + '.tbi'

    if bgzip and not file_exists(gzipped_fpath):
        info('BGzipping VCF')
        cmdline = '{bgzip} -c {vcf_fpath}'.format(**locals())
        call_subprocess(cnf, cmdline, None, gzipped_fpath, exit_on_error=False)

    if tabix and not file_exists(tbi_fpath):
        info('Tabixing VCF')
        cmdline = '{tabix} -f -p vcf {gzipped_fpath}'.format(**locals())
        call_subprocess(cnf, cmdline, None, tbi_fpath, exit_on_error=False)

    return gzipped_fpath, tbi_fpath


def md5_for_file(f, block_size=2**20):
    md5 = hashlib.md5()

    while True:
        data = f.read(block_size)
        if not data:
            break
        md5.update(data)

    return md5.hexdigest()


def check_file_changed(cnf, new, in_work):
    if not file_exists(in_work):
        cnf['reuse_intermediate'] = False

    if cnf.get('reuse_intermediate'):
        if (basename(in_work) != basename(new) or
            md5_for_file(open(in_work, 'rb')) !=
            md5_for_file(open_gzipsafe(new, 'rb'))):

            info('Input file %s changed, setting "reuse_intermediate" '
                'to False.' % str(new))
            cnf['reuse_intermediate'] = False


# def check_inputs_changed(cnf, new_inputs):
#     prev_input_fpath = join(cnf['work_dir'], 'prev_inputs.txt')
#
#     new_inp_hashes = {realpath(fn): md5_for_file(fn) for fn in new_inputs}
#
#     if cnf.get('reuse_intermediate'):
#         if not file_exists(prev_input_fpath):
#             info('File %s does not exist, setting "reuse_intermediate" to '
#                       'False.' % str(prev_input_fpath))
#             cnf['reuse_intermediate'] = False
#
#         else:
#             prev_inp_hashes = dict()
#
#             with open(prev_input_fpath) as f:
#                 for l in f:
#                     input_fname, md5 = l.strip().split('\t')
#                     prev_inp_hashes[input_fname] = md5
#
#             if len(new_inp_hashes) != len(prev_inp_hashes):
#                 info('Number of input files changed, setting "reuse_intermediate" to False.')
#                 cnf['reuse_intermediate'] = False
#
#             for inp_fpath, inp_hash in new_inp_hashes.items():
#                 if inp_fpath not in prev_inp_hashes:
#                     info('Input changed, setting "reuse_intermediate" to False.')
#                     cnf['reuse_intermediate'] = False
#
#                 if inp_hash != prev_inp_hashes[inp_fpath]:
#                     info('Input %s changed, setting "reuse_intermediate" '
#                               'to False.' % str(inp_fpath))
#                     cnf['reuse_intermediate'] = False
#
#     with open(prev_input_fpath, 'w') as f:
#         for inp_fpath, inp_hash in new_inp_hashes.items():
#             f.write(inp_fpath + '\t' + inp_hash + '\n')


def get_tool_cmdline(cnf, tool_name, interpreter='',
                     extra_warning='', suppress_warn=False):
    if interpreter:
        return get_script_cmdline(
            cnf, tool_name, interpreter,
            extra_warning, suppress_warn)

    # IN SYSTEM CONFIG?
    if (cnf.resources and
        tool_name.lower() in cnf.resources and
        'path' in cnf.resources[tool_name.lower()]):

        tool_path = cnf.resources[tool_name.lower()]['path']
        return verify_file(tool_path, tool_name)

    # IN PROJECT ROOT DIR?
    project_base_path = dirname(dirname(abspath(__file__)))
    tool_path = join(project_base_path, tool_name)
    if file_exists(tool_path):
        return verify_file(tool_path, tool_name)

    # SOME OF BASIC SCRIPTS?
    tool_path += '.py'
    if file_exists(tool_path):
        return verify_file(tool_path, tool_name)

    # IN PATH?
    tool_path = which(tool_name)
    if file_exists(tool_path):
        return verify_file(tool_path, tool_name)

    if not suppress_warn:
        err(tool_name + ' was not found. '
            'You may either specify path in the system config, '
            'or load into your PATH environment variable.')
    if extra_warning:
        err(extra_warning)
    return None


def get_java_tool_cmdline(cnf, script,
         extra_warning='', suppress_warn=False):

    if (cnf.resources and
        script in cnf.resources and
        'jvm_opts' in cnf.resources[script]):
        jvm_opts = cnf.resources[script]['jvm_opts']
    else:
        jvm_opts = ['']

    return get_script_cmdline(cnf, 'java', script,
                              (' '.join(jvm_opts) + ' -jar'),
                              extra_warning='', suppress_warn=False)


def get_script_cmdline(cnf, interpreter, script, interpreter_params='',
                       extra_warning='', suppress_warn=False):
    interp_path = get_tool_cmdline(cnf, interpreter)
    if not interp_path:
        return None

    tool_path = get_tool_cmdline(cnf, script)
    if not tool_path:
        return None

    return interp_path + ' ' + interpreter_params + ' ' + tool_path


def call_pipe(cnf, cmdline, *args, **kwargs):
    return call_subprocess(cnf, cmdline, return_proc=True, *args, **kwargs)


def call_check_output(cnf, cmdline, *args, **kwargs):
    return call_subprocess(cnf, cmdline, check_output=True, *args, **kwargs)


def call(cnf, cmdline, output_fpath=None, *args, **kwargs):
    return call_subprocess(cnf, cmdline, None, output_fpath, *args, **kwargs)


def call_subprocess(cnf, cmdline, input_fpath_to_remove=None, output_fpath=None,
         stdout_to_outputfile=True, to_remove=list(), output_is_dir=False,
         stdin_fpath=None, exit_on_error=True, silent=False,

         overwrite=True, check_output=False, return_proc=False):
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
    if output_fpath and cnf.reuse_intermediate and overwrite is False:
        if file_exists(output_fpath):
            info(output_fpath + ' exists, reusing')
            return output_fpath
    if output_fpath and file_exists(output_fpath):
        if output_is_dir:
            shutil.rmtree(output_fpath)
        else:
            os.remove(output_fpath)

    # ERR FILE TO STORE STDERR. IF SUBPROCESS FAIL, STDERR PRINTED
    err_fpath = None
    _, err_fpath = make_tmpfile(cnf, prefix='err_tmp')
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

            if check_output:
                res = subprocess.check_output(
                    cmdl, shell=True, stderr=stderr, stdin=open(stdin_fpath) if stdin_fpath else None)
                clean()
                return res

            proc = subprocess.Popen(cmdl, shell=True, stdout=stdout, stderr=stderr,
                                    stdin=open(stdin_fpath) if stdin_fpath else None)
            if return_proc:
                # TODO: make this yield (as well as other returns), move cleaning to finally, and use this function from with statement
                clean()
                return proc

            else:
                # PRINT STDOUT AND STDERR
                if proc.stdout:
                    for line in iter(proc.stdout.readline, ''):
                        info('   ' + line.strip())
                elif proc.stderr:
                    for line in iter(proc.stderr.readline, ''):
                        info('   ' + line.strip())

            # CHECK RES CODE
            ret_code = proc.wait()
            if ret_code != 0:
                for to_remove_fpath in to_remove:
                    if to_remove_fpath and isfile(to_remove_fpath):
                        os.remove(to_remove_fpath)
                err('Command returned status ' + str(ret_code) +
                    ('. Log in ' + cnf['log'] if 'log' in cnf else '.'))
                if exit_on_error:
                    clean()
                    sys.exit(1)

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
                err('Command returned status ' + str(ret_code) +
                    ('. Log in ' + cnf.log if 'log' in cnf else '.'))
                if exit_on_error:
                    clean()
                    sys.exit(1)
            else:
                if cnf.log and err_fpath:
                    with open(err_fpath) as err_f, \
                         open(cnf.log, 'a') as log_f:
                        log_f.write('')
                        log_f.write(err_f.read())
                        log_f.write('')

    if output_fpath and not output_is_dir:
        with file_transaction(cnf.tmp_dir, output_fpath) as tx_out_fpath:
            do(cmdline, tx_out_fpath)
    else:
        res = do(cmdline)
        if res is not None:
            clean()
            return res

    clean()

    if output_fpath and not output_is_dir:
        info('Saved to ' + output_fpath)

    return output_fpath


def join_parent_conf(child_conf, parent_conf):
    bc = parent_conf.copy()
    bc.update(child_conf)
    child_conf.update(bc)
    return child_conf


def intermediate_fname(cnf, fname, suf):
    output_fname = add_suffix(fname, suf)
    return join(cnf['work_dir'], basename(output_fname))


def dots_to_empty_cells(config, tsv_fpath):
    """Put dots instead of empty cells in order to view TSV with column -t
    """
    def proc_line(l):
        while '\t\t' in l:
            l = l.replace('\t\t', '\t.\t')
        return l
    return iterate_file(config, tsv_fpath, proc_line, 'dots')


def get_gatk_type(tool_cmdline):
    """Retrieve type of GATK jar, allowing support for older GATK lite.
    Returns either `lite` (targeting GATK-lite 2.3.9) or `restricted`,
    the latest 2.4+ restricted version of GATK.
    """
    if LooseVersion(_gatk_major_version(tool_cmdline)) > LooseVersion("2.3"):
        return "restricted"
    else:
        return "lite"


def _gatk_major_version(config):
    """Retrieve the GATK major version, handling multiple GATK distributions.

    Has special cases for GATK nightly builds, Appistry releases and
    GATK prior to 2.3.
    """
    full_version = _get_gatk_version(config)
    # Working with a recent version if using nightlies
    if full_version.startswith("nightly-"):
        return "2.8"
    parts = full_version.split("-")
    if len(parts) == 4:
        appistry_release, version, subversion, githash = parts
    elif len(parts) == 3:
        version, subversion, githash = parts
    # version was not properly implemented in earlier GATKs
    else:
        version = "2.3"
    if version.startswith("v"):
        version = version[1:]
    return version


def _get_gatk_version(tool_cmdline):
    cmdline = tool_cmdline + ' -version'

    version = None
    with subprocess.Popen(cmdline,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT,
                          shell=True).stdout as stdout:
        out = stdout.read().strip()
        last_line = out.split('\n')[-1].strip()
        # versions earlier than 2.4 do not have explicit version command,
        # parse from error output from GATK
        if out.find("ERROR") >= 0:
            flag = "The Genome Analysis Toolkit (GATK)"
            for line in last_line.split("\n"):
                if line.startswith(flag):
                    version = line.split(flag)[-1].split(",")[0].strip()
        else:
            version = last_line
    if not version:
        info('WARNING: could not determine Gatk version, using 1.0')
        return '1.0'
    if version.startswith("v"):
        version = version[1:]
    return version


def format_integer(name, value, unit=''):
    value = int(value)
    if value is not None:
        return '{name}: {value:,}{unit}'.format(**locals())
    else:
        return '{name}: -'.format(**locals())


def format_decimal(name, value, unit=''):
    if value is not None:
        return '{name}: {value:.2f}{unit}'.format(**locals())
    else:
        return '{name}: -'.format(**locals())


def mean(ints):
    return float(sum(ints)) / len(ints) if len(ints) > 0 else float('nan')

