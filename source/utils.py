import sys
import subprocess
import tempfile
import os
import shutil
import re
from os.path import join, basename, isfile, isdir, getsize, exists, expanduser
from distutils.version import LooseVersion
from datetime import datetime

from source.transaction import file_transaction
from source.bcbio_utils import add_suffix, file_exists, which


def err(log, msg=None):
    if msg is None:
        msg, log = log, None
    if log:
        open(log, 'a').write('\n' + msg + '\n')
    sys.stderr.write('\n' + msg + '\n')
    sys.stderr.flush()


def critical(log, msg=None):
    if msg is None:
        msg, log = log, None
    if log:
        open(log, 'a').write('\n' + msg + '\n')
    exit(msg)


def info(log, msg=None):
    if msg is None:
        msg, log = log, None
    print(msg)
    sys.stdout.flush()
    if log:
        open(log, 'a').write(msg + '\n')


def remove_quotes(s):
    if s and s[0] == '"':
        s = s[1:]
    if s and s[-1] == '"':
        s = s[:-1]
    return s


def _tryint(s):
    try:
        return int(s)
    except:
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


def verify_file(fpath, description=''):
    if not fpath:
        sys.stderr.write((description + ': f' if description else 'F') + 'ile name is empty.\n')
        return False
    fpath = expanduser(fpath)
    if not exists(fpath):
        sys.stderr.write((description + ': ' if description else '') + fpath + ' does not exist.\n')
        return False
    if not isfile(fpath):
        sys.stderr.write((description + ': ' if description else '') + fpath + ' is not a file.\n')
        return False
    if getsize(fpath) <= 0:
        sys.stderr.write((description + ': ' if description else '') + fpath + ' is empty.\n')
        return False
    return True


def verify_dir(fpath, description=''):
    if not fpath:
        sys.stderr.write((description + ': d' if description else 'D') + 'ir name is empty.\n')
        return False
    fpath = expanduser(fpath)
    if not exists(fpath):
        sys.stderr.write((description + ': ' if description else '') + fpath + ' does not exist.\n')
        return False
    if not isdir(fpath):
        sys.stderr.write((description + ': ' if description else '') + fpath + ' is not a directory.\n')
        return False
    return True


def safe_mkdir(dirpath, descriptive_name):
    if not dirpath:
        exit(descriptive_name + ' path is empty.')
    if isfile(dirpath):
        exit(descriptive_name + ' ' + dirpath + ' is a file.')
    if not exists(dirpath):
        try:
            os.mkdir(dirpath)
        except OSError:
            exit('Parent directory for ' + descriptive_name +
                 ' ' + dirpath + ' probably does not exist.')


def iterate_file(cnf, input_fpath, proc_line_fun, work_dir, suffix=None,
                 keep_original_if_not_keep_intermediate=False):
    output_fpath = intermediate_fname(work_dir, input_fpath, suf=suffix or 'tmp')

    if suffix and cnf.get('reuse_intermediate'):
        if file_exists(output_fpath):
            info(cnf['log'], output_fpath + ' exists, reusing')
            return output_fpath

    with file_transaction(output_fpath) as tx_fpath:
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
        if (not cnf.get('keep_intermediate') and
            not keep_original_if_not_keep_intermediate and
                input_fpath):
            os.remove(input_fpath)
    return output_fpath


def get_tool_cmdline(sys_cnf, tool_name, extra_warn=''):
    tool_path = which(tool_name) or None

    if not 'resources' in sys_cnf \
            or tool_name not in sys_cnf['resources'] \
            or 'path' not in sys_cnf['resources'][tool_name]:
        if tool_path:
            return tool_path
        else:
            err(tool_name + ' executable was not found. '
                'You can either specify path in the system config, or load into your '
                'PATH environment variable.')
            if extra_warn:
                err(extra_warn)
            return None

    tool_path = sys_cnf['resources'][tool_name]['path']
    if verify_file(tool_path, tool_name):
        return tool_path
    else:
        err(tool_path + ' for ' + tool_name + ' does not exist or is not a file.')
        return None


def call(cnf, cmdline, input_fpath_to_remove, output_fpath,
         stdout_to_outputfile=True, to_remove=None, output_is_file=True,
         stdin=None):
    to_remove = to_remove or []

    # MAYBE REUSE?
    if output_fpath and cnf.get('reuse_intermediate'):
        if file_exists(output_fpath):
            info(cnf.get('log'), output_fpath + ' exists, reusing')
            return output_fpath
    if output_fpath and file_exists(output_fpath):
        if output_is_file:
            os.remove(output_fpath)
        else:
            shutil.rmtree(output_fpath)

    # ERR FILE TO STORE STDERR. IF SUBPROCESS FAIL, STDERR PRINTED
    err_fpath = None
    if cnf.get('work_dir'):
        _, err_fpath = tempfile.mkstemp(dir=cnf.get('work_dir'), prefix='err_tmp')
        to_remove.append(err_fpath)

    # RUN AND PRINT OUTPUT
    def do(cmdline, tx_out_fpath=None):
        stdout = subprocess.PIPE
        stderr = subprocess.STDOUT

        if cnf['verbose']:
            if tx_out_fpath:
                # STDOUT TO PIPE OR TO FILE
                if stdout_to_outputfile:
                    info(cnf.get('log'), cmdline + ' > ' + tx_out_fpath + (' < ' + stdin if stdin else ''))
                    stdout = open(tx_out_fpath, 'w')
                    stderr = subprocess.PIPE
                else:
                    cmdline = cmdline.replace(output_fpath, tx_out_fpath)
                    info(cnf.get('log'), cmdline + (' < ' + stdin if stdin else ''))
                    stdout = subprocess.PIPE
                    stderr = subprocess.STDOUT

            proc = subprocess.Popen(cmdline, shell=True, stdout=stdout, stderr=stderr,
                                    stdin=open(stdin) if stdin else None)

            # PRINT STDOUT AND STDERR
            if proc.stdout:
                for line in iter(proc.stdout.readline, ''):
                    info(cnf.get('log'), '   ' + line.strip())
            elif proc.stderr:
                for line in iter(proc.stderr.readline, ''):
                    info(cnf.get('log'), '   ' + line.strip())

            # CHECK RES CODE
            ret_code = proc.wait()
            if ret_code != 0:
                for fpath in to_remove:
                    if fpath and isfile(fpath):
                        os.remove(fpath)
                critical(cnf.get('log'), 'Command returned status ' + str(ret_code) +
                         ('. Log in ' + cnf['log'] if 'log' in cnf else '.'))

        else:  # NOT VERBOSE, KEEP STDERR TO ERR FILE
            if tx_out_fpath:
                # STDOUT TO PIPE OR TO FILE
                if stdout_to_outputfile:
                    info(cnf.get('log'), cmdline + ' > ' + tx_out_fpath + (' < ' + stdin if stdin else ''))
                    stdout = open(tx_out_fpath, 'w')
                    stderr = open(err_fpath, 'a') if err_fpath else open('/dev/null')
                else:
                    cmdline = cmdline.replace(output_fpath, tx_out_fpath)
                    info(cnf.get('log'), cmdline + (' < ' + stdin if stdin else ''))
                    stdout = open(err_fpath, 'a') if err_fpath else open('/dev/null')
                    stderr = subprocess.STDOUT

            res = subprocess.call(cmdline, shell=True, stdout=stdout, stderr=stderr,
                                  stdin=open(stdin) if stdin else None)

            # PRINT STDOUT AND STDERR
            if res != 0:
                with open(err_fpath) as err_f:
                    info(cnf.get('log'), '')
                    info(cnf.get('log'), err_f.read())
                    info(cnf.get('log'), '')
                for fpath in to_remove:
                    if fpath and isfile(fpath):
                        os.remove(fpath)
                critical(cnf.get('log'), 'Command returned status ' + str(res) +
                         ('. Log in ' + cnf['log'] if 'log' in cnf else '.'))
            else:
                if cnf.get('log') and err_fpath:
                    with open(err_fpath) as err_f, \
                         open(cnf.get('log'), 'a') as log_f:
                        log_f.write('')
                        log_f.write(err_f.read())
                        log_f.write('')

    if output_fpath and output_is_file:
        with file_transaction(output_fpath) as tx_out_fpath:
            do(cmdline, tx_out_fpath)
    else:
        do(cmdline)

    # REMOVE UNNESESSARY
    for fpath in to_remove:
        if fpath and isfile(fpath):
            os.remove(fpath)

    if not cnf.get('keep_intermediate') and input_fpath_to_remove:
        os.remove(input_fpath_to_remove)

    if output_fpath and output_is_file:
        info(cnf.get('log'), 'Saved to ' + output_fpath)
    return output_fpath


def get_java_tool_cmdline(cnf, name):
    cmdline_template = get_script_cmdline_template(cnf, 'java', name)
    jvm_opts = cnf['resources'][name].get('jvm_opts', []) + ['']
    return cmdline_template % (' '.join(jvm_opts) + ' -jar')


def get_script_cmdline_template(cnf, executable, script_name):
    if not which(executable):
        exit(executable + ' executable required, maybe you need '
             'to run "module load ' + executable + '"?')
    if 'resources' not in cnf:
        critical(cnf['log'], 'System config yaml must contain resources section with '
                 + script_name + ' path.')
    if script_name not in cnf['resources']:
        critical(cnf['log'], 'System config resources section must contain '
                 + script_name + ' info (with a path to the tool).')
    tool_config = cnf['resources'][script_name]
    if 'path' not in tool_config:
        critical(script_name + ' section in the system config must contain a path to the tool.')
    tool_path = tool_config['path']
    if not verify_file(tool_path, script_name):
        exit(1)
    return executable + ' %s ' + tool_path


def join_parent_conf(child_conf, parent_conf):
    bc = parent_conf.copy()
    bc.update(child_conf)
    child_conf.update(bc)
    return child_conf


def timestamp():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S  ")


def step_greetings(cnf, name=None):
    if name is None:
        name = cnf
        cnf = dict()
    if name is None:
        name = ''
    if cnf is None:
        cnf = dict()

    info(cnf.get('log'), '')
    info(cnf.get('log'), '-' * 70)
    info(cnf.get('log'), timestamp() + name)
    info(cnf.get('log'), '-' * 70)


def intermediate_fname(work_dir, fname, suf):
    output_fname = add_suffix(fname, suf)
    return join(work_dir, basename(output_fname))


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

