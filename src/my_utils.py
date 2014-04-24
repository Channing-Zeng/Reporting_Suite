import sys
import subprocess
import os
from os.path import join, basename, isfile, getsize, exists
from distutils.version import LooseVersion
from annotate import iterate_file
from src.utils import add_suffix


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


def get_gatk_type(tool_cmdline):
    """Retrieve type of GATK jar, allowing support for older GATK lite.
    Returns either `lite` (targeting GATK-lite 2.3.9) or `restricted`,
    the latest 2.4+ restricted version of GATK.
    """
    if LooseVersion(_gatk_major_version(tool_cmdline)) > LooseVersion("2.3"):
        return "restricted"
    else:
        return "lite"


def err(log, msg=None):
    if not msg:
        msg, log = log, None
    if log:
        open(log, 'a').write('\n' + msg + '\n')
    sys.stderr.write('\n' + msg + '\n')


def critical(log, msg=None):
    if not msg:
        msg, log = log, None
    if log:
        open(log, 'a').write('\n' + msg + '\n')
    exit(msg)


def info(log, msg=None):
    if not msg:
        msg, log = log, None
    print(msg)
    if log:
        open(log, 'a').write(msg + '\n')


def remove_quotes(s):
    if s and s[0] == '"':
        s = s[1:]
    if s and s[-1] == '"':
        s = s[:-1]
    return s


def verify_file(fpath, description=''):
    if not fpath:
        sys.stderr.write((description + ': f' if description else 'F') + 'ile name is empty.\n')
        return False
    if not exists(fpath):
        sys.stderr.write((description + ': ' if description else '') + fpath + ' does not exist.\n')
        return False
    if not isfile(fpath):
        sys.stderr.write((description + ': ' if description else '') + fpath + ' not a file.\n')
        return False
    if getsize(fpath) <= 0:
        sys.stderr.write((description + ': ' if description else '') + fpath + ' is empty.\n')
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
