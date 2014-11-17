from genericpath import isfile, isdir, exists
import os
import sys
import subprocess
from os.path import join, dirname, abspath, pardir
from distutils.version import LooseVersion
from source.file_utils import verify_file, code_base_path, adjust_system_path, verify_dir, verify_obj_by_path

from source.logger import info, err
from source.file_utils import file_exists, which


def tool_cmdline(*args, **kwargs):
    cmdline = get_system_path(*args, **kwargs)
    if not cmdline:
        exit(1)
    return cmdline


def get_system_path(cnf, interpreter, name_in_sys_cnf=None,
                    extra_warning='', suppress_warn=False):
    if name_in_sys_cnf is None:
        name_in_sys_cnf = interpreter
        interpreter = None

    if interpreter:
        if interpreter == 'java':
            return get_java_tool_cmdline(cnf, name_in_sys_cnf, extra_warning, suppress_warn)

        return get_script_cmdline(
            cnf, interpreter, name_in_sys_cnf,
            extra_warning=extra_warning, suppress_warn=suppress_warn)

    # IN SYSTEM CONFIG?
    if (cnf.resources is not None and
        name_in_sys_cnf.lower() in cnf.resources and
        'path' in cnf.resources[name_in_sys_cnf.lower()]):

        tool_path = cnf.resources[name_in_sys_cnf.lower()]['path']
        tool_path = adjust_system_path(tool_path)
        return verify_obj_by_path(tool_path, name_in_sys_cnf)

    # IN PROJECT ROOT DIR? IN EXTERNAL?
    for dirpath in [code_base_path, join(code_base_path, 'external')]:
        tool_path = join(dirpath, name_in_sys_cnf)
        if exists(tool_path):
            return verify_obj_by_path(tool_path, name_in_sys_cnf)

    # IN PATH?
    tool_path = which(name_in_sys_cnf)
    if tool_path and exists(tool_path):
        return verify_obj_by_path(tool_path, name_in_sys_cnf)


    if not suppress_warn:
        err(name_in_sys_cnf + ' was not found. '
            'You may either specify path in the system config, '
            'or load into your PATH environment variable.')
    if extra_warning:
        err(extra_warning)
    return None


def get_java_tool_cmdline(cnf, script, extra_warning='', suppress_warn=False):

    if (cnf.resources and
        script in cnf.resources and
        'jvm_opts' in cnf.resources[script]):
        jvm_opts = cnf.resources[script]['jvm_opts']
    else:
        jvm_opts = ['-Xms750m', '-Xmx3g']

    return get_script_cmdline(
        cnf, 'java', script,
        interpreter_params=(' '.join(jvm_opts) + ' -jar'),
        extra_warning='', suppress_warn=False)


def get_script_cmdline(cnf, interpreter, name_in_sys_cnf, script_fname=None,
                       interpreter_params='', extra_warning='', suppress_warn=False):
    interp_path = get_system_path(cnf, interpreter)
    if not interp_path:
        return None

    tool_path = get_system_path(cnf, name_in_sys_cnf)
    if not tool_path:
        return None

    if script_fname:
        tool_path = join(tool_path, script_fname)
        if not verify_file(tool_path):
            return None

    return interp_path + ' ' + interpreter_params + ' ' + tool_path


def get_gatk_cmdline(cnf):
    executable = get_java_tool_cmdline(cnf, 'gatk')
    if not executable:
        sys.exit(1)
    gatk_opts_line = ''
    if cnf.gatk:
        if 'options' in cnf.gatk:
            gatk_opts_line = ' '.join(cnf.gatk['options'])
    if 'threads' in cnf and ' -nt ' not in gatk_opts_line:
        gatk_opts_line += ' -nt ' + str(cnf.threads)
    return executable + ' ' + gatk_opts_line


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


def get_qualimap_type(tool_cmdline):
    """Qualimap supports multi-bamqc functionality only starting from v.2.0
    """
    if LooseVersion(_get_qualimap_version(tool_cmdline)) >= LooseVersion("2.0"):
        return "full"
    else:
        return "limited"


def _get_qualimap_version(tool_cmdline):
    cmdline = tool_cmdline + ' -version'  # actually, Qualimap doesn't have -version option

    version = None
    with subprocess.Popen(cmdline,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT,
                          shell=True).stdout as stdout:
        out = stdout.read().strip()
        flag = "QualiMap v."
        if out.startswith(flag) >= 0:
            version = out.split(flag)[-1].strip()
    if not version:
        info('WARNING: could not determine Qualimap version, using 1.0')
        return '1.0'
    if version.split('.') > 2:  # only major version
        version = '.'.join(version.split('.')[:2])
    return version