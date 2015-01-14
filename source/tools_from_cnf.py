from genericpath import isfile, isdir, exists
import os
import sys
import subprocess
from os.path import join, dirname, abspath, pardir
from distutils.version import LooseVersion
from source.file_utils import verify_file, code_base_path, adjust_system_path, verify_dir, verify_obj_by_path

from source.logger import info, err
from source.file_utils import file_exists, which


def get_system_path(cnf, interpreter, name=None,
                    extra_warning='', suppress_warn=False):
    """ "name" can be:
        - key in system_into.yaml
        - relative path in the project (e.g. external/...)
        - anything in system path
    """
    if name is None:
        name = interpreter
        interpreter = None

    if interpreter:
        if interpreter == 'java':
            return get_java_tool_cmdline(cnf, name, extra_warning, suppress_warn)

        return get_script_cmdline(
            cnf, interpreter, name,
            extra_warning=extra_warning, suppress_warn=suppress_warn)

    # IN SYSTEM CONFIG?
    if (cnf.resources is not None and
        name.lower() in cnf.resources and
        'path' in cnf.resources[name.lower()]):

        tool_path = cnf.resources[name.lower()]['path']
        tool_path = adjust_system_path(tool_path)
        return verify_obj_by_path(tool_path, name)

    # IN PROJECT ROOT DIR? IN EXTERNAL?
    for dirpath in [code_base_path]:
        tool_path = join(dirpath, name)
        if exists(tool_path):
            return verify_obj_by_path(tool_path, name)

    # IN PATH?
    tool_path = which(name)
    if tool_path and exists(tool_path):
        return verify_obj_by_path(tool_path, name)

    if not suppress_warn:
        err(name + ' was not found. '
            'You may either specify path in the system config, '
            'or load into your PATH environment variable.')
    if extra_warning:
        err(extra_warning)
    return None


def system_path(*args, **kwargs):
    cmdline = get_system_path(*args, **kwargs)
    if not cmdline:
        sys.exit(1)
    return cmdline


def get_script_cmdline(cnf, interpreter, script,
                       interpreter_params='', extra_warning='', suppress_warn=False):
    interp_path = get_system_path(cnf, interpreter)
    if not interp_path:
        return None

    tool_path = get_system_path(cnf, script)
    if not tool_path:
        return None

    return interp_path + ' ' + interpreter_params + ' ' + tool_path


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