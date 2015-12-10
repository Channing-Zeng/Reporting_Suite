import subprocess
from genericpath import exists
from os.path import join
from distutils.version import LooseVersion
import sys

from source.file_utils import code_base_path, adjust_system_path, verify_obj_by_path
from source.logger import info, err, critical
from source.file_utils import which


def get_system_path(cnf, interpreter_or_name, name=None, extra_warning='',
                    suppress_warn=False, is_critical=False):
    """ "name" can be:
        - key in system_into.yaml
        - relative path in the project (e.g. external/...)
        - anything in system path
    """
    interpreter = interpreter_or_name
    if name is None:
        name = interpreter_or_name
        interpreter = None

    if interpreter:
        if interpreter == 'java':
            return get_java_tool_cmdline(cnf, name, extra_warning, suppress_warn, is_critical=is_critical)

        return get_script_cmdline(cnf, interpreter, name,
            extra_warning=extra_warning, suppress_warn=suppress_warn, is_critical=is_critical)

    # IN SYSTEM CONFIG?
    if cnf and (cnf.resources is not None and
        name.lower() in cnf.resources and
        'path' in cnf.resources[name.lower()]):

        tool_path = cnf.resources[name.lower()]['path']
        tool_path = adjust_system_path(tool_path)
        return verify_obj_by_path(tool_path, name, is_critical=is_critical)

    # IN PROJECT ROOT DIR? IN EXTERNAL?
    for dirpath in [code_base_path]:
        tool_path = join(dirpath, name)
        if exists(tool_path):
            return verify_obj_by_path(tool_path, name, is_critical=is_critical)

    # IN PATH?
    tool_path = which(name)
    if tool_path and exists(tool_path):
        return verify_obj_by_path(tool_path, name, is_critical=is_critical)

    msg = (name + ' was not found. You may either specify path in the system '
        'config, or load into your PATH environment variable. ' + extra_warning)
    if not suppress_warn:
        err(msg)
    if is_critical:
        critical(msg)
    return None


def system_path(*args, **kwargs):
    cmdline = get_system_path(*args, is_critical=True, **kwargs)
    return cmdline


def get_script_cmdline(cnf, interpreter, script, interpreter_params='',
                       extra_warning='', suppress_warn=False, is_critical=False):
    # if interpreter == 'python':
    #     interp_path = sys.executable
    # else:
    interp_path = get_system_path(cnf, interpreter, is_critical=is_critical)
    if not interp_path:
        return None

    tool_path = get_system_path(cnf, script, is_critical=is_critical)
    if not tool_path:
        return None

    return interp_path + ' ' + interpreter_params + ' ' + tool_path


def get_java_tool_cmdline(cnf, script, extra_warning='', suppress_warn=False, is_critical=False):
    jvm_opts = None
    if cnf and (cnf.resources and
        script in cnf.resources and
        'jvm_opts' in cnf.resources[script]):
        jvm_opts = cnf.resources[script]['jvm_opts']
    else:
        jvm_opts = ['-Xms750m', '-Xmx3g']
    if not '-Djava.io.tmpdir=' + join(cnf.work_dir) in jvm_opts:
        jvm_opts.append('-Djava.io.tmpdir=' + join(cnf.work_dir))

    return get_script_cmdline(
        cnf, 'java', script,
        interpreter_params=(' '.join(jvm_opts) + ' -jar'),
        extra_warning='', suppress_warn=False, is_critical=is_critical)


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


def get_snpeff_type(tool_cmdline):
    """New options format since v4.2
    """
    if LooseVersion(_get_snpeff_version(tool_cmdline)) >= LooseVersion("4.2"):
        return "new"
    else:
        return "old"


def _get_snpeff_version(tool_cmdline):
    cmdline = tool_cmdline + ' -version'

    with subprocess.Popen(cmdline,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT,
                          shell=True).stdout as stdout:
        version = stdout.read().strip()
    return version