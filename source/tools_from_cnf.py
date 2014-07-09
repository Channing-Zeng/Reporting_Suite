import sys
import subprocess
from os.path import join, dirname, abspath
from distutils.version import LooseVersion
from source.file_utils import verify_file

from source.logger import info, err
from source.utils_from_bcbio import file_exists, which


def get_tool_cmdline(cnf, tool_name, interpreter=None,
                     extra_warning='', suppress_warn=False):
    if interpreter:
        return get_script_cmdline(
            cnf, tool_name, interpreter,
            extra_warning, suppress_warn)
    else:
        interpreter = ''

    # IN SYSTEM CONFIG?
    if (cnf.resources is not None and
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
