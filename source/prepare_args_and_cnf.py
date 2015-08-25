from distutils import file_util
import getpass
import hashlib
import sys
import os
from os.path import join, pardir, isfile, isdir, expanduser, dirname, abspath, basename
from os import getcwd
import datetime
from source import logger

from source.file_utils import verify_dir, safe_mkdir, adjust_path, verify_file, adjust_system_path, verify_obj_by_path, \
    remove_quotes, file_exists
from source.config import defaults, Config
from source.logger import info, critical, warn, err
from source.file_utils import which, file_exists, safe_mkdir
from source.ngscat.bed_file import verify_bam, verify_bed
from source.utils import is_uk, is_us, is_cloud
from source.utils import is_china
from source.utils import is_local


def add_cnf_t_reuse_prjname_reuse_marker_genome(parser):
    parser.add_option('--sys-cnf', '--sys-info', '--sys-cfg', dest='sys_cnf', help='System configuration yaml with paths to external tools and genome resources (see default one %s)' % defaults['sys_cnf'])
    parser.add_option('--run-cnf', '--run-info', '--run-cfg', dest='run_cnf', help='Run configuration yaml (see default one %s)' % defaults['run_cnf_exome_seq'])
    # parser.add_option('-v', dest='verbose', action='store_true', help='Verbose')
    parser.add_option('-t', dest='threads', type='int', help='Max number of threads', default=10)
    parser.add_option('--reuse', dest='reuse_intermediate', action='store_true', help='Reuse intermediate non-empty files in the work dir from previous run')
    # parser.add_option('--runner', dest='qsub_runner', help='Bash script that takes command line as the 1st argument. This script will be submitted to GRID. Default: ' + defaults['qsub_runner'])
    parser.add_option('--project-name', '--project', dest='project_name')
    parser.add_option('--done-marker', dest='done_marker')
    parser.add_option('--work-dir', dest='work_dir', metavar='DIR', help='Default is temporary directory')
    # parser.add_option('--email', dest='email')
    parser.add_option('--genome', dest='genome', help='Genome build.', default='hg19')
    parser.add_option('--debug', dest='debug', help='Debug mode; keep work directory.', action='store_true', default=False)


def check_genome_resources(cnf):
    if not cnf.genomes:
        critical('"genomes" section is not specgetified in system config.')

    info('Checking paths in the genomes sections in ' + cnf.sys_cnf)
    info()

    info('Genome: ' + str(cnf.genome))

    for build_name, genome_cnf in cnf.genomes.items():
        for key in genome_cnf.keys():
            if isinstance(genome_cnf[key], basestring):
                genome_cnf[key] = adjust_system_path(genome_cnf[key])

                if not verify_obj_by_path(genome_cnf[key], key, silent=True):
                    if not genome_cnf[key].endswith('.gz') and verify_file(genome_cnf[key] + '.gz', silent=True):
                        gz_fpath = genome_cnf[key] + '.gz'
                        if verify_file(gz_fpath, silent=True):
                            genome_cnf[key] = gz_fpath
                            if build_name == cnf.genome:
                                info(key + ': ' + gz_fpath)
                    else:
                        if build_name == cnf.genome:
                            err('   Err: no ' + genome_cnf[key] + (' and .gz' if not genome_cnf[key].endswith('gz') else ''))
                        else:
                            warn('   Warn: no ' + genome_cnf[key] + (' and .gz' if not genome_cnf[key].endswith('gz') else ''))
                # else:
                #     info(key + ': ' + genome_cnf[key])
        info()
        genome_cnf['name'] = build_name

    cnf.genome = cnf.genomes[cnf.genome]

    info('Checked genome resources.')
    info('*' * 70)
    info()


def check_keys(cnf, required_keys):
    errors = []

    for key in required_keys:
        if key not in cnf or not cnf[key]:
            to_exit = True
            errors.append('Error: "' + key + '" must be provided in options or '
                'in ' + cnf.run_cnf + '.')
    return errors


def check_inputs(cnf, file_keys=list(), dir_keys=list()):
    errors = []

    def _verify_input_file(_key):
        cnf[_key] = adjust_path(cnf[_key])
        if not verify_file(cnf[_key], _key):
            return False
        if 'bam' in _key and not verify_bam(cnf[_key]):
            return False
        if 'bed' in _key and not verify_bed(cnf[_key]):
            return False
        return True

    for key in file_keys:
        if key and key in cnf and cnf[key]:
            if not _verify_input_file(key):
                errors.append('File ' + cnf[key] + ' is empty or cannot be found')
            else:
                cnf[key] = adjust_path(cnf[key])

    for key in dir_keys:
        if key and key in cnf and cnf[key]:
            cnf[key] = adjust_system_path(cnf[key])
            if not verify_dir(cnf[key], key):
                errors.append('Directory ' + cnf[key] + ' is empty or cannot be found')
            else:
                cnf[key] = abspath(expanduser(cnf[key]))

    return errors


def input_fpaths_from_cnf(cnf, required_inputs, optional_inputs):
    input_fpaths = []
    for key in required_inputs:
        input_fpaths.append(cnf[key])
    for key in optional_inputs:
        if key in cnf:
            input_fpaths.append(cnf[key])
    return input_fpaths


def check_system_resources(cnf, required=list(), optional=list()):
    to_exit = False

    for program in required:
        if not which(program):
            if cnf.resources is None:
                critical('No "resources" section in system config.')

            data = cnf.resources.get(program)
            if data is None:
                err(program + ' is required. Specify path in system config or in your environment.')
                to_exit = True
            else:
                if 'module' in data:
                    os.system('module load ' + data['module'])
                    # if 'path' not in data:
                    #     data['path'] = program
                elif 'path' in data:
                    data['path'] = adjust_system_path(data['path'])
                    if not isdir(data['path']) and not file_exists(data['path']):
                        err(data['path'] + ' does not exist.')
                        to_exit = True

    for program in optional:
        resources = cnf.get('resources')
        if not resources:
            break

        data = resources.get(program)
        if data is None:
            continue
        else:
            data['path'] = adjust_system_path(data['path'])
            if not isdir(data['path']) and not file_exists(data['path']):
                err(data['path'] + ' does not exist.')
                to_exit = True

    if to_exit:
        exit()


def set_up_dirs(cnf):
    """ Creates output_dir, work_dir; sets up log
    """
    cnf.output_dir = adjust_path(cnf.output_dir)
    safe_mkdir(cnf.output_dir, 'output_dir')
    info('Saving into ' + cnf.output_dir)

    set_up_work_dir(cnf)

    if not cnf.log_dir:
        cnf.log_dir = join(cnf.work_dir, 'log')
    safe_mkdir(cnf.log_dir)

    set_up_log(cnf)


def set_up_work_dir(cnf):
    # timestamp = str(datetime.datetime.now())
    # user_prid = getpass.getuser()
    # hasher = hashlib.sha1( + timestamp)
    # path_hash = base64.urlsafe_b64encode(hasher.digest()[0:4])[:-1]

    if not cnf.work_dir:
        work_dir_name = 'work' + ('_' + cnf.sample if cnf.sample else '')
        cnf.work_dir = join(cnf.output_dir, work_dir_name)
        # if not cnf.reuse_intermediate and isdir(cnf.work_dir):
        #     rmtree(cnf.work_dir)
    else:
        cnf.work_dir = adjust_path(cnf.work_dir)

    safe_mkdir(cnf.work_dir, 'working directory')


def set_up_log(cnf, proc_name=None, project_name=None, project_fpath=None, output_dir=None):
    logger.proc_name = proc_name
    logger.project_name = project_name
    logger.project_fpath = project_fpath or output_dir
    logger.address = remove_quotes(cnf.email) if cnf.email else ''
    logger.smtp_host = cnf.smtp_host

    log_fname = (proc_name + '_' if proc_name else '') + (cnf.sample + '_' if cnf.sample else '') + '_log.txt'
    log_fpath = join(cnf.log_dir, log_fname)

    if not proc_name:
        i = 1
        if file_exists(log_fpath):
            bak_fpath = log_fpath + '.' + str(i)
            while isfile(bak_fpath):
                bak_fpath = log_fpath + '.' + str(i)
                i += 1
            try:
                os.rename(log_fpath, bak_fpath)
            except OSError:
                pass
        elif isfile(log_fpath):
            try:
                os.remove(log_fpath)
            except OSError:
                pass

    info('log_fpath: ' + log_fpath)
    logger.log_fpath = cnf.log = log_fpath


def determine_sys_cnf(opts):
    if opts.sys_cnf:
        return verify_file(opts.sys_cnf, is_critical=True)
    else:
        opts.sys_cnf = verify_file(detect_sys_cnf_by_location(), is_critical=True)

    info('Using ' + opts.sys_cnf)
    return opts.sys_cnf


def determine_run_cnf(opts, is_wgs=False, is_targeteq=False):
    if opts.run_cnf:
        opts.run_cnf = adjust_path(opts.run_cnf)
    elif is_wgs:
        opts.run_cnf = defaults['run_cnf_wgs']
    elif is_targeteq:
        opts.run_cnf = defaults['run_cnf_deep_seq']
    else:
        opts.run_cnf = defaults['run_cnf_exome_seq']

    verify_file(opts.run_cnf, is_critical=True)
    info('Using ' + opts.run_cnf)
    return opts.run_cnf


def detect_sys_cnf_by_location():
    sys_cnf = defaults['sys_cnfs']['us']
    if is_uk():
        sys_cnf = defaults['sys_cnfs']['uk']
    elif is_china():
        sys_cnf = defaults['sys_cnfs']['china']
    elif is_local():
        sys_cnf = defaults['sys_cnfs']['local']
    elif is_us():
        sys_cnf = defaults['sys_cnfs']['us']
    elif is_cloud():
        sys_cnf = defaults['sys_cnfs']['cloud']
    return sys_cnf
