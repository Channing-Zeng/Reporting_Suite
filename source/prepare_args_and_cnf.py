import socket
from distutils import file_util
import getpass
import hashlib
import sys
import os
from optparse import SUPPRESS_HELP
from os.path import join, pardir, isfile, isdir, expanduser, dirname, abspath, basename
from os import getcwd
import datetime
import tempfile

from az.reference_data import get_canonical_transcripts

import bcbio_postproc
from source import logger
from source.file_utils import verify_dir, safe_mkdir, adjust_path, verify_file, adjust_system_path, verify_obj_by_path, \
    remove_quotes, file_exists
from source.config import defaults, Config
from source.logger import info, critical, warn, err, debug
from source.file_utils import which, file_exists, safe_mkdir
from source.targetcov.bam_and_bed_utils import verify_bam, verify_bed
from source.utils import is_uk, is_us, is_cloud, is_sweden, is_ace, is_china, is_local, is_chihua


def add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser, threads=defaults['bcbio_postproc_threads']):
    parser.add_option('--sys-cnf', '--sys-info', '--sys-cfg', dest='sys_cnf',
                      help='System configuration yaml with paths to external tools and genome resources (see default one %s)' % defaults['sys_cnf'])
    parser.add_option('--run-cnf', '--run-info', '--run-cfg', dest='run_cnf',
                      help='Run configuration yaml (see default one %s)' % defaults['run_cnf_exome_seq'])
    parser.add_option('-t', dest='threads', type='int', default=threads, help='Max number of slots, default is %d' % threads)
    parser.add_option('--reuse', dest='reuse_intermediate', action='store_true', help='Reuse intermediate non-empty files in the work dir from previous run')
    parser.add_option('--no-check', dest='no_check', action='store_true', help=SUPPRESS_HELP)
    parser.add_option('--project-name', '--project', dest='project_name', help='Project name. If not set, it gets parsed from JIRA or from the location path.')
    parser.add_option('--genome', '-g', dest='genome', help='Genome build')
    parser.add_option('--done-marker', dest='done_marker', help=SUPPRESS_HELP)
    parser.add_option('--work-dir', dest='work_dir', metavar='DIR', help=SUPPRESS_HELP)  # Default is temporary directory
    parser.add_option('--debug', dest='debug', help=SUPPRESS_HELP, action='store_true', default=False)
    parser.add_option('--queue', dest='queue', help=SUPPRESS_HELP)  # Queue name for qsub
    parser.add_option('--priority', dest='qsub_priority', help='Priority to submit jobs to qsub')  # Queue name for qsub


def check_genome_resources(cnf):
    if cnf.genome is None:
        critical('Please, specify genome build (one of available in ' + cnf.sys_cnf +
                 ') using the --genome option (e.g., --genome hg38).')

    if not cnf.genomes:
        critical('"genomes" section is not specified in system config ' + cnf.sys_cnf)

    info('Genome: ' + str(cnf.genome.name))

    for key in cnf.genome.keys():
        if key != 'name' and isinstance(cnf.genome[key], basestring):
            cnf.genome[key] = adjust_system_path(cnf.genome[key])

            if not verify_obj_by_path(cnf.genome[key], key, silent=True):
                if not cnf.genome[key].endswith('.gz') and verify_file(cnf.genome[key] + '.gz', silent=True):
                    gz_fpath = cnf.genome[key] + '.gz'
                    if verify_file(gz_fpath, silent=True):
                        cnf.genome[key] = gz_fpath

    if not cnf.genome.features or not cnf.genome.bed_annotation_features or not cnf.genome.cds:
        warn('Warning: features and bed_annotation_features and cds in the system config (' + cnf.sys_cnf + ') must be specified.')

    if not cnf.transcripts_fpath:
        cnf.transcripts_fpath = cnf.transcripts_fpath or get_canonical_transcripts(cnf.genome.name, ensembl=True)
        # g = 'hg19' if (('hg19' or 'GRCh37') in cnf.genome.name) else 'hg38'
        # cnf.transcripts_fpath = verify_file(join(bcbio_postproc.project_dir,
        #     'reference_data', 'canonical_transcripts', 'canonical_transcripts_' + g + '.txt'), description='Canonical transcripts')


def check_keys_presence(cnf, required_keys):
    errors = []

    for key in required_keys:
        if key not in cnf or not cnf[key]:
            to_exit = True
            errors.append('Error: "' + key + '" must be provided in options or '
                'in ' + cnf.run_cnf + '.')
    return errors


def check_dirs_and_files(cnf, file_keys=list(), dir_keys=list()):
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
            cnf[key] = adjust_path(cnf[key])
            if not verify_dir(cnf[key], key):
                errors.append('Directory ' + cnf[key] + ' is empty or cannot be found')
            else:
                cnf[key] = adjust_path(cnf[key])

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


def set_up_dirs(cnf, log_dir_name='log'):
    """ Creates output_dir, work_dir; sets up log
    """
    if cnf.output_dir:
        cnf.output_dir = adjust_path(cnf.output_dir)
        safe_mkdir(cnf.output_dir, 'output_dir')
        info('Saving into ' + cnf.output_dir)

    set_up_work_dir(cnf)

    if cnf.log_dir == '-':
        cnf.log_dir = None
    else:
        if not cnf.log_dir:
            cnf.log_dir = join(cnf.work_dir, log_dir_name)
        safe_mkdir(cnf.log_dir)
        info('Created log dir ' + cnf.log_dir)

    set_up_log(cnf)


def set_up_work_dir(cnf):
    # timestamp = str(datetime.datetime.now())
    # user_prid = getpass.getuser()
    # hasher = hashlib.sha1( + timestamp)
    # path_hash = base64.urlsafe_b64encode(hasher.digest()[0:4])[:-1]

    if not cnf.work_dir:
        if cnf.output_dir:
            work_dir_name = 'work' + ('_' + cnf.sample if cnf.sample else '')
            cnf.work_dir = join(cnf.output_dir, work_dir_name)
            info('Work dir: ' + cnf.work_dir)
            # if not cnf.reuse_intermediate and isdir(cnf.work_dir):
            #     rmtree(cnf.work_dir)
        else:
            cnf.work_dir = tempfile.mkdtemp()
            info('Creating temprorary directory for work dir: ' + cnf.work_dir)
    else:
        cnf.work_dir = adjust_path(cnf.work_dir)
        info('Work dir: ' + cnf.work_dir)

    safe_mkdir(cnf.work_dir, 'working directory')


def set_up_log(cnf, proc_name=None, project_name=None, project_fpath=None, output_dir=None):
    logger.proc_name = proc_name
    logger.project_name = project_name
    logger.project_fpath = project_fpath or output_dir
    logger.cnf_address = remove_quotes(cnf.email) if cnf.email else ''
    logger.smtp_host = cnf.smtp_host

    if cnf.log_dir:
        log_fname = (proc_name + '_' if proc_name else '') + (cnf.sample + '_' if cnf.sample else '') + 'log.txt'
        log_fpath = join(cnf.log_dir, log_fname)

        if file_exists(log_fpath):
            timestamp = datetime.datetime.fromtimestamp(os.stat(log_fpath).st_mtime)
            mv_log_fpath = log_fpath + '.' + timestamp.strftime("%Y-%m-%d_%H-%M-%S")
            try:
                if isfile(mv_log_fpath):
                    os.remove(mv_log_fpath)
                if not isfile(mv_log_fpath):
                    os.rename(log_fpath, mv_log_fpath)
            except OSError:
                pass
        info('log_fpath: ' + log_fpath)
        info()
        logger.log_fpath = cnf.log = log_fpath


def determine_sys_cnf(opts):
    if 'sys_cnf' in opts.__dict__ and opts.sys_cnf:
        return verify_file(opts.sys_cnf, is_critical=True)
    else:
        opts.__dict__['sys_cnf'] = verify_file(detect_sys_cnf_by_location(), is_critical=True)

    debug('Using system configuration ' + opts.sys_cnf)
    return opts.sys_cnf


def determine_run_cnf(opts, is_wgs=False, is_targetseq=False):
    if opts.run_cnf:
        opts.run_cnf = adjust_path(opts.run_cnf)
    elif is_wgs:
        opts.run_cnf = defaults['run_cnf_wgs']
    elif is_targetseq:
        opts.run_cnf = defaults['run_cnf_deep_seq']
    else:
        opts.run_cnf = defaults['run_cnf_exome_seq']

    verify_file(opts.run_cnf, is_critical=True)
    debug('Using run configuration ' + opts.run_cnf)
    return opts.run_cnf


def detect_sys_cnf_by_location():
    if is_uk():
        res = defaults['sys_cnfs']['uk']
    elif is_sweden():
        res = defaults['sys_cnfs']['sweden']
    elif is_china():
        res = defaults['sys_cnfs']['china']
    elif is_us():
        res = defaults['sys_cnfs']['us']
    elif is_cloud():
        res = defaults['sys_cnfs']['cloud']
    elif is_local():
        res = defaults['sys_cnfs']['local']
    elif is_ace():
        res = defaults['sys_cnfs']['ace']
    elif is_chihua():
        res = defaults['sys_cnfs']['chihua']
    else:
        warn('Warning: could not detect location by hostname: ' + socket.gethostname() + '. Using local')
        res = defaults['sys_cnfs']['local']
    return res
