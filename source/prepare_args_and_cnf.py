from distutils import file_util
import sys
import os
from os.path import join, pardir, isfile, isdir, expanduser, dirname, abspath, basename
from os import getcwd
from source import logger

from source.file_utils import verify_dir, safe_mkdir, adjust_path, verify_file, adjust_system_path, verify_obj_by_path, \
    remove_quotes, file_exists
from source.config import defaults, Config
from source.logger import info, critical, warn, err
from source.file_utils import which, file_exists, safe_mkdir
from source.ngscat.bed_file import verify_bam, verify_bed


def add_post_bcbio_args(parser):
    parser.add_option('--sys-cnf', '--sys-info', '--sys-cfg', dest='sys_cnf', help='System configuration yaml with paths to external tools and genome resources (see default one %s)' % defaults['sys_cnf'])
    parser.add_option('--run-cnf', '--run-info', '--run-cfg', dest='run_cnf', help='Run configuration yaml (see default one %s)' % defaults['run_cnf'])
    # parser.add_option('-v', dest='verbose', action='store_true', help='Verbose')
    parser.add_option('-t', dest='threads', type='int', help='Number of threads for each process', default=1)
    parser.add_option('--reuse', dest='reuse_intermediate', action='store_true', help='Reuse intermediate non-empty files in the work dir from previous run')
    # parser.add_option('--runner', dest='qsub_runner', help='Bash script that takes command line as the 1st argument. This script will be submitted to GRID. Default: ' + defaults['qsub_runner'])
    parser.add_option('--project-name', '--project', dest='project_name')
    # parser.add_option('--email', dest='email')
    parser.add_option('--bed', dest='bed', help='BED file to run targetSeq and Seq2C analysis on.')
    parser.add_option('--exons', '--exome', dest='exons', help='Exons BED file to make targetSeq exon/amplicon regions reports.')
    parser.add_option('--genome', dest='genome', help='Genome build.', default='hg19')
    parser.add_option('-f', '--freq', '--min-freq', dest='min_freq', type='float', help='Minimum allele frequency for the filtering. Default %f' % defaults['default_min_freq'])


def check_genome_resources(cnf):
    if not cnf.genomes:
        critical('"genomes" section is not specified in system config.')

    info('Checking paths in the genomes sections in ' + cnf.sys_cnf)
    info()

    for build_name, genome_cnf in cnf.genomes.items():
        info(build_name)
        for key in genome_cnf.keys():
            if isinstance(genome_cnf[key], basestring):
                genome_cnf[key] = adjust_system_path(genome_cnf[key])

            if not verify_obj_by_path(genome_cnf[key], key):
                if not genome_cnf[key].endswith('.gz') and verify_file(genome_cnf[key] + '.gz'):
                    gz_fpath = genome_cnf[key] + '.gz'
                    if verify_file(gz_fpath):
                        info(key + ': ' + gz_fpath)
                        genome_cnf[key] = gz_fpath
                else:
                    err('   err: no ' + genome_cnf[key] + (' and .gz' if not genome_cnf[key].endswith('gz') else ''))
            else:
                info(key + ': ' + genome_cnf[key])
        info()
        genome_cnf['name'] = build_name

    cnf.genome = cnf.genomes[cnf.genome]

    info('Checked genome resources.')
    info('*' * 70)
    info()


def check_keys(cnf, required_keys):
    to_exit = False

    for key in required_keys:
        if key not in cnf or not cnf[key]:
            to_exit = True
            err('Error: "' + key + '" must be provided in options or '
                'in ' + cnf.run_cnf + '.')
    return not to_exit


def check_inputs(cnf, file_keys=list(), dir_keys=list()):
    to_exit = False

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
                to_exit = True
            else:
                cnf[key] = adjust_path(cnf[key])

    for key in dir_keys:
        if key and key in cnf and cnf[key]:
            cnf[key] = adjust_system_path(cnf[key])
            if not verify_dir(cnf[key], key):
                to_exit = True
            else:
                cnf[key] = abspath(expanduser(cnf[key]))

    return not to_exit


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
    if not cnf.work_dir:
        work_dir_name = 'work' + ('_' + cnf.name if cnf.name else '')
        cnf.work_dir = join(cnf.output_dir, work_dir_name)
        # if not cnf.reuse_intermediate and isdir(cnf.work_dir):
        #     rmtree(cnf.work_dir)
    else:
        cnf.work_dir = adjust_path(cnf.work_dir)

    safe_mkdir(cnf.work_dir, 'working directory')


def set_up_log(cnf):
    logger.proc_name = cnf.proc_name
    logger.project_name = cnf.project_name
    logger.project_fpath = cnf.project_fpath or cnf.output_dir
    logger.address = remove_quotes(cnf.email) if cnf.email else ''
    logger.smtp_host = cnf.smtp_host

    log_fname = (cnf.proc_name + '_' if cnf.proc_name else '') + cnf.name + '_log.txt'
    log_fpath = join(cnf.log_dir, log_fname)

    if not cnf.proc_name:
        i = 1
        if file_exists(log_fpath):
            bak_fpath = log_fpath + '.' + str(i)
            while isfile(bak_fpath):
                bak_fpath = log_fpath + '.' + str(i)
                i += 1
            os.rename(log_fpath, bak_fpath)
        elif isfile(log_fpath):
            try:
                os.remove(log_fpath)
            except OSError:
                pass

    info('log_fpath: ' + log_fpath)
    logger.log_fpath = cnf.log = log_fpath


def determine_cnf_files(opts):
    opts.sys_cnf = adjust_path(opts.sys_cnf) if opts.sys_cnf else detect_sys_cnf(opts)
    if not verify_file(opts.sys_cnf): sys.exit(1)
    info('Using ' + opts.sys_cnf)

    opts.run_cnf = adjust_path(opts.run_cnf) if opts.run_cnf else defaults['run_cnf']
    if not verify_file(opts.run_cnf): sys.exit(1)
    info('Using ' + opts.run_cnf)


def detect_sys_cnf(opts):
    import socket
    hostname = socket.gethostname()
    info('hostname: ' + hostname)
    opts.sys_cnf = defaults['sys_cnfs']['us']
    if 'ukap' in hostname:
        opts.sys_cnf = defaults['sys_cnfs']['uk']
    elif 'cniclhpc' in hostname:
        opts.sys_cnf = defaults['sys_cnfs']['china']
    elif 'local' in hostname or 'Home' in hostname:
        opts.sys_cnf = defaults['sys_cnfs']['local']
    elif any(name in hostname for name in ['rask', 'blue', 'chara', 'usbod']):
        opts.sys_cnf = defaults['sys_cnfs']['us']
    return opts.sys_cnf
