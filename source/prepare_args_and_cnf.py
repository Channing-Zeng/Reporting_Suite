from optparse import OptionParser
from distutils import file_util
import sys
import os
from os.path import join, pardir, isfile, isdir, expanduser, dirname, abspath, basename
from os import getcwd

from source.bcbio_structure import BCBioStructure, load_bcbio_cnf
from source.file_utils import verify_dir, safe_mkdir, adjust_path, verify_file, adjust_system_path, verify_obj_by_path
from source.config import defaults, Config
from source.logger import info, critical, warn, err
from source.ngscat.bed_file import verify_bam, verify_bed


def add_post_bcbio_args(parser):
    parser.add_option('--sys-cnf', '--sys-info', '--sys-cfg', dest='sys_cnf', help='System configuration yaml with paths to external tools and genome resources (see default one %s)' % defaults['sys_cnf'])
    parser.add_option('--run-cnf', '--run-info', '--run-cfg', dest='run_cnf', help='Run configuration yaml (see default one %s)' % defaults['run_cnf'])
    # parser.add_option('-v', dest='verbose', action='store_true', help='Verbose')
    parser.add_option('-t', dest='threads', type='int', help='Number of threads for each process')
    parser.add_option('--reuse', dest='reuse_intermediate', action='store_true', help='Reuse intermediate non-empty files in the work dir from previous run')
    # parser.add_option('--runner', dest='qsub_runner', help='Bash script that takes command line as the 1st argument. This script will be submitted to GRID. Default: ' + defaults['qsub_runner'])
    parser.add_option('--project-name', '--project', dest='project_name')
    # parser.add_option('--email', dest='email')
    parser.add_option('--bed', dest='bed', help='BED file to run targetSeq and Seq2C analysis on.')
    parser.add_option('--exons', '--exome', dest='exons', help='Exons BED file to make targetSeq exon/amplicon regions reports.')
    parser.add_option('--genome', dest='genome', help='Genome build.', default='hg19')
    parser.add_option('-f', '--freq', '--min-freq', dest='min_freq', type='float', help='Minimum allele frequency for the filtering. Default %f' % defaults['default_min_freq'])


def process_post_bcbio_args(parser):
    (opts, args) = parser.parse_args()

    dir_arg = args[0] if len(args) > 0 else getcwd()
    dir_arg = adjust_path(dir_arg)
    if not verify_dir(dir_arg):
        sys.exit(1)

    bcbio_project_dirpath, final_dirpath, config_dirpath = _set_bcbio_dirpath(dir_arg)

    _set_sys_config(config_dirpath, opts)

    _set_run_config(config_dirpath, opts)

    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)

    if 'qsub_runner' in cnf:
        cnf.qsub_runner = adjust_system_path(cnf.qsub_runner)
    if not check_inputs(cnf, file_keys=['qsub_runner']):
        sys.exit(1)

    bcbio_cnf = load_bcbio_cnf(cnf, config_dirpath)

    check_genome_resources(cnf)

    return cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath


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


def summary_script_proc_params(name, dir_name=None, description=None, extra_opts=None):
    description = description or 'This script generates project-level summaries based on per-sample ' + name + ' reports.'
    parser = OptionParser(description=description)
    add_post_bcbio_args(parser)

    parser.add_option('--log-dir', dest='log_dir')
    parser.add_option('--dir', dest='dir_name', default=dir_name, help='Optional - to distinguish VarQC_summary and VarQC_after_summary')
    parser.add_option('--name', dest='name', default=name, help='Procedure name')
    for args, kwargs in extra_opts or []:
        parser.add_option(*args, **kwargs)

    cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath = process_post_bcbio_args(parser)

    bcbio_structure = BCBioStructure(cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath, cnf.name)

    cnf.output_dir = join(bcbio_structure.date_dirpath, cnf.dir_name) if cnf.dir_name else None
    cnf.work_dir = bcbio_structure.work_dir
    set_up_work_dir(cnf)

    info('*' * 70)
    info()

    return cnf, bcbio_structure


def determine_cnf_files(opts):
    opts.sys_cnf = adjust_path(opts.sys_cnf) if opts.sys_cnf else detect_sys_cnf(opts)
    if not verify_file(opts.sys_cnf): sys.exit(1)
    info('Using ' + opts.sys_cnf)

    opts.run_cnf = adjust_path(opts.run_cnf) if opts.run_cnf else defaults['run_cnf']
    if not verify_file(opts.run_cnf): sys.exit(1)
    info('Using ' + opts.run_cnf)


def check_keys(cnf, required_keys):
    to_exit = False

    for key in required_keys:
        if key not in cnf or not cnf[key]:
            to_exit = True
            err('Error: "' + key + '" must be provided in options or '
                'in ' + cnf.run_cnf + '.')
    return not to_exit


def set_up_work_dir(cnf):
    if not cnf.work_dir:
        work_dir_name = 'work_' + cnf.name
        cnf.work_dir = join(cnf.output_dir, work_dir_name)
        # if not cnf.reuse_intermediate and isdir(cnf.work_dir):
        #     rmtree(cnf.work_dir)
    else:
        cnf.work_dir = adjust_path(cnf.work_dir)

    safe_mkdir(cnf.work_dir, 'working directory')


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


def _set_bcbio_dirpath(dir_arg):
    final_dirpath, bcbio_project_dirpath, config_dirpath = None, None, None

    if isdir(join(dir_arg, 'config')):
        bcbio_project_dirpath = dir_arg
        final_dirpath = None
        config_dirpath = join(bcbio_project_dirpath, 'config')

    elif isdir(abspath(join(dir_arg, pardir, 'config'))):
        bcbio_project_dirpath = abspath(join(dir_arg, pardir))
        final_dirpath = dir_arg
        config_dirpath = join(bcbio_project_dirpath, 'config')

    else:
        critical(
            'No config directory ' + join(dir_arg, 'config') + ' or ' + abspath(join(dir_arg, pardir, 'config')) +
            ', check if you provided a correct path to the bcbio directory.'
            '\nIt can be provided by the first argument for the script, or by changing to it.')

    info('BCBio project directory: ' + bcbio_project_dirpath)
    if final_dirpath: info('final directory: ' + final_dirpath)
    info('Config directory: ' + config_dirpath)
    return bcbio_project_dirpath, final_dirpath, config_dirpath


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


def _set_sys_config(config_dirpath, opts):
    provided_cnf_fpath = adjust_path(opts.sys_cnf)
    # provided in commandline?
    if provided_cnf_fpath:
        if not verify_file(provided_cnf_fpath):
            sys.exit(1)
        # alright, in commandline.
        opts.sys_cnf = provided_cnf_fpath

    else:
        # or probably in config dir?
        fpaths_in_config = [
            abspath(join(config_dirpath, fname))
            for fname in os.listdir(config_dirpath)
            if fname.startswith('system_info') and fname.endswith('.yaml')]

        if len(fpaths_in_config) > 1:
            critical('More than one YAML file containing run_info in name found in the config '
                     'directory ' + config_dirpath + ': ' + ' '.join(fpaths_in_config))

        elif len(fpaths_in_config) == 1:
            opts.sys_cnf = fpaths_in_config[0]
            if not verify_file(opts.sys_cnf):
                sys.exit(1)
            # alright, in config dir.

        else:
            # detect system and use default system config
            detect_sys_cnf(opts)

    info('Using ' + opts.sys_cnf)


def _set_run_config(config_dirpath, opts):
    provided_cnf_fpath = adjust_path(opts.run_cnf)
    # provided in commandline?
    if provided_cnf_fpath:
        if not verify_file(provided_cnf_fpath):
            sys.exit(1)

        # alright, in commandline. copying over to config dir.
        opts.run_cnf = provided_cnf_fpath
        project_run_cnf_fpath = adjust_path(join(config_dirpath, basename(provided_cnf_fpath)))
        if provided_cnf_fpath != project_run_cnf_fpath:
            info('Using ' + provided_cnf_fpath + ', copying to ' + project_run_cnf_fpath)
            if isfile(project_run_cnf_fpath):
                try:
                    os.remove(project_run_cnf_fpath)
                except OSError:
                    pass

            if not isfile(project_run_cnf_fpath):
                run_info_fpaths_in_config = [
                    abspath(join(config_dirpath, fname))
                    for fname in os.listdir(config_dirpath)
                    if fname.startswith('run_info') and fname.endswith('.yaml')]

                if len(run_info_fpaths_in_config) > 0:
                    warn('Warning: there are run_info files in config directory ' + config_dirpath + '. '
                         'Provided config will be copied there and can cause ambigity in future.')

                file_util.copy_file(provided_cnf_fpath, project_run_cnf_fpath, preserve_times=False)

    else:  # no configs provided in command line options
        run_info_fpaths_in_config = [
            abspath(join(config_dirpath, fname))
            for fname in os.listdir(config_dirpath)
            if fname.startswith('run_info') and fname.endswith('.yaml')]

        if len(run_info_fpaths_in_config) > 1:
            critical('More than one YAML file containing run_info in name found in the config '
                     'directory ' + config_dirpath + ': ' + ' '.join(run_info_fpaths_in_config))

        elif len(run_info_fpaths_in_config) == 1:
            opts.run_cnf = run_info_fpaths_in_config[0]
            if not verify_file(opts.run_cnf):
                sys.exit(1)
            # alright, in config dir.

        elif len(run_info_fpaths_in_config) == 0:
            info('No YAMLs containing run_info in name found in the config directory ' +
                 config_dirpath + ', using the default one.')

            # using default one.
            opts.run_cnf = defaults['run_cnf']
            project_run_cnf_fpath = adjust_path(join(config_dirpath, basename(opts.run_cnf)))
            info('Using ' + opts.run_cnf + ', copying to ' + project_run_cnf_fpath)
            if isfile(project_run_cnf_fpath):
                try:
                    os.remove(project_run_cnf_fpath)
                except OSError:
                    pass
            if not isfile(project_run_cnf_fpath):
                file_util.copy_file(opts.run_cnf, project_run_cnf_fpath, preserve_times=False)

    info('Using ' + opts.run_cnf)


