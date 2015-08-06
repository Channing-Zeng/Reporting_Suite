import os
import shutil
import sys
import re
from dircache import listdir
from genericpath import isdir, isfile
from collections import defaultdict, OrderedDict
from os.path import join, abspath, exists, pardir, splitext, basename, islink, dirname, realpath
from optparse import OptionParser
from distutils import file_util

import source
from source import logger, BaseSample
from source.logger import info, err, critical, warn
from source.calling_process import call
from source.config import load_yaml_config, Config, defaults
from source.file_utils import verify_dir, verify_file, adjust_path, remove_quotes, adjust_system_path
from source.ngscat.bed_file import verify_bed, verify_bam
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_reuse_marker_genome, set_up_log, set_up_work_dir, \
    detect_sys_cnf_by_location, check_genome_resources, check_inputs
from source.tools_from_cnf import get_system_path
from source.file_utils import file_exists, safe_mkdir
from source.utils import OrderedDefaultDict


def summary_script_proc_params(proc_name, proc_dir_name=None, description=None, extra_opts=None):
    description = description or 'This script generates project-level summaries based on per-sample ' + proc_name + ' reports.'
    parser = OptionParser(description=description)
    add_cnf_t_reuse_prjname_reuse_marker_genome(parser)

    parser.add_option('--log-dir', dest='log_dir')
    parser.add_option('--varqc-dir', dest='proc_dir_name', default=proc_dir_name, help='Optional - to distinguish VarQC_summary and VarQC_after_summary')
    parser.add_option('--varqc-name', dest='proc_name', default=proc_name, help='Procedure name')
    parser.add_option('-o', dest='output_dir', metavar='DIR')
    parser.add_option('--email', dest='email')

    for args, kwargs in extra_opts or []:
        parser.add_option(*args, **kwargs)

    cnf, bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths, tags = process_post_bcbio_args(parser)

    cnf_project_name = cnf.project_name
    if len(bcbio_project_dirpaths) > 1:
        cnf.project_name = None

    bcbio_structures = []
    for bcbio_project_dirpath, bcbio_cnf, final_dirpath in zip(bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths):
        bcbio_structures.append(BCBioStructure(cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath, cnf.proc_name))

    # Single project, running as usually
    if len(bcbio_structures) == 1:
        bcbio_structure = bcbio_structures[0]
        cnf.output_dir = join(bcbio_structure.date_dirpath, cnf.proc_dir_name) if cnf.proc_dir_name else None
        cnf.work_dir = cnf.work_dir or join(bcbio_structure.work_dir, cnf.proc_name)
        set_up_work_dir(cnf)

        info('*' * 70)
        info()

        return cnf, bcbio_structure

    # Special case: multiple projects input. Need to join summary reports instead.
    else:
        if cnf_project_name:
            cnf.project_name = cnf_project_name
        else:
            cnf.project_name == '_'.join([bs.project_name for bs in bcbio_structures])

        if cnf.output_dir is None:
            cnf.output_dir = join(os.getcwd(), cnf.project_name)

        return cnf, bcbio_structures


def process_post_bcbio_args(parser):
    (opts, args) = parser.parse_args()

    bcbio_dirpaths = []
    tags = []

    for dir_arg in args or [os.getcwd()]:
        # /ngs/oncology/Analysis/bioscience/Bio_0038_KudosCellLinesExomes/Bio_0038_150521_D00443_0159_AHK2KTADXX/bcbio,Kudos159 /ngs/oncology/Analysis/bioscience/Bio_0038_KudosCellLinesExomes/Bio_0038_150521_D00443_0160_BHKWMNADXX/bcbio,Kudos160
        dirpath = adjust_path(dir_arg.split(',')[0])
        bcbio_dirpaths.append(dirpath)
        if len(dir_arg.split(',')) > 1:
            tags.append(dir_arg.split(',')[1])
        else:
            tags.append(None)

    for dirpath, tag in zip(bcbio_dirpaths, tags):
        verify_dir(dirpath, is_critical=True, description='Path to bcbio project' + (' ' if tag else ''))

    cnf = None
    bcbio_project_dirpaths = []
    bcbio_cnfs = []
    final_dirpaths = []

    for dirpath in bcbio_dirpaths:
        bcbio_project_dirpath, final_dirpath, config_dirpath = _detect_bcbio_dirpath(dirpath)
        bcbio_project_dirpaths.append(bcbio_project_dirpath)
        final_dirpaths.append(final_dirpath)

        if cnf is None:
            _detect_sys_config(config_dirpath, opts)
            _detect_move_run_config(config_dirpath, opts)

            cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)
            check_genome_resources(cnf)

            if 'qsub_runner' in cnf:
                cnf.qsub_runner = adjust_system_path(cnf.qsub_runner)
            errors = check_inputs(cnf, file_keys=['qsub_runner'])
            if errors:
                critical(errors)

        bcbio_cnfs.append(load_bcbio_cnf(cnf, config_dirpath))

    return cnf, bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths, tags


def _detect_bcbio_dirpath(dir_arg):
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


def _detect_sys_config(config_dirpath, opts):
    provided_cnf_fpath = adjust_path(opts.sys_cnf)
    # provided in commandline?
    if provided_cnf_fpath:
        verify_file(provided_cnf_fpath, is_critical=True)
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
            verify_file(opts.sys_cnf, is_critical=True)
            # alright, in config dir.

        else:
            # detect system and use default system config
            opts.sys_cnf = detect_sys_cnf_by_location()

    info('Using ' + opts.sys_cnf)


def _detect_move_run_config(config_dirpath, opts):
    provided_cnf_fpath = adjust_path(opts.run_cnf)

    if not provided_cnf_fpath and opts.deep_seq:  # TODO: if there is a run_info in dir, ignore the deep_seq opt
        provided_cnf_fpath = defaults['run_cnf_deep_seq']

    # provided in commandline?
    if provided_cnf_fpath:
        info(provided_cnf_fpath + ' was provided in the command line options')

        verify_file(provided_cnf_fpath, is_critical=True)

        # alright, in commandline. copying over to config dir.
        opts.run_cnf = provided_cnf_fpath
        project_run_cnf_fpath = adjust_path(join(config_dirpath, basename(provided_cnf_fpath)))

        if not isfile(project_run_cnf_fpath) or not os.path.samefile(provided_cnf_fpath, project_run_cnf_fpath):
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
                    if (fname.startswith('run_info') or 'post' in fname) and fname.endswith('.yaml')]

                if len(run_info_fpaths_in_config) > 0:
                    warn('Warning: there are run_info files in config directory ' + config_dirpath + '. '
                         'Provided config will be copied there and can cause ambigity in future.')

                file_util.copy_file(provided_cnf_fpath, project_run_cnf_fpath, preserve_times=False)

    else:
        info('no configs provided in command line options')
        run_info_fpaths_in_config = [
            abspath(join(config_dirpath, fname))
            for fname in os.listdir(config_dirpath)
            if (fname.startswith('run_info') or 'post' in fname) and fname.endswith('.yaml')]

        if len(run_info_fpaths_in_config) > 1:
            critical('More than one YAML file containing run_info in name found in the config '
                     'directory ' + config_dirpath + ': ' + ' '.join(run_info_fpaths_in_config))

        elif len(run_info_fpaths_in_config) == 1:
            opts.run_cnf = run_info_fpaths_in_config[0]
            verify_file(opts.run_cnf, is_critical=True)
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


class BCBioSample(BaseSample):
    def __init__(self, name, final_dir, **kwargs):
        self.dirpath = join(final_dir, name)
        BaseSample.__init__(self, name, self.dirpath, '{dirpath}/{name}/', **kwargs)
        self.sv_bed = None
        self.ngscat_html_fpath             = self.make_fpath('{dirpath}/qc/{name}/captureQC.html', name=source.ngscat_name)
        self.qualimap_html_fpath           = self.make_fpath('{dirpath}/qc/{name}/qualimapReport.html', name=source.qualimap_name)
        self.qualimap_genome_results_fpath = self.make_fpath('{dirpath}/qc/{name}/genome_results.txt', name=source.qualimap_name)
        self.qualimap_ins_size_hist_fpath  = self.make_fpath('{dirpath}/qc/{name}/raw_data_qualimapReport/insert_size_histogram.txt', name=source.qualimap_name)
        self.fastqc_html_fpath             = self.make_fpath('{dirpath}/qc/{name}/fastqc_report.html', name=source.fastqc_name)
        self.project_tag = None

    # ----------
    def annotated_vcfs_dirpath(self):
        return join(self.dirpath, BCBioStructure.varannotate_dir)

    def get_filtered_vcfs_dirpath(self):
        return join(self.dirpath, BCBioStructure.varfilter_dir)

    # raw variants
    def find_raw_vcf_by_callername(self, callername):
        fpath = self.get_raw_vcf_fpath_by_callername(callername, gz=True)
        if not isfile(fpath):
            fpath = self.get_raw_vcf_fpath_by_callername(callername, gz=False)
        if not isfile(fpath):
            fpath = self.get_rawest_vcf_fpath_by_callername(callername, gz=True)
        if not isfile(fpath):
            fpath = self.get_rawest_vcf_fpath_by_callername(callername, gz=False)
        return verify_file(fpath)

    def get_raw_vcf_fpath_by_callername(self, callername, gz):
        return join(self.dirpath, BCBioStructure.var_dir,
                    self.name + '-' + callername + '.vcf' + ('.gz' if gz else ''))

    def get_rawest_vcf_fpath_by_callername(self, callername, gz):
        return join(self.dirpath, self.name + '-' + callername + '.vcf' + ('.gz' if gz else ''))

    # annotated vcf
    def find_anno_vcf_by_callername(self, callername):
        fpath = self.get_anno_vcf_fpath_by_callername(callername, gz=True)
        if not verify_file(fpath):
            fpath = self.get_anno_vcf_fpath_by_callername(callername, gz=False)
        return verify_file(fpath)

    def get_anno_vcf_fpath_by_callername(self, callername, gz):
        return join(self.dirpath, BCBioStructure.varannotate_dir,
                    self.name + '-' + callername + '.' + BCBioStructure.anno_vcf_ending +
                    ('.gz' if gz else ''))

    # varqc
    def find_varqc_fpath_by_callername(self, callername):
        fpath = self.get_varqc_fpath_by_callername(callername)
        return verify_file(fpath)

    def get_varqc_fpath_by_callername(self, callername):
        return join(self.dirpath, BCBioStructure.varqc_dir,
                    self.name + '-' + callername + '.' + BCBioStructure.varqc_name + '.html')

    # filtered
    def find_filt_vcf_by_callername(self, callername):
        fpath = self.get_filt_vcf_fpath_by_callername(callername, gz=True)
        if not verify_file(fpath):
            fpath = self.get_filt_vcf_fpath_by_callername(callername, gz=False)
        return verify_file(fpath)

    def get_filt_vcf_fpath_by_callername(self, callername, gz):
        return join(self.dirpath, BCBioStructure.varfilter_dir,
                    self.name + '-' + callername + BCBioStructure.filt_vcf_ending +
                    ('.gz' if gz else ''))

    # varqc after filtering
    def find_varqc_after_fpath_by_callername(self, callername):
        fpath = self.get_varqc_after_fpath_by_callername(callername)
        return verify_file(fpath)

    def get_varqc_after_fpath_by_callername(self, callername):
        return join(self.dirpath, BCBioStructure.varqc_after_dir,
                    self.name + '-' + callername + '.' + BCBioStructure.varqc_after_name + '.html')

    # filtered passed
    def find_pass_filt_vcf_by_callername(self, callername):
        fpath = self.get_pass_filt_vcf_fpath_by_callername(callername, gz=True)
        if not verify_file(fpath):
            fpath = self.get_pass_filt_vcf_fpath_by_callername(callername, gz=False)
        return verify_file(fpath)

    def get_pass_filt_vcf_fpath_by_callername(self, callername, gz):
        return join(self.dirpath, BCBioStructure.varfilter_dir,
                    self.name + '-' + callername + BCBioStructure.pass_filt_vcf_ending +
                    ('.gz' if gz else ''))

    # filtered TSV
    def get_filt_tsv_fpath_by_callername(self, callername):
        return join(self.dirpath, BCBioStructure.varfilter_dir,
                    self.name + '-' + callername + BCBioStructure.filt_tsv_ending)

    # ...other
    def for_json(self):
        return dict((k, v) for k, v in self.__dict__.items() if k != 'bcbio_structure')
    #     return dict(
    #         (k, (v if k != 'vcf_by_caller' else (dict((c.name, v) for c, v in v.items()))))
    #         for k, v in self.__dict__.items())

    @staticmethod
    def load(data, bcbio_structure=None):
        if bcbio_structure:
            data['bcbio_structure'] = bcbio_structure
        sample = BCBioSample(**data)
        sample.__dict__ = data
        return sample


class VariantCaller:
    def __init__(self, name, bcbio_structure=None):
        self.name = self.suf = name
        self.bcbio_structure = bcbio_structure
        self.samples = []

        self.summary_qc_report = None
        self.summary_qc_rep_fpaths = []

        self.single_mut_res_fpath = None
        self.paired_mut_res_fpath = None
        self.single_vcf2txt_res_fpath = None
        self.paired_vcf2txt_res_fpath = None

    def find_fpaths_by_sample(self, dir_name, name, ext, final_dirpaths=None):
        return self._find_files_by_sample(dir_name, '.' + name + '.' + ext, final_dirpaths)

    def find_anno_vcf_by_sample(self):
        return self._find_files_by_sample(BCBioStructure.varannotate_dir, BCBioStructure.anno_vcf_ending)

    def get_filt_vcf_by_sample(self):
        return self._find_files_by_sample(BCBioStructure.varfilter_dir, BCBioStructure.filt_vcf_ending)

    def find_pass_filt_vcf_by_sample(self):
        return self._find_files_by_sample(BCBioStructure.varfilter_dir, BCBioStructure.pass_filt_vcf_ending)

    def _find_files_by_sample(self, dir_name, ending, final_dirpaths=None):
        if final_dirpaths is None:
            final_dirpaths = [self.bcbio_structure.final_dirpath]

        if isinstance(final_dirpaths, basestring):
            final_dirpaths = [final_dirpaths]

        files_by_sample = OrderedDict()

        for s in self.samples:
            for final_dirpath in final_dirpaths:
                fpath = join(
                    final_dirpath,
                    s.name,
                    dir_name,
                    s.name + '-' + self.suf + ending + '.gz')

                if isfile(fpath):
                    if verify_file(fpath):
                        files_by_sample[s.name] = fpath
                else:
                    fpath = fpath[:-3]
                    if isfile(fpath):
                        if verify_file(fpath):
                            files_by_sample[s.name] = fpath
                    elif s.phenotype != 'normal':
                        info('Warning: no ' + fpath + ' for ' + s.name + ', ' + self.name)

        return files_by_sample

    def __str__(self):
        return self.name

    def for_json(self):
        return {k: v for k, v in self.__dict__
                if k not in ['bcbio_structure', 'samples']}


class Batch:
    def __init__(self, name=None):
        self.name = name
        self.normal = None
        self.tumor = []

    def __str__(self):
        return self.name


class BCBioStructure:
    varfilter_name   = varfilter_dir                           = 'varFilter'
    varannotate_name = varannotate_dir                         = 'varAnnotate'
    targetseq_name   = targetseq_dir = targetseq_summary_dir   = 'targetSeq'
    cnv_dir                          = cnv_summary_dir         = 'cnv'

    varqc_name               = 'varQC'
    varqc_after_name         = 'varQC_postVarFilter'
    ngscat_name              = 'ngscat'
    qualimap_name            = 'qualimap'
    targqc_name              = 'targQC'
    fastqc_name              = 'fastqc'

    varqc_repr               = 'Var QC'
    varqc_after_repr         = 'Var QC after filtering'
    ngscat_repr              = 'Ngscat'
    qualimap_repr            = 'Qualimap'
    targqc_repr              = 'Target QC'
    fastqc_repr              = 'FastQC'

    varqc_dir        = varqc_summary_dir       = join('qc', varqc_name)
    varqc_after_dir  = varqc_after_summary_dir = join('qc', varqc_after_name)
    ngscat_dir       = ngscat_summary_dir      = join('qc', ngscat_name)
    qualimap_dir     = qualimap_summary_dir    = join('qc', qualimap_name)
    targqc_summary_dir                         = join('qc', targqc_name)
    fastqc_dir       = fastqc_summary_dir      = join('qc', fastqc_name)

    seq2c_name = 'Seq2C'
    seq2c_seq2cov_ending = 'seq2c_seq2cov.txt'

    var_dir = 'var'
    anno_vcf_ending = '.anno.vcf'
    filt_vcf_ending = '.anno.filt.vcf'
    pass_filt_vcf_ending = '.anno.filt.pass.vcf'
    filt_tsv_ending = '.anno.filt.tsv'

    def __init__(self, cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath=None, proc_name=None):
        self._set_final_dir(bcbio_cnf, bcbio_project_dirpath, final_dirpath)

        self.bcbio_cnf = bcbio_cnf
        self.cnf = cnf
        self.batches = OrderedDefaultDict(Batch)
        self.paired = False
        self.samples = []
        self.variant_callers = OrderedDict()

        self.bed = None
        self.sv_bed = None
        self.project_name = None
        if cnf.project_name:
            self.project_name = cnf.project_name
        else:
            # path is like /ngs/oncology/analysis/bioscience/Bio_0031_Heme_MRL_DLBCL_IRAK4/bcbio_Dev_0079
            if '/ngs/oncology/analysis/' in realpath(bcbio_project_dirpath):
                short_path = realpath(bcbio_project_dirpath).split('/ngs/oncology/analysis/')[1]  # bioscience/Bio_0031_Heme_MRL_DLBCL_IRAK4/bcbio_Dev_0079
                self.project_name = '_'.join(short_path.split('/')[1:])  # Bio_0031_Heme_MRL_DLBCL_IRAK4_bcbio_Dev_0079

        bcbio_project_dirname = basename(bcbio_project_dirpath)  # bcbio_Dev_0079
        bcbio_project_parent_dirname = basename(dirname(bcbio_project_dirpath))  # Bio_0031_Heme_MRL_DLBCL_IRAK4
        if not self.project_name:
            self.project_name = bcbio_project_parent_dirname + '_' + bcbio_project_dirname

        if 'fc_date' not in bcbio_cnf:
            critical('Error: fc_date not in bcbio config!')

        if 'fc_name' in bcbio_cnf:
            if not self.project_name:
                self.project_name = bcbio_cnf['fc_name']
            # Date dirpath is from bcbio and named after fc_name, not our own project name
            self.date_dirpath = join(self.final_dirpath, bcbio_cnf['fc_date'] + '_' + bcbio_cnf['fc_name'])
        else:
            self.date_dirpath = join(self.final_dirpath, bcbio_cnf['fc_date'] + '_' + self.project_name)

        if not verify_dir(self.date_dirpath):
            err('Warning: no project directory of format {fc_date}_{fc_name}, creating ' + self.date_dirpath)
        safe_mkdir(self.date_dirpath)

        cnf.project_name = self.project_name
        info('Project name: ' + self.project_name)
        # self.cnf.name = proc_name or self.project_name

        self.log_dirpath = self.cnf.log_dir = self.cnf.log_dir or join(self.date_dirpath, 'log')
        info('log_dirpath: ' + self.log_dirpath)
        safe_mkdir(self.log_dirpath)
        # info('cnf.name: ' + str(self.cnf.name))

        set_up_log(self.cnf, proc_name, self.project_name, self.final_dirpath)

        self.work_dir = self.cnf.work_dir = self.cnf.work_dir or abspath(join(self.final_dirpath, pardir, 'work', 'post_processing'))
        set_up_work_dir(cnf)

        self.var_dirpath = join(self.date_dirpath, BCBioStructure.var_dir)

        # Moving raw variants in the date dir to var
        for fname in os.listdir(self.date_dirpath):
            if '.vcf' in fname and '.anno.filt' not in fname:
                if not isdir(self.var_dirpath):
                    safe_mkdir(self.var_dirpath)
                src_fpath = join(self.date_dirpath, fname)
                dst_fpath = join(self.var_dirpath, fname)
                info('Moving ' + src_fpath + ' to ' + self.var_dirpath)
                try:
                    os.rename(src_fpath, dst_fpath)
                except OSError:
                    pass

        # cleaning date dir
        for fname in listdir(self.date_dirpath):
            if fname.endswith('.log') or fname in ['project-summary.yaml', 'programs.txt']:
                os.rename(join(self.date_dirpath, fname),
                          join(self.log_dirpath, fname))

        info(' '.join(sys.argv))
        info()
        info('-' * 70)

        # reading sampels
        for sample in [self._read_sample_details(sample_info) for sample_info in bcbio_cnf['details']]:
            if sample.dirpath is None:
                err('For sample ' + sample.name + ', directory does not exist. Thus, skipping that sample.')
            else:
                self.samples.append(sample)

        if not self.samples:
            critical('No directory for any sample. Exiting.')

        # sorting samples
        self.samples.sort(key=lambda _s: _s.key_to_sort())
        for caller in self.variant_callers.values():
            caller.samples.sort(key=lambda _s: _s.key_to_sort())

        self.project_level_report_fpath = join(self.date_dirpath, self.project_name + '.html')
        self.fastqc_summary_fpath =       join(self.date_dirpath, BCBioStructure.fastqc_summary_dir,      BCBioStructure.fastqc_name + '.html')
        self.targqc_summary_fpath =       join(self.date_dirpath, BCBioStructure.targqc_summary_dir,      BCBioStructure.targqc_name + '.html')
        self.varqc_report_fpath =         join(self.date_dirpath, BCBioStructure.varqc_summary_dir,       BCBioStructure.varqc_name + '.html')
        self.varqc_after_report_fpath =   join(self.date_dirpath, BCBioStructure.varqc_after_summary_dir, BCBioStructure.varqc_name + '.html')
        self.varqc_report_fpath_by_caller =       OrderedDict([(k, join(dirname(self.varqc_report_fpath),       k + '.' + basename(self.varqc_report_fpath)))       for k in self.variant_callers.keys()])
        self.varqc_after_report_fpath_by_caller = OrderedDict([(k, join(dirname(self.varqc_after_report_fpath), k + '.' + basename(self.varqc_after_report_fpath))) for k in self.variant_callers.keys()])

        # setting bed files for samples
        if cnf.bed:
            verify_file(cnf.bed)
            self.bed = cnf.bed = adjust_path(cnf.bed)
            self.sv_bed = adjust_path(cnf.bed)
        else:
            bed_files_used = [s.bed for s in self.samples]
            if len(set(bed_files_used)) > 2:
                critical('Error: more than 1 BED file found: ' + str(set(bed_files_used)))
            self.bed = bed_files_used[0] if bed_files_used else None

            sv_bed_files_used = [s.sv_bed for s in self.samples]
            if len(set(sv_bed_files_used)) > 2:
                critical('Error: more than 1 SV BED file found: ' + str(set(sv_bed_files_used)))
            self.sv_bed = sv_bed_files_used[0] if sv_bed_files_used else self.bed

        # setting up batch properties
        for b in self.batches.values():
            if b.normal and b.tumor:
                b.paired = True
                info('Batch ' + b.name + ' is paired')
            else:
                b.paired = False
                info('Batch ' + b.name + ' is single')

        for b in self.batches.values():
            for t_sample in b.tumor:
                t_sample.normal_match = b.normal

        if not self.cnf.verbose:
            info('', ending='')

        # all_variantcallers = set()
        # for s_info in self.bcbio_cnf.details:
        #     all_variantcallers |= set(s_info['algorithm'].get('variantcaller')) or set()

        # samples_fpath = abspath(join(self.cnf.work_dir, 'samples.txt'))
        # with open(samples_fpath, 'w') as f:
        #     for sample_info in self.bcbio_cnf.details:
        #         sample = sample_info['description']
        #         f.write(sample + '\n')

        if not self.cnf.verbose:
            print ''
        else:
            info('Done loading BCBio structure.')

    @staticmethod
    def move_vcfs_to_var(sample):
        fpaths_to_move = []
        for fname in os.listdir(sample.dirpath):
            if any(fname.endswith(ending) for ending in
                   [BCBioStructure.filt_tsv_ending,
                    BCBioStructure.filt_vcf_ending,
                    BCBioStructure.filt_vcf_ending + '.gz',
                    BCBioStructure.filt_vcf_ending + '.idx']):
                continue

            if 'vcf' in fname.split('.') and not (islink(fname) and '.anno.filt' in fname):
                fpaths_to_move.append([sample, fname])

        if fpaths_to_move:
            if not exists(sample.var_dirpath):
                info('Creating "var" directory ' + sample.var_dirpath)
                safe_mkdir(sample.var_dirpath)

        for sample, fname in fpaths_to_move:
            src_fpath = join(sample.dirpath, fname)
            dst_fpath = join(sample.var_dirpath, fname)
            if exists(dst_fpath):
                try:
                    os.remove(dst_fpath)
                except OSError:
                    critical('Cannot move ' + src_fpath + ' to ' + dst_fpath + ': dst exists, and permissions do now allow to remove it.')
                    continue
            safe_mkdir(sample.var_dirpath)
            info('Moving ' + src_fpath + ' to ' + dst_fpath)
            os.rename(src_fpath, dst_fpath)

    def _read_sample_details(self, sample_info):
        sample = BCBioSample(name=str(sample_info['description']).replace('.', '_'), final_dir=self.final_dirpath)

        info('Sample "' + sample.name + '"')
        if not self.cnf.verbose: info(ending='')

        sample.dirpath = adjust_path(join(self.final_dirpath, sample.name))
        if not verify_dir(sample.dirpath):
            sample.dirpath = None
            return sample

        self._set_bam_file(sample)

        if 'min_allele_fraction' in sample_info['algorithm']:
            sample.min_af = float(sample_info['algorithm']['min_allele_fraction']) / 100

        sample.phenotype = None

        sample.genome = sample_info.get('genome_build') or 'hg19'
        if sample.genome not in self.cnf.genomes:
            critical('No section in genomes for ' + sample.genome + ' in ' + self.cnf.sys_cnf)
        # sample.genome = self.cnf.genomes[sample.genome]

        self._set_bed_file(sample, sample_info)

        if 'metadata' in sample_info:
            sample.phenotype = sample_info['metadata'].get('phenotype') or 'tumor'
            info('Phenotype: ' + str(sample.phenotype))

            batch_names = sample_info['metadata']['batch']
            if isinstance(batch_names, basestring):
                batch_names = [batch_names]

            for batch_name in batch_names:
                self.batches[batch_name].name = batch_name

                if sample.phenotype == 'normal':
                    if self.batches[batch_name].normal:
                        critical('Multiple normal samples for batch ' + batch_name)
                    self.batches[batch_name].normal = sample

                elif sample.phenotype == 'tumor':
                    self.batches[batch_name].tumor.append(sample)

        sample.var_dirpath = adjust_path(join(sample.dirpath, BCBioStructure.var_dir))
        # self.move_vcfs_to_var(sample)  # moved to filtering.py

        variantcallers = sample_info['algorithm'].get('variantcaller') or []
        if 'ensemble' in sample_info['algorithm'] and len(variantcallers) >= 2:
            variantcallers.append('ensemble')

        if isinstance(variantcallers, basestring):
            variantcallers = [variantcallers]

        for caller_name in variantcallers:
            info(caller_name)
            vcf_fpath = self._set_vcf_file(caller_name, sample)

            if vcf_fpath:
                caller = self.variant_callers.get(caller_name)
                if not caller:
                    self.variant_callers[caller_name] = VariantCaller(caller_name, self)

                self.variant_callers[caller_name].samples.append(sample)
                sample.vcf_by_callername[caller_name] = vcf_fpath

        info()
        return sample

    def _set_final_dir(self, bcbio_cnf, bcbio_project_dirpath, final_dirpath=None):
        if final_dirpath:
            self.final_dirpath = final_dirpath
        elif 'upload' in bcbio_cnf and 'dir' in bcbio_cnf['upload']:
            final_dirname = bcbio_cnf['upload']['dir']
            self.final_dirpath = adjust_path(join(bcbio_project_dirpath, 'config', final_dirname))
            verify_dir(self.final_dirpath, 'upload directory specified in the bcbio config', is_critical=True)
        else:
            self.final_dirpath = join(bcbio_project_dirpath, 'final')
            if not verify_dir(self.final_dirpath):
                critical('If final directory it is not named "final", please, specify it in the bcbio config.')
        info('Final dirpath: ' + self.final_dirpath)

    def _set_bed_file(self, sample, sample_info):
        bed = None
        if self.cnf.bed:  # Custom BED provided in command line?
            sample.sv_bed = sample.bed = verify_bed(self.cnf.bed, is_critical=True)
            info('BED file for ' + sample.name + ': ' + str(sample.bed))

        else:
            if sample_info['algorithm'].get('variant_regions'):  # Variant regions?
                bed = adjust_path(sample_info['algorithm']['variant_regions'])
                if bed and bed.endswith('.bed'):
                    verify_bed(bed, is_critical=True)
                    sample.bed = bed
                    info('SV BED file for ' + sample.name + ': ' + str(sample.bed))
                else:
                    warn('sv_regions file for ' + sample.name + ' is not BED: ' + str(bed))

            if sample_info['algorithm'].get('sv_regions'):  # SV regions?
                sv_bed = adjust_path(sample_info['algorithm']['sv_regions'])
                if sv_bed and sv_bed.endswith('.bed'):
                    verify_bed(sv_bed, is_critical=True)
                    sample.sv_bed = sv_bed
                    info('SV BED file for ' + sample.name + ': ' + str(sample.sv_bed))
                else:
                    warn('sv_regions file for ' + sample.name + ' is not BED: ' + str(sv_bed))
            else:
                sample.sv_bed = sample.bed

        if not sample.bed and not sample.sv_bed:
            info('No BED file and no sv_regions bed file for ' + sample.name + '. Assuming WGS')

    def _set_bam_file(self, sample):
        bam = adjust_path(join(sample.dirpath, sample.name + '-ready.bam'))
        if isfile(bam) and verify_bam(bam):
            sample.bam = bam
            info('BAM file for ' + sample.name + ': ' + sample.bam)
        else:
            sample.bam = None
            err('No BAM file for ' + sample.name)

    def _set_vcf_file(self, caller_name, sample):
        vcf_fname = sample.name + '-' + caller_name + '.vcf'

        vcf_fpath_gz = adjust_path(join(sample.dirpath, vcf_fname + '.gz'))  # in var
        var_vcf_fpath_gz = adjust_path(join(sample.var_dirpath, vcf_fname + '.gz'))  # in var
        vcf_fpath = adjust_path(join(sample.dirpath, vcf_fname))
        var_vcf_fpath = adjust_path(join(sample.var_dirpath, vcf_fname))  # in var

        if isfile(vcf_fpath_gz):
            verify_file(vcf_fpath_gz, is_critical=True)
            info('Found VCF ' + vcf_fpath_gz)
            return vcf_fpath_gz

        if isfile(var_vcf_fpath_gz):
            verify_file(var_vcf_fpath_gz, is_critical=True)
            info('Found VCF in the var/ dir ' + var_vcf_fpath_gz)
            return var_vcf_fpath_gz

        if isfile(vcf_fpath):
            verify_file(vcf_fpath, is_critical=True)
            info('Found uncompressed VCF ' + var_vcf_fpath)
            return vcf_fpath

        if isfile(var_vcf_fpath):
            verify_file(var_vcf_fpath, is_critical=True)
            info('Found uncompressed VCF in the var/ dir ' + var_vcf_fpath)
            return var_vcf_fpath

        if sample.phenotype != 'normal':
            warn('Warning: no VCF found for ' + sample.name + ', ' + caller_name + ', gzip or uncompressed version in and outsize '
                'the var directory. Phenotype is ' + str(sample.phenotype))
        else:
            info('Notice: no VCF file for ' + sample.name + ', ' + caller_name + ', phenotype ' + str(sample.phenotype))

        return None

    def clean(self):
        for sample in self.samples:
            info('Sample ' + sample.name)

            for dic in [sample.filtered_vcf_by_callername,
                        sample.filtered_tsv_by_callername,
                        sample.filtered_maf_by_callername]:
                for c, fpath in dic.items():
                    try:
                        os.unlink(fpath)
                        info('Removed symlink ' + fpath)
                    except OSError:
                        pass

            for fname in listdir(sample.var_dirpath):
                if not fname.startswith('.'):
                    fpath = join(sample.var_dirpath, fname)
                    os.rename(fpath, join(sample.dirpath, fname))

            for dirname in [BCBioStructure.varannotate_dir,
                            BCBioStructure.varfilter_dir,
                            BCBioStructure.varqc_dir,
                            BCBioStructure.varqc_after_dir,
                            BCBioStructure.ngscat_dir,
                            BCBioStructure.qualimap_dir,
                            BCBioStructure.targetseq_dir,
                            BCBioStructure.var_dir]:
                dirpath = join(sample.dirpath, dirname)
                if isdir(dirpath):
                    info('  removing ' + dirpath)
                    shutil.rmtree(dirpath)
            info()

        for dirname in [BCBioStructure.targetseq_summary_dir,
                        BCBioStructure.cnv_summary_dir,
                        BCBioStructure.varqc_summary_dir,
                        BCBioStructure.varqc_after_summary_dir,
                        BCBioStructure.ngscat_summary_dir,
                        BCBioStructure.qualimap_summary_dir]:
            dirpath = join(self.date_dirpath, dirname)
            if isdir(dirpath):
                info('  removing ' + dirpath)
                shutil.rmtree(dirpath)


def load_bcbio_cnf(cnf, config_dirpath):
    yaml_files_in_config_dir = [
        join(config_dirpath, fname)
        for fname in listdir(config_dirpath)
        if fname.endswith('.yaml')]

    if len(yaml_files_in_config_dir) == 0:
        critical('No YAML file in the config directory.')

    config_fpaths = [
        fpath for fpath in yaml_files_in_config_dir
        if not any(n in fpath for n in ['run_info', 'system_info'])]
    if not config_fpaths:
        critical('No BCBio YAMLs in the config directory ' + config_dirpath +
                 ' (only ' + ', '.join(map(basename, yaml_files_in_config_dir)) + ')')

    yaml_fpath = config_fpaths[0]
    if len(config_fpaths) > 1:
        some_yaml_files = [f for f in config_fpaths if splitext(basename(f))[0] in config_dirpath]
        if len(some_yaml_files) == 0:
            critical('More than one YAML file in the config directory ' +
                     config_dirpath + ': ' + ' '.join(config_fpaths) +
                     ', and no YAML file named after the project.')
        yaml_fpath = some_yaml_files[0]

    yaml_file = abspath(yaml_fpath)

    info('Using bcbio YAML config ' + yaml_file)

    return load_yaml_config(yaml_file)


def _normalize(name):
    return name.lower().replace('_', '').replace('-', '')


def ungzip_if_needed(cnf, fpath):
    if fpath.endswith('.gz'):
        fpath = fpath[:-3]
    if not file_exists(fpath) and file_exists(fpath + '.gz'):
        gz_fpath = fpath + '.gz'
        gunzip = get_system_path(cnf, 'gunzip')

        cmdline = '{gunzip} -c {gz_fpath}'.format(**locals())
        res = call(cnf, cmdline, output_fpath=fpath, exit_on_error=False)
        info()
        if not res:
            return None
    return fpath


# def get_trailing_number(string):
#     m = re.search(r'\d+$', string)
#     return int(m.group()) if m else None


