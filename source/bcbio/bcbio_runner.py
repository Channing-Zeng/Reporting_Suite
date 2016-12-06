import getpass
import os
import hashlib
import base64
import shutil
import traceback
from collections import defaultdict
from os.path import join, exists, dirname, abspath, expanduser, pardir, isfile, isdir, islink, getsize, basename
import datetime
from time import sleep
from traceback import format_exc

import variant_filtering
from yaml import dump
try:
    from yaml import CDumper as Dumper, CLoader as Loader
except ImportError:
    from yaml import Dumper, Loader

import source
from scripts.post.qualimap import get_qualimap_max_mem
from source.bcbio.bcbio_filtering import finish_filtering_for_bcbio
from source.bcbio.bcbio_structure import BCBioStructure
from source.calling_process import call
from source.file_utils import verify_file, add_suffix, symlink_plus, remove_quotes, verify_dir, adjust_path, \
    file_transaction
from source.bcbio.project_level_report import make_report_metadata, get_oncoprints_link, make_multiqc_report
from source.qsub_utils import del_jobs
from source.targetcov.summarize_targetcov import get_bed_targqc_inputs
from source.tools_from_cnf import get_system_path
from source.file_utils import safe_mkdir
from source.logger import info, err, critical, send_email, warn, is_local, CriticalError
from source.targetcov.bam_and_bed_utils import verify_bam, prepare_beds, extract_gene_names_and_filter_exons, verify_bed, \
    check_md5
from source.utils import is_us, md5, is_uk, is_sweden
from source.variants import summarize_qc
from source.variants.filtering import make_vcf2txt_cmdl_params
from source.variants.vcf_processing import verify_vcf
from source.webserver.exposing import sync_with_ngs_server, convert_gpfs_path_to_url
from source.config import defaults, with_cnf
from tools.add_jbrowse_tracks import add_project_files_to_jbrowse


class Step:
    def __init__(self, cnf, run_id, name, script, dir_name=None,
                 interpreter=None, short_name=None, paramln='', env_vars=None,
                 log_fpath_template=None, run_on_chara=False):
        self.name = name
        self.dir_name = dir_name
        self.cnf = cnf
        assert run_id
        self.run_id_ = run_id
        self.short_name = short_name or name
        self.param_line = paramln
        self.run_id = None
        self.script = script
        self.interpreter = interpreter
        self.env_vars = env_vars
        self.log_fpath_template = log_fpath_template
        self.run_on_chara = run_on_chara

    def job_name(self, sample=None, caller=None):
        jn = self.short_name.upper() + '_' + self.run_id_ + \
               ('_' + sample if sample else '') + \
               ('_' + caller if caller else '')
        if not jn[0].isalpha():
            jn = 'j_' + jn
        return jn

    def __repr__(self):
        return self.name


class Steps(list):
    def __init__(self):
        super(Steps, self).__init__()

    @staticmethod
    def normalize(name):
        return name.lower().replace('_', '').replace('-', '')

    @staticmethod
    def contains(name_list, step_name):
        return Steps.normalize(step_name) in [Steps.normalize(name) for name in name_list]

    def __contains__(self, step):
        if step is None:
            return False
        if isinstance(step, Step):
            return Steps.contains([s.name for s in self], step.name)
        else:
            return Steps.contains([s.name for s in self], step)

    def append(self, step):
        if step and not self.__contains__(step.name):
            super(Steps, self).append(step)

    def extend(self, iterable):
        for step in iterable:
            self.append(step)


class JobRunning:
    def __init__(self, step, job_id, sample_name, caller_suf, log_fpath,
                 qsub_cmdline, done_marker_fpath, error_marker_fpath, threads, not_wait=False):
        self.step = step
        self.job_id = job_id
        self.sample_name = sample_name
        self.caller_suf = caller_suf
        self.log_fpath = log_fpath
        self.qsub_cmdline = qsub_cmdline
        self.done_marker = done_marker_fpath
        self.error_marker = error_marker_fpath
        self.repr = step.name
        self.threads = threads
        if sample_name:
            self.repr += ' for ' + sample_name
        if caller_suf:
            self.repr += ', ' + caller_suf
        self.not_wait = not_wait
        self.is_done = False
        self.has_errored = False


# noinspection PyAttributeOutsideInit
class BCBioRunner:
    def __init__(self, cnf, bcbio_structure, bcbio_cnf):
        self.bcbio_structure = bcbio_structure

        self.final_dir = bcbio_structure.final_dirpath
        self.bcbio_cnf = bcbio_cnf
        self.cnf = cnf
        cnf.work_dir = bcbio_structure.work_dir

        filt_cnf_fpath = join(self.bcbio_structure.work_dir, 'filt_cnf.yaml')
        with open(filt_cnf_fpath, 'w') as f:
            dump(cnf.variant_filtering.__dict__, f, default_flow_style=False, Dumper=Dumper)

        user_prid = getpass.getuser()
        timestamp = str(datetime.datetime.now())
        self.run_id = BCBioRunner.__generate_run_id(self.final_dir, bcbio_structure.project_name, user_prid, timestamp)
        info('User PRID: ' + user_prid + ', run_id: ' + self.run_id)

        self.qsub_runner = abspath(expanduser(cnf.qsub_runner))

        self.max_threads = self.cnf.threads
        total_samples_num = len(self.bcbio_structure.samples)
        total_callers_num = total_samples_num * len(self.bcbio_structure.variant_callers)
        self.threads_per_sample = 1  # max(self.max_threads / total_samples_num, 1)

        target_bed, exons_bed, exons_no_genes_bed, genes_fpath, seq2c_bed, original_bed = self._prep_bed()
        cnf.bed = target_bed

        self._init_steps(cnf, self.run_id, target_bed, exons_bed, exons_no_genes_bed,
                         genes_fpath, seq2c_bed, original_bed, filt_cnf_fpath)

        if not cnf.steps:
            cnf.steps = []

        self.steps = Steps()
        if 'Variants' in cnf.steps:
            self.steps.extend([self.varannotate, self.varfilter])
        if Steps.contains(cnf.steps, 'VarAnnotate'):
            self.steps.extend([self.varannotate])
        if Steps.contains(cnf.steps, 'VarFilter'):
            self.steps.extend([self.varfilter])

        if not self.bcbio_structure.is_rnaseq and Steps.contains(cnf.steps, 'TargQC'):
            self.steps.extend([self.targqc])
        # if Steps.contains(cnf.steps, 'AbnormalCovReport'):
        #    self.steps.append(self.abnormal_regions)

        # if Steps.contains(cnf.steps, 'ClinicalReport') or \
        #         Steps.contains(cnf.steps, 'ClinicalReports') or \
        #         Steps.contains(cnf.steps, source.clinreport_name):
        if Steps.contains(cnf.steps, 'ClinicalReport'):
            self.steps.extend([self.clin_report])

        if self.bcbio_structure.is_rnaseq and Steps.contains(cnf.steps, 'Expression'):
            self.steps.extend([self.gene_expression])

        from sys import platform as _platform
        if 'linux' in _platform and not self.cnf.no_bam2bigwig:
            self.steps.append(self.bw_converting)
        if is_us() or is_local() and Steps.contains(cnf.steps, 'Exac'):
            if cnf.bed:
                self.steps.append(self.evaluate_capture)
            self.steps.append(self.prepare_for_exac)

        # self.vardict_steps.extend(
        #     [s for s in [
        #         self.vardict,
        #         self.testsomatic,
        #         self.var_to_vcf_somatic,
        #         self.varqc,
        #         self.varannotate,
        #         self.varfilter_all,
        #         self.varqc_after,
        #         self.varqc_summary,
        #         self.varqc_after_summary
        # ] if contains(s.name, cnf.vardict_steps)])

        info('Final list of steps to run:')
        for s in self.steps:
            info('  ' + s.name)

        self.jobs_running = []

    def _prep_bed(self):
        if self.bcbio_structure.is_rnaseq and not self.bcbio_structure.bed:
            return None, None, None, None, None, None

        info()
        info('Checking BED files')
        self.is_wgs = self.bcbio_structure.is_wgs

        target_bed, exons_bed, genes_fpath = get_bed_targqc_inputs(self.cnf, self.bcbio_structure.bed)
        original_target_bed = target_bed
        seq2c_bed = self.bcbio_structure.sv_bed
        original_seq2c_bed = seq2c_bed

        reuse = self.cnf.reuse_intermediate
        if reuse and target_bed:
            reuse = check_md5(self.cnf.work_dir, adjust_path(target_bed), 'target')
            if reuse:
                info('Target ' + target_bed + ' didn\'t change')
                if exons_bed:
                    reuse = check_md5(self.cnf.work_dir, exons_bed, 'features')
                    if reuse:
                        info('Features ' + exons_bed + ' didn\'t change')

        with with_cnf(self.cnf, reuse_intermediate=reuse) as cnf:
            exons_bed, exons_no_genes_bed, target_bed, seq2c_bed = prepare_beds(cnf, exons_bed, target_bed, seq2c_bed)
            _, _, target_bed, exons_bed, exons_no_genes_bed = \
                extract_gene_names_and_filter_exons(cnf, target_bed, exons_bed, exons_no_genes_bed)

        # exposing
        if target_bed and original_target_bed:
            ready_target_bed = join(self.bcbio_structure.date_dirpath, basename(original_target_bed))
            try:
                shutil.copy(target_bed, ready_target_bed)
            except OSError:
                err(traceback.format_exc())
            info('Target BED file is saved in ' + ready_target_bed)

        if seq2c_bed and original_seq2c_bed:
            cnv_dirpath = join(self.bcbio_structure.date_dirpath, BCBioStructure.cnv_dir)
            safe_mkdir(cnv_dirpath)
            ready_seq2c_bed = join(cnv_dirpath, basename(original_seq2c_bed))
            try:
                shutil.copy(seq2c_bed, ready_seq2c_bed)
            except OSError:
                err(traceback.format_exc())
            info('Seq2C BED file is saved in ' + ready_seq2c_bed)

        return target_bed, exons_bed, exons_no_genes_bed, genes_fpath, seq2c_bed, original_target_bed

    @staticmethod
    def __generate_run_id(final_dir, project_name, prid='', timestamp=''):
        hasher = hashlib.sha1(final_dir + prid + timestamp)
        path_hash = base64.urlsafe_b64encode(hasher.digest()[0:4])[:-2]
        return project_name + '_' + path_hash

    def _init_steps(self, cnf, run_id, target_bed, exons_bed, exons_no_genes_bed,
                    genes_fpath, seq2c_bed, original_bed, filt_cnf_fpath):
        basic_params = \
            ' --sys-cnf ' + self.cnf.sys_cnf + \
            ' --run-cnf ' + self.cnf.run_cnf + \
            ' --project-name ' + self.bcbio_structure.project_name + ' '

        summaries_cmdline_params = \
            basic_params + \
           (' --reuse ' if self.cnf.reuse_intermediate else '') + \
            ' --log-dir -'

        # Params for those who doesn't call bcbio_structure
        params_for_one_sample = \
            basic_params + \
            ' -t ' + str(self.threads_per_sample) + \
           (' --reuse ' if self.cnf.reuse_intermediate else '') + \
           (' --debug ' if self.cnf.debug else '') + \
            ' --log-dir -' + \
            ' --genome {cnf.genome.name}'

        if cnf.email:
            summaries_cmdline_params += ' --email ' + remove_quotes(self.cnf.email) + ' '
            params_for_one_sample += ' --email ' + remove_quotes(self.cnf.email) + ' '

        anno_paramline = params_for_one_sample + ('' +
            ' --vcf \'{vcf}\' {bam_cmdline} {normal_match_cmdline} ' +
            '-o \'{output_dir}\' -s \'{sample}\' -c {caller} --qc ' +
            '--work-dir \'' + join(cnf.work_dir, BCBioStructure.varannotate_name) + '_{sample}_{caller}\' ')
        self.varannotate = Step(cnf, run_id,
            name=BCBioStructure.varannotate_name, short_name='va',
            interpreter='python',
            script=join('scripts', 'post', 'varannotate.py'),
            dir_name=BCBioStructure.varannotate_dir,
            log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', BCBioStructure.varannotate_name + '-{caller}.log'),
            paramln=anno_paramline,
        )
        # self.varqc = Step(cnf, run_id,
        #     name=BCBioStructure.varqc_name, short_name='vq',
        #     interpreter='python',
        #     script=join('scripts', 'post', 'varqc.py'),
        #     dir_name=BCBioStructure.varqc_dir,
        #     log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', BCBioStructure.varqc_name + '-{caller}.log'),
        #     paramln=params_for_one_sample + ' --vcf \'{vcf}\' -o \'{output_dir}\' -s \'{sample}\' -c {caller} '
        #             '--work-dir \'' + join(cnf.work_dir, BCBioStructure.varqc_name) + '_{sample}_{caller}\'',
        # )
        # self.varqc_after = Step(cnf, run_id,
        #     name=BCBioStructure.varqc_after_name, short_name='vqa',
        #     interpreter='python',
        #     script=join('scripts', 'post', 'varqc.py'),
        #     dir_name=BCBioStructure.varqc_after_dir,
        #     log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', BCBioStructure.varqc_after_name + '-{caller}.log'),
        #     paramln=params_for_one_sample + ' --vcf \'{vcf}\' -o \'{output_dir}\' -s \'{sample}\' -c {caller} '
        #             '--work-dir \'' + join(cnf.work_dir, BCBioStructure.varqc_after_name) + '_{sample}_{caller}\' ' +
        #             '--proc-name ' + BCBioStructure.varqc_after_name
        # )


        # self.vcf2txt_single = Step(cnf, run_id,
        #     name='vcf2txt_single', short_name='vcf2txt_single',
        #     interpreter='perl',
        #     script='vcf2txt',
        #     dir_name=BCBioStructure.var_dir,
        #     log_fpath_template=join(self.bcbio_structure.log_dirpath, 'vcf2txt-single-{caller}.log'),
        #     paramln='{paramln}',
        # )
        # self.vcf2txt_paired = Step(cnf, run_id,
        #     name='vcf2txt_paired', short_name='vcf2txt_paired',
        #     interpreter='perl',
        #     script='vcf2txt',
        #     dir_name=BCBioStructure.var_dir,
        #     log_fpath_template=join(self.bcbio_structure.log_dirpath, 'vcf2txt-paired-{caller}.log'),
        #     paramln='{paramln}',
        # )

        # if self.is_wgs:
            # call varfilter.py scripts
            # collect vardict.txt (separately for single/paired), merge
            # collect vardict.PASS.txt (separately for single/paired), merge
            # symlink

        # else:
            # collect vcfs separately for single/paired
            # prepare vcf2txt.pl cmdlines
            # add vcf2txt.pl script to GRID
            # for each samples, add varfilter.py with input as output from vcf2txt.pl and depending on vcf2txt.pl
            # collect vardict.PASS.txt (separately for single/paired), merge
            # symlink
        ##### END FILTERING #####

        if self.bcbio_structure.is_rnaseq:
            self.clin_report = None
            self.targqc = None
            self.abnormal_regions = None
            self.varfilter = None
            self.evaluate_capture = None
            self.bw_converting = None
            self.prepare_for_exac = None
            gene_expression_cmdl = summaries_cmdline_params + ' --genome {cnf.genome.name} ' + self.final_dir
            self.gene_expression = Step(
                cnf, run_id,
                name=BCBioStructure.expression_dir, short_name='expr',
                interpreter='python',
                script=join('scripts', 'post_bcbio', 'gene_expression_summary.py'),
                dir_name=BCBioStructure.expression_dir,
                log_fpath_template=join(self.bcbio_structure.log_dirpath, BCBioStructure.expression_dir  + '.log'),
                paramln=gene_expression_cmdl
            )
        else:
            varfilter_paramline = (' ' +
                ' -o {output_dir} ' +
                ' --output-file {output_file} ' +
                ' --sample {sample} ' +
                ' --caller {caller} ' +
                ' --vcf {vcf} ' +
                ' {vcf2txt_cmdl} ' +
                ' --debug ' +
                ' --project-name ' + self.bcbio_structure.project_name + ' ' +
                ' --genome {cnf.genome.name}' +
                ' --work-dir ' + join(cnf.work_dir, BCBioStructure.varfilter_name) + '_{sample}_{caller} ' +
                ' --dbsnp-multi-mafs ' + cnf.genome.dbsnp_multi_mafs +
                ' --run-info ' + cnf.run_cnf
            )
            if cnf.min_freq is not None:
                varfilter_paramline += ' --min-freq ' + str(cnf.min_freq)

            self.varfilter = Step(cnf, run_id,
                name=BCBioStructure.varfilter_name, short_name='vf',
                interpreter='python',
                script='varfilter',
                dir_name=BCBioStructure.varfilter_dir,
                log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', BCBioStructure.varfilter_name + '-{caller}.log'),
                paramln=varfilter_paramline,
            )

            # targqc *.bam --bed target.bed -g hg19 -o targqc_results -t 3 -s sge -q batch.q -r pename=smp

            targqc_params = (
                ' {bams}' +
                ' --project-name ' + self.bcbio_structure.project_name +
                ' -t ' + str(self.max_threads) +
              (' --reuse ' if self.cnf.reuse_intermediate else '') +
              (' --debug ' if self.cnf.debug else '') +
               ' --genome ' + cnf.genome.name +
               ' -o ' + join(self.bcbio_structure.date_dirpath, BCBioStructure.targqc_dir) +
               ' --work-dir ' + join(self.bcbio_structure.work_dir, BCBioStructure.targqc_name)
            )
            if cnf.queue:
                targqc_params += ' -s sge -q ' + cnf.queue
            if target_bed:
                targqc_params += ' --bed ' + target_bed
            self.targqc = Step(cnf, run_id,
                name=BCBioStructure.targqc_name, short_name='tc',
                script='targqc',
                dir_name=BCBioStructure.targqc_dir,
                log_fpath_template=join(self.bcbio_structure.log_dirpath, BCBioStructure.targqc_name + '.log'),
                paramln=targqc_params
            )

            bcbio_summary = verify_file(join(self.bcbio_structure.date_dirpath, 'project-summary.txt'), is_critical=True)
            clinreport_paramline = ('' +
                ' --sample {sample}' +
                ' -o {output_dir}' +
                ' {bam_cmdl}' +
                ' --bcbio-summary ' + bcbio_summary +
                ' {mutations_cmdl}' +
                ' {var_cmdl}' +
                ' {match_cmdl}' +
                ' {seq2c_cmdl}' +
                ' {sv_cmdl}' +
              ((' --bed ' + target_bed) if target_bed else '') +
               (' --jira ' + self.cnf.jira if self.cnf.jira else '') +
                ' --debug' +
                ' --project-report {project_report_path}' +
                ' -g ' + cnf.genome.name +
                ' --target-type ' + self.bcbio_structure.target_type +
                ' --filt-cnf ' + filt_cnf_fpath +
                ' --project-name ' + self.bcbio_structure.project_name +
                ' --threads ' + str(self.threads_per_sample)
            )

            self.clin_report = Step(cnf, run_id,
                name=source.clinreport_name, short_name='clin',
                interpreter='python',
                script=join('ngs_reporting'),
                dir_name=source.clinreport_dir,
                log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', source.clinreport_name  + '.log'),
                paramln=clinreport_paramline + ' --work-dir ' +
                        join(self.bcbio_structure.work_dir, '{sample}_' + source.clinreport_name)
            )

            self.bw_converting = Step(cnf, run_id,
                name='bam_to_bigwig', short_name='bamtobw',
                interpreter='python',
                script=join('scripts', 'post', 'bam_to_bigwig.py'),
                log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', BCBioStructure.bigwig_name + '.log'),
                paramln=params_for_one_sample + ' -s \'{sample}\' --bam \'{bam}\''
                   ' --work-dir ' + join(self.bcbio_structure.work_dir, '{sample}_' + BCBioStructure.bigwig_name)
            )

            evaluate_capture_cmdline = (summaries_cmdline_params + ' ' + self.final_dir + ' -o ' +
                join(self.bcbio_structure.date_dirpath, 'qc', 'bad_coverage.{min_depth}'))
            self.evaluate_capture = Step(cnf, run_id,
                name='evaluate_capture_target', short_name='capture_eval',
                interpreter='python',
                script=join('tools', 'evaluate_capture_target.py'),
                log_fpath_template=join(self.bcbio_structure.log_dirpath, 'evaluate_capture.{min_depth}.log'),
                paramln=evaluate_capture_cmdline + ' --exac-only-filtering --tricky-regions --min-depth {min_depth}' +
                   ' --work-dir ' + join(self.bcbio_structure.work_dir, 'evaluate_capture.{min_depth}'))

            exac_cmdline = summaries_cmdline_params + ' ' + self.final_dir
            if target_bed:
                exac_cmdline += ' --bed ' + target_bed
            self.prepare_for_exac = Step(cnf, run_id,
                name='prepare_data_for_exac', short_name='exac',
                interpreter='python',
                script=join('tools', 'prepare_data_for_exac.py'),
                log_fpath_template=join(self.bcbio_structure.log_dirpath, 'exac.log'),
                paramln=exac_cmdline + ' --genome {cnf.genome.name}' +
                   ' --work-dir ' + join(self.bcbio_structure.work_dir, 'prepare_data_for_exac')
            )

    def step_log_marker_and_output_paths(self, step, sample_name, caller=None, **kwargs):
        if sample_name:
            base_output_dirpath = abspath(join(self.final_dir, sample_name))
        else:
            base_output_dirpath = abspath(self.bcbio_structure.date_dirpath)

        output_dirpath = join(base_output_dirpath, step.dir_name or '')

        kwargs['caller'] = caller
        log_fpath = step.log_fpath_template.format(**kwargs)
        safe_mkdir(dirname(log_fpath))

        done_markers_dirpath = join(self.bcbio_structure.work_dir, 'done_markers')
        safe_mkdir(done_markers_dirpath)

        marker = join(done_markers_dirpath,
             (step.name + '_' + step.run_id_ +
              ('_' + sample_name if sample_name else '') +
              ('_' + caller if caller else '')))
        return output_dirpath, log_fpath, marker + '.done', marker + '.error'


    def _submit_job(self, step, sample_name='', caller=None, create_dir=True,
                    log_out_fpath=None, wait_for_steps=None, threads=1, mem_m=None, not_wait=False, **kwargs):
        job_name = step.job_name(sample_name, caller)

        for job_id_to_wait in wait_for_steps or []:
            job_to_wait = next((j for j in self.jobs_running if j.job_id == job_id_to_wait), None)
            if not job_to_wait or job_to_wait.has_errored:
                if job_to_wait:
                    warn('Job ' + job_to_wait.job_id + ' has failed, and it required to run this job ' + job_name)
                    warn()

        if sum(j.threads for j in self.jobs_running if not j.is_done) >= self.max_threads:
            self.wait_for_jobs(self.max_threads / 2)  # maximum nubmer of jobs were submitted; waiting for half them to finish

        output_dirpath, log_err_fpath, done_marker_fpath, error_marker_fpath = \
            self.step_log_marker_and_output_paths(step, sample_name, caller, **kwargs)

        if output_dirpath and not isdir(output_dirpath) and create_dir:
            safe_mkdir(join(output_dirpath, pardir))
            safe_mkdir(output_dirpath)

        log_out_fpath = log_out_fpath or log_err_fpath
        safe_mkdir(dirname(log_out_fpath))

        if isfile(log_out_fpath):
            try:
                os.remove(log_out_fpath)
            except OSError:
                err('Warning: cannot remove log stdout file ' + log_out_fpath + ', probably permission denied.')

        if log_err_fpath and isfile(log_err_fpath):
            try:
                os.remove(log_err_fpath)
            except OSError:
                err('Warning: cannot remove log stderr file ' + log_err_fpath + ', probably permission denied.')

        # interpreter = get_system_path(self.cnf, step.interpreter, is_critical=True)
        tool_cmdline = get_system_path(self.cnf, step.interpreter, step.script, is_critical=True)
        if not tool_cmdline: critical('Cannot find: ' + ', '.join(filter(None, [step.interpreter, step.script])))
        params = dict({'output_dir': output_dirpath}.items() + self.__dict__.items() + kwargs.items())
        params['caller'] = caller
        cmdline = tool_cmdline + ' ' + step.param_line.format(**params)

        hold_jid_line = '-hold_jid ' + ','.join(wait_for_steps or ['_'])
        qsub = get_system_path(self.cnf, 'qsub')
        queue = self.cnf.queue
        runner_script = self.qsub_runner
        bash = get_system_path(self.cnf, 'bash')
        extra_qsub_opts = ''
        if step.run_on_chara and is_us():
            extra_qsub_opts += '-l h="chara|rask" '

        mem_opts = ''
        if mem_m and not is_local() and not is_sweden():
            mem_m = min(max(mem_m, 200), 90 * 1024)
            mem = str(int(mem_m)) + 'M'
            if mem_m < 1:
                mem_opts = ''
            else:
                mem_opts = '-l mem_free="' + mem + '" '

        priority = 0
        if self.cnf.qsub_priority:
            priority = self.cnf.qsub_priority

        qsub_cmdline = (
            '{qsub} -pe smp {threads} {mem_opts} {extra_qsub_opts} -S {bash} -q {queue} -p {priority} '
            '-j n -o {log_err_fpath} -e {log_err_fpath} {hold_jid_line} '
            '-N {job_name} {runner_script} {done_marker_fpath} {error_marker_fpath} "{cmdline}"'.format(**locals()))
        # print qsub_cmdline

        if self.cnf.verbose:
            info(step.name + (' - ' + sample_name if sample_name else '') + (' - ' + caller if caller else ''))
            info(qsub_cmdline)
        else:
            print step.name,

        if isfile(done_marker_fpath): os.remove(done_marker_fpath)
        if isfile(error_marker_fpath): os.remove(error_marker_fpath)
        job = JobRunning(step, job_name, sample_name, caller, log_err_fpath, qsub_cmdline,
                         done_marker_fpath, error_marker_fpath, threads=threads, not_wait=not_wait)
        self.jobs_running.append(job)
        call(self.cnf, qsub_cmdline, silent=True, env_vars=step.env_vars, exit_on_error=is_local())

        if self.cnf.verbose: info()
        return output_dirpath


    # def _qualimap_bed(self, bed_fpath):
    #     if self.qualimap in self.steps and bed_fpath:
    #         qualimap_bed_fpath = join(self.cnf.work_dir, 'tmp_qualimap.bed')
    #
    #         fix_bed_for_qualimap(bed_fpath, qualimap_bed_fpath)
    #
    #         return qualimap_bed_fpath
    #     else:
    #         return bed_fpath


    def post_jobs(self):
        self._symlink_cnv()
        info()

        error_msg = []
        try:
            targqc_wait_for_steps = []
            for sample in self.bcbio_structure.samples:
                if not (any(step in self.steps for step in [self.varannotate])):
                    continue

                info(sample.name)
                # Processing VCFs: QC, annotation
                for caller in self.bcbio_structure.variant_callers.values():
                    vcf_fpath = sample.vcf_by_callername.get(caller.name)
                    info()
                    if not vcf_fpath:
                        if sample.phenotype != 'normal':
                            err('VCF does not exist: sample ' + sample.name + ', caller ' + caller.name + '.')
                    else:
                        self._process_vcf(sample, sample.bam, vcf_fpath, caller, threads=self.threads_per_sample)

                info('-' * 70)

            if self.varfilter in self.steps:
                info('Filtering')
                info('Per-sample variant filtering')
                for sample in self.bcbio_structure.samples:
                    if sample.phenotype != 'normal':
                        info('  sample ' + sample.name)
                        for caller in self.bcbio_structure.variant_callers.values():
                            # if not sample.vcf_by_callername.get(caller.name):
                            #     info('    no VCF found for ' + sample.name + ' in ' + str(sample.vcf_by_callername))
                            # else:
                            anno_vcf_fpath = sample.get_anno_vcf_fpath_by_callername(caller.name, gz=True)
                            vcf2txt_cmdl = ''
                            # if self.cohort_mode:
                            #     if sample.normal_match:
                            #         vcf2txt_fpath = caller.paired_vcf2txt_res_fpath
                            #     else:
                            #         vcf2txt_fpath = caller.single_vcf2txt_res_fpath
                            #     vcf2txt_cmdl = ' --vcf2txt ' + vcf2txt_fpath  #sample.get_vcf2txt_by_callername(caller_name)
                            # else:
                            if self.varannotate not in self.steps:
                                if not verify_vcf(anno_vcf_fpath):
                                    critical('Error: VarAnnotate is not in steps, and annotated VCF does not exist: ' + anno_vcf_fpath)

                            wait_for_steps = []
                            # if self.cohort_mode:
                            #     if caller.paired_anno_vcf_by_sample:
                            #         wait_for_steps.append(self.vcf2txt_paired.job_name(caller=caller.name))
                            #     if caller.single_anno_vcf_by_sample:
                            #         wait_for_steps.append(self.vcf2txt_single.job_name(caller=caller.name))
                            # else:
                            wait_for_steps.append(self.varannotate.job_name(sample=sample.name, caller=caller.name))

                            self._submit_job(
                                self.varfilter, sample.name, caller_suf=caller.name,
                                vcf=anno_vcf_fpath, vcf2txt_cmdl=vcf2txt_cmdl, output_file=sample.get_vcf2txt_by_callername(caller.name),
                                threads=1, sample=sample.name, caller=caller.name, genome=sample.genome,
                                wait_for_steps=wait_for_steps)

            # BAMS
            # if self.targqc in self.steps:
            #     info('TargQC')
            #     bams_cmdl = ' '.join(s.bam for s in self.bcbio_structure.samples if s.bam and verify_bam(s.bam))
            #     self._submit_job(
            #         self.targqc,
            #         bams=bams_cmdl,
            #         wait_for_steps=targqc_wait_for_steps)

            # if self.seq2c and self.seq2c in self.steps:
            #     self._submit_job(
            #         self.seq2c,
            #         wait_for_steps=[self.targetcov.job_name(s.name) for s in self.bcbio_structure.samples if self.targetcov in self.steps],
            #         genome=self.bcbio_structure.samples[0].genome)

            if (is_uk() or is_us()) and self.cnf.genome.name.startswith('hg') and self.bw_converting in self.steps and not self.bcbio_structure.is_rnaseq:
                for sample in self.bcbio_structure.samples:
                    if sample.bam and isfile(sample.bam):
                        self._submit_job(self.bw_converting, sample_name=sample.name,
                            sample=sample.name, genome=sample.genome, bam=sample.bam,
                            not_wait=True, mem_m=getsize(sample.bam) * 1.1 / 1024 / 1024 + 500)

            if self.bcbio_structure.is_rnaseq and self.gene_expression in self.steps:
                self._submit_job(self.gene_expression)

            if not self.cnf.verbose:
                print ''
            if self.cnf.verbose:
                info('The following jobs were submitted:')

            if not self.jobs_running:
                info()
                info('No jobs submitted.')
            else:
                msg = ['Submitted jobs for the project ' + self.bcbio_structure.project_name + '. '
                       'Log files for each jobs to track:']
                # lengths = []
                # for job in self.jobs_running:
                #     lengths.append(len(job.name))
                # max_length = max(lengths)

                for job in self.jobs_running:
                    # msg += ' ' * (max_length - len(job.name)) + job.log_fpath)
                    info('  ' + job.repr)

            self.wait_for_jobs()
            info('Finished. Jobs done: ' + str(len([j for j in self.jobs_running if j.is_done])) +
                 ', jobs errored: ' + str(len([j for j in self.jobs_running if j.has_errored])) +
                 ', jobs didn\'t run: ' + str(len([j for j in self.jobs_running if not j.is_done])) +
                 ', total was: ' + str(len([j for j in self.jobs_running])))
            info()

            if self.varfilter in self.steps:
                finish_filtering_for_bcbio(self.cnf, self.bcbio_structure, self.bcbio_structure.variant_callers.values())
                info()

            if self.clin_report in self.steps:
                clinical_report_caller = \
                    self.bcbio_structure.variant_callers.get('vardict') or \
                    self.bcbio_structure.variant_callers.get('vardict-java')
                if not clinical_report_caller:
                    warn('No vardict or vardict-java in the variants callers: ' + ', '.join(self.bcbio_structure.variant_callers.keys()))

                for sample in self.bcbio_structure.samples:
                    var_cmdl = ''
                    mutation_cmdl = ''
                    match_cmdl = ''
                    bam_cmdl = ''
                    sv_cmdl = ''

                    wait_for_steps = []
                    wait_for_steps += [self.targqc.job_name()] if self.targqc in self.steps else []

                    if not sample.phenotype or sample.phenotype != 'normal':
                        match_cmdl = ' --match ' + sample.normal_match.name if sample.normal_match else ''
                        if clinical_report_caller:
                            wait_for_steps += [self.varannotate.job_name(sample.name, caller=clinical_report_caller.name)] if self.varannotate in self.steps else []
                            wait_for_steps += [self.varfilter.job_name(sample.name, caller=clinical_report_caller.name)] if self.varfilter in self.steps else []

                            var_cmdl += ' --vcf ' + sample.get_anno_vcf_fpath_by_callername(clinical_report_caller.name, gz=True)
                            var_cmdl += ' --filt-vcf ' + sample.get_filt_vcf_fpath_by_callername(clinical_report_caller.name, gz=True)

                            variants_fpath = sample.get_vcf2txt_by_callername(clinical_report_caller.name)
                            mut_fpath = add_suffix(variants_fpath, variant_filtering.mut_pass_suffix)
                            if self.varfilter in self.steps or verify_file(mut_fpath):
                                mutation_cmdl = ' --mutations ' + mut_fpath + ' --circos-mutations ' + variants_fpath

                    if sample.bam:
                        bam_cmdl = ' --bam ' + sample.bam
                    # targqc_dirpath = join(self.final_dir, sample.name, BCBioStructure.targqc_dir)
                    # if self.targqc in self.steps or verify_dir(targqc_dirpath):
                    #     targqc_cmdl = ' --targqc ' + join(self.final_dir, sample.name, BCBioStructure.targqc_dir)

                    sv_fpath = sample.find_sv_fpath()
                    if sv_fpath:
                        sv_cmdl = ' --sv ' + sv_fpath
                    sample_dirpath = join(self.bcbio_structure.final_dirpath, sample.name)
                    sample_cnv_dirpath = join(sample_dirpath, BCBioStructure.cnv_dir)
                    if isdir(sample_cnv_dirpath):
                        for fname in os.listdir(sample_cnv_dirpath):
                            if '-manta' in fname and fname.endswith('.vcf.gz'):
                                sv_cmdl += ' --sv-vcf ' + join(sample_cnv_dirpath, fname)

                    seq2c_cmdl = ''
                    seq2c_fpath = join(sample.dirpath, 'cnv', sample.name + '-seq2c.tsv')
                    if verify_file(seq2c_fpath, silent=True):
                        seq2c_cmdl = ' --seq2c ' + seq2c_fpath

                    self._submit_job(
                        self.clin_report,
                        sample.name,
                        sample=sample.name, genome=sample.genome,
                        match_cmdl=match_cmdl, mutations_cmdl=mutation_cmdl,
                        var_cmdl=var_cmdl, bam_cmdl=bam_cmdl,
                        seq2c_cmdl=seq2c_cmdl, sv_cmdl=sv_cmdl,
                        project_report_path=self.bcbio_structure.multiqc_fpath,
                        wait_for_steps=wait_for_steps,
                        threads=self.threads_per_sample)

            self.wait_for_jobs()
            info('NGS oncology reports jobs finished. Jobs done: ' + str(len([j for j in self.jobs_running if j.is_done])) +
                 ', jobs errored: ' + str(len([j for j in self.jobs_running if j.has_errored])) +
                 ', jobs didn\'t run: ' + str(len([j for j in self.jobs_running if not j.is_done])) +
                 ', total was: ' + str(len([j for j in self.jobs_running])))
            info()

            if self.evaluate_capture in self.steps:
                for depth in self.cnf.coverage_reports.exac_depth_thresholds:
                    self._submit_job(
                        self.evaluate_capture,
                        min_depth=str(depth))
            self.wait_for_jobs()

            if self.prepare_for_exac in self.steps:
                info('Exposing to ExAC browser...')
                info()
                self._submit_job(
                    self.prepare_for_exac,
                    not_wait=True)

            if is_us() or is_uk():
                info('Exposing to jBrowse')
                try:
                    add_project_files_to_jbrowse(self.cnf, self.bcbio_structure)
                except:
                    traceback.print_exc()
                    err('Error: cannot export to jBrowse')
                info()

            oncoprints_link = None
            if is_us() and not self.bcbio_structure.is_rnaseq:
                oncoprints_link = get_oncoprints_link(self.cnf,
                    self.bcbio_structure, self.bcbio_structure.project_name)

            info()
            info('Preparing AZ specific metadata for MultiQC')
            metadata_fpath = make_report_metadata(self.cnf,
                bcbio_structure=self.bcbio_structure,
                oncoprints_link=oncoprints_link)
            self.bcbio_structure.multiqc_fpath = make_multiqc_report(self.cnf, self.bcbio_structure, metadata_fpath)

            html_report_fpath = self.bcbio_structure.multiqc_fpath

            html_report_url = None
            if html_report_fpath:
                html_report_url = sync_with_ngs_server(self.cnf,
                    jira_url=self.cnf.jira,
                    project_name=self.bcbio_structure.project_name,
                    sample_names=[s.name for s in self.bcbio_structure.samples],
                    bcbio_final_dirpath=self.bcbio_structure.final_dirpath,
                    summary_report_fpath=html_report_fpath)
                if not html_report_url:
                    if is_us():
                        for key in ['analysis', 'datasets', 'scratch']:
                            if '/' + key + '/' in html_report_fpath:
                                rel_url = html_report_fpath.split('/' + key + '/')[1]
                                html_report_url = join('http://blue.usbod.astrazeneca.net/~klpf990/' + key + '/' + rel_url)
            _final_email_notification(self.cnf, html_report_url or html_report_fpath,
                                      self.cnf.jira, self.bcbio_structure)
            if html_report_url:
                info()
                info('HTML report url: ' + html_report_url)
        except KeyboardInterrupt:
            warn('Interrupted.')
        except SystemExit:
            warn('Interrupted.')
        except CriticalError as e:
            warn('Finished with errors.')
            error_msg = e.args[0]
        finally:
            info('Deleting running jobs...')
            del_jobs(self.cnf, self.jobs_running)
            info('Finishing...')
            if isdir(self.bcbio_structure.final_dirpath):
                change_permissions(self.cnf, self.bcbio_structure.final_dirpath)
            if isdir(self.bcbio_structure.work_dir):
                change_permissions(self.cnf, self.bcbio_structure.work_dir)
            if isdir(join(self.bcbio_structure.work_dir, '..', 'config')):
                change_permissions(self.cnf, join(self.bcbio_structure.work_dir, '..', 'config'))
            info()
            if error_msg:
                err('Done post-processing with errors:')
                err('-' * 70)
                err(error_msg)
            else:
                info('Done post-processing.')
                bcbio_work_dirpath = dirname(self.bcbio_structure.work_dir)

                scratch_root_dirpath = None
                analysis_root_dirpath = None
                if is_us() and '/ngs/oncology/analysis/' in bcbio_work_dirpath:
                    scratch_root_dirpath = '/ngs/scratch/'
                    analysis_root_dirpath = '/ngs/oncology/analysis/'
                elif is_local() and '/Users/vlad/googledrive/az/analysis/' in bcbio_work_dirpath:
                    scratch_root_dirpath = '/Users/vlad/scratch/'
                    analysis_root_dirpath = '/Users/vlad/googledrive/az/analysis/'
                if scratch_root_dirpath:
                    info()
                    work_scratch_dirpath = bcbio_work_dirpath.replace(analysis_root_dirpath, scratch_root_dirpath)
                    if not exists(work_scratch_dirpath):
                        assert work_scratch_dirpath != bcbio_work_dirpath, (work_scratch_dirpath, bcbio_work_dirpath)
                        safe_mkdir(dirname(work_scratch_dirpath))
                        info('Moving work directory to scratch: ' + bcbio_work_dirpath + ' -> ' + work_scratch_dirpath)
                        shutil.move(bcbio_work_dirpath, work_scratch_dirpath)
                        os.symlink(work_scratch_dirpath, bcbio_work_dirpath)
                        info('Symlinked work directory ' + bcbio_work_dirpath + ' -> ' + work_scratch_dirpath)


    def wait_for_jobs(self, number_of_jobs_allowed_to_left_running=0):
        info()
        num_occupied_slots = sum(j.threads for j in self.jobs_running if not j.is_done)
        num_occupied_jobs = sum(1 for j in self.jobs_running if not j.is_done)
        info('Waiting for ' + str(num_occupied_slots - number_of_jobs_allowed_to_left_running) + ' slots to free '
                              'out of ' + str(num_occupied_slots) + ' occupied (by ' + str(num_occupied_jobs) + ' jobs)')
        is_waiting = False  # just we don't want to print info that we are waiting if nothing changed
        while True:
            # set flags for all done jobs
            for j in self.jobs_running:
                if not j.is_done:
                    if isfile(j.done_marker):
                        j.is_done = True
                        if is_waiting: info('', print_date=False)
                        info('Done ' + j.repr)
                    if isfile(j.error_marker):
                        j.is_done = True
                        j.has_errored = True
                        if is_waiting: info('', print_date=False)
                        err('Finished with error: ' + j.repr + '. Please, check the log: ' + str(j.log_fpath))
                    if j.is_done:
                        is_waiting = False

            # check flags and wait if not all are done
            if sum(1 for j in self.jobs_running if not j.is_done and not j.not_wait) <= number_of_jobs_allowed_to_left_running:
                break
            else:
                if not is_waiting:
                    is_waiting = True
                    strs = set()
                    for j in self.jobs_running:
                        if not j.is_done:
                            l = sum(1 for j2 in self.jobs_running if not j2.is_done and j2.step.name == j.step.name)
                            strs.add(j.step.name + ' (' + str(l) + ')')
                    info('Waiting for the jobs to be processed on the cluster (monitor with qstat). '
                         'Jobs running: ' + ', '.join(strs))
                    info('', print_date=True, ending='')
                sleep(20)
                info('.', print_date=False, ending='')


    def _process_vcf(self, sample, bam_fpath, vcf_fpath, caller, threads, steps=None, job_names_to_wait=None):
        steps = steps or self.steps

        bam_cmdline = '--bam ' + bam_fpath if bam_fpath else ''

        normal_match_cmdline = ''
        if sample.normal_match:
            normal_match_cmdline = ' --match-normal-sample-name ' + sample.normal_match.name + ' '

        if self.varannotate in steps:
            self._submit_job(
                self.varannotate, sample.name, caller=caller.name, vcf=vcf_fpath, threads=threads,
                bam_cmdline=bam_cmdline, sample=sample.name,
                genome=sample.genome, normal_match_cmdline=normal_match_cmdline,
                wait_for_steps=job_names_to_wait)

        # if self.varqc in steps:
        #     self._submit_job(
        #         self.varqc, sample.name, caller_suf=caller.name, vcf=sample.get_anno_vcf_fpath_by_callername(caller.name, gz=True),
        #         threads=1, sample=sample.name, caller=caller.name, genome=sample.genome,
        #         wait_for_steps=[self.varannotate.job_name(sample.name, caller.name)] if self.varannotate in self.steps else [])

        # anno_dirpath, _ = self.step_output_dir_and_log_paths(self.varannotate, sample_name, caller=caller_name)
        # annotated_vcf_fpath = join(anno_dirpath, basename(add_suffix(vcf_fpath, 'anno')))

            # wait_for_callers_steps = []
            # for caller in self.bcbio_structure.variant_callers.values():
            #     info('varFilter for ' + caller.name)
            #     if self.is_wgs:  # WGS
            #         wait_for_callers_steps.append(self.varfilter.job_name(caller.name))

    def _symlink_cnv(self):
        cnv_summary_dirpath = join(self.bcbio_structure.date_dirpath, BCBioStructure.cnv_summary_dir)

        for sample in self.bcbio_structure.samples:
            sample_dirpath = join(self.bcbio_structure.final_dirpath, sample.name)
            sample_cnv_dirpath = join(sample_dirpath, BCBioStructure.cnv_dir)

            for cnv_caller in ['-seq2c', '-seq2c-coverage.tsv', '-cnvkit', '-cn_mops',
                               '-lumpy', '-manta', '-wham', '-sv-prioritize', '-metasv']:
                for fname in os.listdir(sample_dirpath):
                    if cnv_caller in fname:
                        # Copy to <sample>/cnv
                        try:
                            safe_mkdir(cnv_summary_dirpath)
                            os.rename(join(sample_dirpath, fname), join(sample_cnv_dirpath, fname))
                        except OSError:
                            err(format_exc())
                            info()

                if isdir(sample_cnv_dirpath):
                    for fname in os.listdir(sample_cnv_dirpath):
                        if cnv_caller in fname:
                            # Symlink to <datestamp>/cnv/<cnvcaller>
                            dst_dirpath = cnv_summary_dirpath
                            dst_fname = fname
                            if sample.name not in fname:
                                dst_fname = sample.name + '.' + dst_fname
                            dst_fpath = join(dst_dirpath, dst_fname)
                            try:
                                safe_mkdir(dst_dirpath)
                                if islink(dst_fpath):
                                    os.unlink(dst_fpath)
                                symlink_plus(join(sample_cnv_dirpath, fname), dst_fpath)
                            except OSError:
                                pass

        for sample in self.bcbio_structure.samples:
            seq2c_fpath = join(cnv_summary_dirpath, sample.name + '-seq2c.tsv')
            if isfile(seq2c_fpath):
                sample.seq2c_fpath = seq2c_fpath

        # Merging all Seq2C into one
        if any(s.seq2c_fpath for s in self.bcbio_structure.samples):
            merged_seq2c_fpath = join(safe_mkdir(cnv_summary_dirpath), BCBioStructure.seq2c_name + '.tsv')
            if not (isfile(merged_seq2c_fpath) and self.cnf.reuse_intermediate):
                with file_transaction(None, merged_seq2c_fpath) as tx:
                    with open(tx, 'w') as out:
                        for file_index, sample in enumerate(self.bcbio_structure.samples):
                            if sample.seq2c_fpath and isfile(sample.seq2c_fpath):
                                with open(sample.seq2c_fpath) as inp:
                                    for line_index, l in enumerate(inp):
                                        if file_index == 0 or line_index > 0:
                                            out.write(l)
            self.bcbio_structure.seq2c_fpath = merged_seq2c_fpath


def _final_email_notification(cnf, html_report_url, jira_url, bs):
    subj = bs.small_project_path or bs.project_name
    txt = 'Post-processing finished for ' + bs.project_name + '\n'
    txt += '\n'
    txt += 'Path: ' + bs.final_dirpath + '\n'
    txt += 'URL: ' + convert_gpfs_path_to_url(bs.final_dirpath) + '\n'
    txt += 'Report: ' + html_report_url + '\n'
    if jira_url:
        txt += 'Jira: ' + jira_url
    send_email(cnf, txt, subj)


def change_permissions(cnf, path):
    try:
        call(cnf, 'chmod -R g+w ' + path, silent=True, exit_on_error=False)
    except:
        pass
