import getpass
import os
import hashlib
import base64
from os.path import join, dirname, abspath, expanduser, pardir, isfile, isdir, islink
import datetime
from time import sleep
from traceback import format_exc

import source
from source.bcbio.bcbio_structure import BCBioStructure
from source.calling_process import call
from source.file_utils import verify_file, add_suffix, symlink_plus, remove_quotes
from source.bcbio.project_level_report import make_project_level_report
from source.qsub_utils import del_jobs
from source.tools_from_cnf import get_system_path
from source.file_utils import safe_mkdir
from source.logger import info, err, critical, send_email, warn, is_local
from source.targetcov.bam_and_bed_utils import verify_bam
from source.utils import is_us
from source.webserver.exposing import sync_with_ngs_server
from source.config import defaults


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
        if isinstance(step, Step):
            return Steps.contains([s.name for s in self], step.name)
        else:
            return Steps.contains([s.name for s in self], step)

    def add_step(self, step):
        if not self.__contains__(step.name):
            self.append(step)

    def extend(self, iterable):
        for step in iterable:
            self.add_step(step)


class JobRunning:
    def __init__(self, step, job_id, sample_name, caller_suf, log_fpath,
                 qsub_cmdline, done_marker_fpath, error_marker_fpath, threads):
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

        user_prid = getpass.getuser()
        timestamp = str(datetime.datetime.now())
        self.run_id = BCBioRunner.__generate_run_id(self.final_dir, bcbio_structure.project_name, user_prid, timestamp)
        info('User PRID: ' + user_prid + ', run_id: ' + self.run_id)

        self.qsub_runner = abspath(expanduser(cnf.qsub_runner))

        self.max_threads = self.cnf.threads
        total_samples_num = len(self.bcbio_structure.samples)
        total_callers_num = total_samples_num * len(self.bcbio_structure.variant_callers)
        self.filtering_threads = min(self.max_threads, total_samples_num)
        if not is_us():
            self.filtering_threads = min(self.max_threads, total_samples_num, 10)
        self.threads_per_sample = 1  # max(self.max_threads / total_samples_num, 1)

        self._init_steps(cnf, self.run_id)

        if not cnf.steps:
            cnf.steps = []

        self.steps = Steps()
        if 'Variants' in cnf.steps:
            self.steps.extend([
                self.varannotate,
                self.varqc,
                self.varqc_summary,
                self.varfilter,
                self.varqc_after,
                self.varqc_after_summary])
        if Steps.contains(cnf.steps, 'VarAnnotate'):
            self.steps.extend([self.varannotate, self.varqc, self.varqc_summary])
        if Steps.contains(cnf.steps, 'VarQC'):
            self.steps.extend([self.varqc, self.varqc_summary, self.varqc_after, self.varqc_after_summary])
        if Steps.contains(cnf.steps, 'VarFilter'):
            self.steps.extend([self.varfilter, self.varqc_after, self.varqc_after_summary])
        if Steps.contains(cnf.steps, 'VarQC_postVarFilter'):
            self.steps.extend([self.varqc_after, self.varqc_after_summary])

        if Steps.contains(cnf.steps, 'TargQC'):
            self.steps.extend([self.targetcov, self.targqc_summary])
        if any(Steps.contains(cnf.steps, name) for name in ['TargetCov', 'TargetSeq']):
            self.steps.extend([self.targetcov, self.targqc_summary])
        # if Steps.contains(cnf.steps, 'Qualimap'):
        #     self.steps.extend([self.qualimap, self.targqc_summary])
        if Steps.contains(cnf.steps, 'ngsCAT'):
            self.steps.extend([self.ngscat, self.targqc_summary])

        if Steps.contains(cnf.steps, 'Seq2C'):
            self.steps.extend([self.seq2c])
        if Steps.contains(cnf.steps, 'AbnormalCovReport'):
            self.steps.append(self.abnormal_regions)

        if Steps.contains(cnf.steps, 'FastQC'):
            self.steps.extend([self.fastqc_summary])

        if Steps.contains(cnf.steps, 'ClinicalReport'):
            self.steps.extend([self.clin_report])

        if Steps.contains(cnf.steps, 'Summary'):
            self.steps.extend([self.varqc_summary, self.varqc_after_summary, self.targqc_summary, self.fastqc_summary])

        # fastqc summary and clinical report -- special case (turn on if user uses default steps)
        if set(defaults['steps']) == set(cnf.steps):
            self.steps.extend([self.fastqc_summary, self.clin_report])

        self.steps.extend([self.varqc_summary, self.varqc_after_summary, self.targqc_summary])

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

    @staticmethod
    def __generate_run_id(final_dir, project_name, prid='', timestamp=''):
        hasher = hashlib.sha1(final_dir + prid + timestamp)
        path_hash = base64.urlsafe_b64encode(hasher.digest()[0:4])[:-2]
        return project_name + '_' + path_hash

    def _init_steps(self, cnf, run_id):
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
            ' --log-dir -' + \
            ' --genome {genome}'

        if cnf.email:
            summaries_cmdline_params += ' --email ' + remove_quotes(self.cnf.email) + ' '
            params_for_one_sample += ' --email ' + remove_quotes(self.cnf.email) + ' '

        anno_paramline = params_for_one_sample + ('' +
            ' --vcf \'{vcf}\' {bam_cmdline} {normal_match_cmdline} ' +
            '-o \'{output_dir}\' -s \'{sample}\' -c {caller} ' +
            '--work-dir \'' + join(cnf.work_dir, BCBioStructure.varannotate_name) + '_{sample}_{caller}\' ')
        # log_fpath = join(self.bcbio_structure.log_dirpath,
        #      (step.name + ('_' + sample_name if sample_name else '') +
        #                   ('_' + caller if caller else '')) + '.log')

        self.varannotate = Step(cnf, run_id,
            name=BCBioStructure.varannotate_name, short_name='va',
            interpreter='python',
            script=join('scripts', 'post', 'varannotate.py'),
            dir_name=BCBioStructure.varannotate_dir,
            log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', BCBioStructure.varannotate_name + '-{caller}.log'),
            paramln=anno_paramline,
        )
        self.varqc = Step(cnf, run_id,
            name=BCBioStructure.varqc_name, short_name='vq',
            interpreter='python',
            script=join('scripts', 'post', 'varqc.py'),
            dir_name=BCBioStructure.varqc_dir,
            log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', BCBioStructure.varqc_name + '-{caller}.log'),
            paramln=params_for_one_sample + ' --vcf \'{vcf}\' -o \'{output_dir}\' -s \'{sample}\' -c {caller} '
                    '--work-dir \'' + join(cnf.work_dir, BCBioStructure.varqc_name) + '_{sample}_{caller}\'',
        )
        self.varqc_after = Step(cnf, run_id,
            name=BCBioStructure.varqc_after_name, short_name='vqa',
            interpreter='python',
            script=join('scripts', 'post', 'varqc.py'),
            dir_name=BCBioStructure.varqc_after_dir,
            log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', BCBioStructure.varqc_after_name + '-{caller}.log'),
            paramln=params_for_one_sample + ' --vcf \'{vcf}\' -o \'{output_dir}\' -s \'{sample}\' -c {caller} '
                    '--work-dir \'' + join(cnf.work_dir, BCBioStructure.varqc_after_name) + '_{sample}_{caller}\' ' +
                    '--proc-name ' + BCBioStructure.varqc_after_name
        )

        targetcov_params = params_for_one_sample + ' --bam \'{bam}\' {bed} -o \'{output_dir}\' ' \
            '-s \'{sample}\' --work-dir \'' + join(cnf.work_dir, BCBioStructure.targqc_name) + '_{sample}\' '
        if cnf.exons:
            targetcov_params += '--exons {cnf.exons} '
        if cnf.reannotate:
            targetcov_params += '--reannotate '
        if cnf.steps and 'AbnormalCovReport' in cnf.steps:
            targetcov_params += '--extended '
        self.targetcov = Step(cnf, run_id,
            name=BCBioStructure.targqc_name, short_name='tc',
            interpreter='python',
            script=join('scripts', 'post', 'targetcov.py'),
            dir_name=BCBioStructure.targqc_dir,
            log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', BCBioStructure.targqc_name + '.log'),
            paramln=targetcov_params,
        )
        self.abnormal_regions = Step(cnf, run_id,
            name='AbnormalCovReport', short_name='acr',
            interpreter='python',
            script=join('scripts', 'post', 'abnormal_regions.py'),
            dir_name=BCBioStructure.targqc_dir,
            log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', 'abnormalRegionsReport.log'),
            paramln=summaries_cmdline_params + ' --mutations {mutations_fpath}' + ' ' + self.final_dir
        )
        self.ngscat = Step(cnf, run_id,
            name=BCBioStructure.ngscat_name, short_name='nc',
            interpreter='python',
            script=join('scripts', 'post', 'ngscat.py'),
            dir_name=BCBioStructure.ngscat_dir,
            log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', BCBioStructure.ngscat_name + '.log'),
            paramln=params_for_one_sample + ' --bam \'{bam}\' --bed \'{bed}\' -o \'{output_dir}\' -s \'{sample}\' '
                    '--saturation y --work-dir \'' + join(cnf.work_dir, BCBioStructure.ngscat_name) + '_{sample}\''
        )
        # self.qualimap = Step(cnf, run_id,
        #     name=BCBioStructure.qualimap_name, short_name='qm',
        #     interpreter='python',
        #     script=join('scripts', 'post', 'qualimap.py'),
        #     dir_name=BCBioStructure.qualimap_dir,
        #     paramln=params_for_one_sample + ' --bam {bam} {bed} -o {output_dir}',
        # )
        #############
        # Summaries #
        self.varqc_summary = Step(cnf, run_id,
            name=BCBioStructure.varqc_name + '_summary', short_name='vqs',
            interpreter='python',
            script=join('scripts', 'post_bcbio', 'varqc_summary.py'),
            dir_name=BCBioStructure.varqc_summary_dir,
            log_fpath_template=join(self.bcbio_structure.log_dirpath, BCBioStructure.varqc_name + '_summary.log'),
            paramln=summaries_cmdline_params + ' ' + self.final_dir
        )
        self.varqc_after_summary = Step(cnf, run_id,
            name=BCBioStructure.varqc_after_name + '_summary', short_name='vqas',
            interpreter='python',
            script=join('scripts', 'post_bcbio', 'varqc_summary.py'),
            dir_name=BCBioStructure.varqc_after_summary_dir,
            log_fpath_template=join(self.bcbio_structure.log_dirpath, BCBioStructure.varqc_after_name + '_summary.log'),
            paramln=summaries_cmdline_params + ' ' + self.final_dir +
                ' --varqc-name ' + BCBioStructure.varqc_after_name +
                ' --varqc-dir ' + BCBioStructure.varqc_after_dir
        )
        varfilter_paramline = summaries_cmdline_params + ' ' + self.final_dir + ' --caller {caller} '
        if cnf.datahub_path:
            varfilter_paramline += ' --datahub-path ' + cnf.datahub_path
        if cnf.min_freq is not None:
            varfilter_paramline += ' --freq ' + str(cnf.min_freq)
        varfilter_paramline += ' -t ' + str(self.filtering_threads) + ' --wgs '

        self.varfilter = Step(cnf, run_id,
            name=BCBioStructure.varfilter_name, short_name='vfs',
            interpreter='python',
            script=join('scripts', 'post_bcbio', 'varfilter.py'),
            dir_name=BCBioStructure.varfilter_dir,
            log_fpath_template=join(self.bcbio_structure.log_dirpath, BCBioStructure.varfilter_name + '-{caller}.log'),
            paramln=varfilter_paramline,
            run_on_chara=True,
        )

        clinreport_paramline = (params_for_one_sample +
           ' --targqc-dir ' + join(self.final_dir, '{sample}', BCBioStructure.targqc_dir) +
           ' --mutations {mutations_fpath} -s {sample}' +
           ' --varqc {varqc}' +
           ' --seq2c ' + self.bcbio_structure.seq2c_fpath +
           ' --target-type ' + self.bcbio_structure.target_type +
          (' --bed ' + self.bcbio_structure.bed if self.bcbio_structure.bed else '') +
           ' -s {sample} -o {output_dir} ' +
           ' {match_cmdl} ' +
           ' --project-level-report {project_report_path}' +
           ' --work-dir ' + join(self.bcbio_structure.work_dir, '{sample}_' + source.clinreport_name))
        self.clin_report = Step(cnf, run_id,
            name=source.clinreport_name, short_name='clin',
            interpreter='python',
            script=join('scripts', 'post', 'clinical_report.py'),
            dir_name=source.clinreport_dir,
            log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', source.clinreport_name  + '.log'),
            paramln=clinreport_paramline
        )

        self.mongo_loader = Step(cnf, run_id,
            name='MongoLoader', short_name='ml',
            interpreter='java',
            script='vcf_loader',
            dir_name='mongo_loader',
            log_fpath_template=join(self.bcbio_structure.log_dirpath, 'mongo_loader.log'),
            paramln='-module loader -project {project} -sample {sample} -path {path} -variantCaller {variantCaller}'
        )

        seq2c_cmdline = summaries_cmdline_params + ' ' + self.final_dir + ' --genome {genome} '
        if self.bcbio_structure.bed:
            seq2c_cmdline += ' --bed ' + self.bcbio_structure.bed
        normal_snames = [b.normal.name for b in self.bcbio_structure.batches.values() if b.normal]
        if normal_snames or cnf.seq2c_controls:
            controls = (normal_snames or []) + (cnf.seq2c_controls.split(':') if cnf.seq2c_controls else [])
            seq2c_cmdline += ' -c ' + ':'.join(controls)
        if cnf.seq2c_opts:
            seq2c_cmdline += ' --seq2c_opts ' + cnf.seq2c_opts
        if cnf.reannotate:
            seq2c_cmdline += ' --reannotate '
        self.seq2c = Step(cnf, run_id,
            name=BCBioStructure.seq2c_name, short_name='seq2c',
            interpreter='python',
            script=join('scripts', 'post_bcbio', 'seq2c.py'),
            log_fpath_template=join(self.bcbio_structure.log_dirpath, BCBioStructure.seq2c_name + '.log'),
            dir_name=BCBioStructure.cnv_summary_dir,
            paramln=seq2c_cmdline,
            run_on_chara=True
        )

        targqc_cmdline = summaries_cmdline_params + ' ' + self.final_dir
        if self.bcbio_structure.bed:
            targqc_cmdline += ' --bed ' + self.bcbio_structure.bed

        self.targqc_summary = Step(cnf, run_id,
            name=BCBioStructure.targqc_name + '_summary', short_name='targqc',
            interpreter='python',
            script=join('scripts', 'post_bcbio', 'targqc_summary.py'),
            log_fpath_template=join(self.bcbio_structure.log_dirpath, BCBioStructure.targqc_name + '_summary.log'),
            dir_name=BCBioStructure.targqc_summary_dir,
            paramln=targqc_cmdline
        )
        self.fastqc_summary = Step(cnf, run_id,
            name=BCBioStructure.fastqc_name, short_name='fastqc',
            interpreter='python',
            script=join('scripts', 'post_bcbio', 'fastqc_summary.py'),
            log_fpath_template=join(self.bcbio_structure.log_dirpath, BCBioStructure.fastqc_name + '_summary.log'),
            dir_name=BCBioStructure.fastqc_summary_dir,
            paramln=summaries_cmdline_params + ' ' + self.final_dir
        )


    def step_log_marker_and_output_paths(self, step, sample_name, caller=None):
        if sample_name:
            base_output_dirpath = abspath(join(self.final_dir, sample_name))
        else:
            base_output_dirpath = abspath(self.bcbio_structure.date_dirpath)

        output_dirpath = join(base_output_dirpath, step.dir_name) if step.dir_name else ''

        log_fpath = step.log_fpath_template.format(sample=sample_name, caller=caller)
        safe_mkdir(dirname(log_fpath))

        done_markers_dirpath = join(self.bcbio_structure.work_dir, 'done_markers')
        safe_mkdir(done_markers_dirpath)

        marker = join(done_markers_dirpath,
             (step.name + '_' + step.run_id_ +
              ('_' + sample_name if sample_name else '') +
              ('_' + caller if caller else '')))
        return output_dirpath, log_fpath, marker + '.done', marker + '.error'


    def _submit_job(self, step, sample_name='', caller_suf=None, create_dir=True,
                    log_out_fpath=None, wait_for_steps=None, threads=1, **kwargs):
        job_name = step.job_name(sample_name, caller_suf)

        for jn in wait_for_steps or []:
            j = next((j for j in self.jobs_running if j.job_id == jn), None)
            if j and j.has_errored:
                warn('Job ' + j.job_id + ' has failed, and it required to run this job ' + job_name)
                return None

        if sum(j.threads for j in self.jobs_running if not j.is_done) >= self.max_threads:
            self.wait_for_jobs(self.max_threads / 2)  # maximum nubmer of jobs were submitted; waiting for half them to finish

        output_dirpath, log_err_fpath, done_marker_fpath, error_marker_fpath = \
            self.step_log_marker_and_output_paths(step, sample_name, caller_suf)

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
        cmdline = tool_cmdline + ' ' + step.param_line.format(**params)

        hold_jid_line = '-hold_jid ' + ','.join(wait_for_steps or ['_'])
        qsub = get_system_path(self.cnf, 'qsub')
        mem = str(threads * 15)
        queue = self.cnf.queue
        runner_script = self.qsub_runner
        bash = get_system_path(self.cnf, 'bash')
        extra_qsub_opts = ''
        if step.run_on_chara and is_us():
            extra_qsub_opts += '-l h="chara|rask" '
        qsub_cmdline = (
            '{qsub} -pe smp {threads} {extra_qsub_opts} -S {bash} -q {queue} -j n '
            '-o {log_err_fpath} -e {log_err_fpath} {hold_jid_line} -N {job_name} '
            '{runner_script} {done_marker_fpath} {error_marker_fpath} '
            '"{cmdline}"'.format(**locals()))
        # print qsub_cmdline

        if self.cnf.verbose:
            info(step.name + (' - ' + sample_name if sample_name else '') + (' - ' + caller_suf if caller_suf else ''))
            info(qsub_cmdline)
        else:
            print step.name,

        if isfile(done_marker_fpath): os.remove(done_marker_fpath)
        if isfile(error_marker_fpath): os.remove(error_marker_fpath)
        job = JobRunning(step, job_name, sample_name, caller_suf, log_err_fpath, qsub_cmdline,
                         done_marker_fpath, error_marker_fpath, threads=threads)
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

        callers = self.bcbio_structure.variant_callers.values()

        try:
            targqc_wait_for_steps = []
            for sample in self.bcbio_structure.samples:
                if not (any(step in self.steps for step in
                            [self.targetcov,
                             self.seq2c,
                             self.ngscat,
                             self.varqc,
                             self.varqc_after,
                             self.varannotate,
                             self.mongo_loader,
                             self.abnormal_regions])):
                    continue

                info('Processing "' + sample.name + '"')
                if not self.cnf.verbose:
                    info(ending='')

                # BAMS
                if any(step in self.steps for step in [
                       self.targetcov,
                       self.ngscat]):
                    if not sample.bam or not verify_bam(sample.bam):
                        err('Cannot run coverage reports (targetcov, qualimap, ngscat) without BAM files.')
                    else:
                        # TargetCov reports
                        if self.targetcov in self.steps:
                            info('Target coverage for "' + sample.name + '"')
                            self._submit_job(
                                self.targetcov, sample.name,
                                bam=sample.bam, bed=(('--bed ' + self.bcbio_structure.bed) if self.bcbio_structure.bed else ''), sample=sample.name, genome=sample.genome,
                                caller_names='', vcfs='', threads=self.threads_per_sample, wait_for_steps=targqc_wait_for_steps)
                            # if not sample.bed:  # WGS
                            #     targqc_wait_for_steps.append(self.targetcov.job_name(sample.name))

                        # ngsCAT reports
                        if (self.ngscat in self.steps) and (not sample.bed or not verify_file(sample.bed)):
                            warn('Warning: no BED file, assuming WGS, thus skipping ngsCAT reports.')
                        else:
                            if self.ngscat in self.steps:
                                self._submit_job(
                                    self.ngscat, sample.name, bam=sample.bam, bed=self.bcbio_structure.bed or self.cnf.genomes[sample.genome].exons,
                                    sample=sample.name, genome=sample.genome, threads=self.threads_per_sample)

                # Processing VCFs: QC, annotation
                for caller in self.bcbio_structure.variant_callers.values():
                    vcf_fpath = sample.vcf_by_callername.get(caller.name)
                    if not vcf_fpath:
                        if sample.phenotype != 'normal':
                            err('VCF does not exist: sample ' + sample.name + ', caller ' + caller.name + '.')
                    else:
                        self._process_vcf(sample, sample.bam, vcf_fpath, caller.name, threads=self.threads_per_sample)

                if self.cnf.verbose:
                    info('-' * 70)
                else:
                    print ''
                    info()

            if not self.cnf.verbose:
                info('', ending='')

            if self.seq2c in self.steps:
                self._submit_job(
                    self.seq2c,
                    wait_for_steps=[self.targetcov.job_name(s.name) for s in self.bcbio_structure.samples if self.targetcov in self.steps],
                    genome=self.bcbio_structure.samples[0].genome)

            if self.targqc_summary in self.steps:
                wait_for_steps = []
                wait_for_steps += [self.targetcov.job_name(s.name) for s in self.bcbio_structure.samples if self.targetcov in self.steps]
                wait_for_steps += [self.ngscat.job_name(s.name) for s in self.bcbio_structure.samples if self.ngscat in self.steps]
                self._submit_job(
                    self.targqc_summary,
                    wait_for_steps=wait_for_steps)

            if self.fastqc_summary in self.steps:
                self._submit_job(self.fastqc_summary)

            if self.varqc_summary in self.steps:
                self._submit_job(
                    self.varqc_summary,
                    wait_for_steps=[
                        self.varqc.job_name(s.name, v.name)
                        for v in self.bcbio_structure.variant_callers.values()
                        for s in v.samples
                        if self.varqc in self.steps])

            if self.varfilter in self.steps:
                wait_for_callers_steps = []
                for caller in self.bcbio_structure.variant_callers.values():
                    info('varFilter for ' + caller.name)
                    self._submit_job(
                        self.varfilter,
                        caller_suf=caller.name, caller=caller.name,
                        wait_for_steps=[
                            self.varannotate.job_name(s.name, caller.name)
                            for s in caller.samples
                            if self.varannotate in self.steps] + wait_for_callers_steps,
                        create_dir=False,
                        threads=self.filtering_threads)
                    if not self.bcbio_structure.bed:  # WGS
                        wait_for_callers_steps.append(self.varfilter.job_name(caller.name))

            if self.clin_report in self.steps:
                clinical_report_caller = \
                    self.bcbio_structure.variant_callers.get('vardict') or \
                    self.bcbio_structure.variant_callers.get('vardict-java')

                if clinical_report_caller:
                    vardict_txt_fname = source.mut_fname_template.format(caller_name=clinical_report_caller.name)
                    vardict_txt_fpath = join(self.bcbio_structure.date_dirpath, vardict_txt_fname)
                    mutations_fpath = add_suffix(vardict_txt_fpath, source.mut_pass_suffix)
                    
                    for sample in clinical_report_caller.samples:
                        wait_for_steps = []
                        wait_for_steps += [self.targetcov.job_name(sample.name)] if self.targetcov in self.steps else []
                        wait_for_steps += [self.varqc.job_name(sample.name, caller=clinical_report_caller.name)] if self.varqc in self.steps else []
                        wait_for_steps += [self.varfilter.job_name(caller=clinical_report_caller.name)] if self.varfilter in self.steps else []
                        wait_for_steps += [self.seq2c.job_name()] if self.seq2c in self.steps else []
                        match_cmdl = ' --match ' + sample.normal_match.name if sample.normal_match else ''
                        self._submit_job(
                            self.clin_report,
                            sample.name,
                            sample=sample.name, match_cmdl=match_cmdl, genome=sample.genome, mutations_fpath=mutations_fpath,
                            varqc=sample.get_varqc_fpath_by_callername(clinical_report_caller.name, ext='.json'),
                            project_report_path=self.bcbio_structure.project_report_html_fpath,
                            wait_for_steps=wait_for_steps,
                            threads=self.threads_per_sample)
                else:
                    warn('Warning: Clinical report cannot be created.')

            # TargetSeq reports
            if self.abnormal_regions in self.steps:
                variant_caller = \
                    self.bcbio_structure.variant_callers.get('vardict') or \
                    self.bcbio_structure.variant_callers.get('vardict-java')
                vardict_txt_fname = source.mut_fname_template.format(caller_name=variant_caller.name)
                vardict_txt_fpath = join(self.bcbio_structure.date_dirpath, vardict_txt_fname)
                mutations_fpath = add_suffix(vardict_txt_fpath, source.mut_pass_suffix)

                for sample in self.bcbio_structure.samples:
                    if not self.cnf.verbose:
                        info(ending='')

                    callers_and_filtered_vcfs = [(c, f) for c, f in ((c.name, c.get_filt_vcf_by_sample().get(sample.name)) for c in callers) if f]
                    if callers_and_filtered_vcfs:
                        caller_names, filtered_vcfs = zip(*callers_and_filtered_vcfs)
                    else:
                        caller_names, filtered_vcfs = [], []

                    wait_for_steps = []
                    if self.varfilter in self.steps:
                        wait_for_steps.extend([self.varfilter.job_name(caller=caller.name)])
                    if self.targetcov in self.steps:
                        wait_for_steps.extend([self.targetcov.job_name(sample.name)])

                    self._submit_job(
                        self.abnormal_regions, sample.name,
                        wait_for_steps=wait_for_steps,
                        sample=sample, threads=self.threads_per_sample, genome=sample.genome, mutations_fpath=mutations_fpath,
                        caller_names='--caller-names ' + ','.join(caller_names) if caller_names else '',
                        vcfs='--vcfs ' + ','.join(filtered_vcfs) if filtered_vcfs else '')

            if self.varqc_after in self.steps:
                info('VarQC_postVarFilter:')
                for caller in self.bcbio_structure.variant_callers.values():
                    info('caller: ' + caller.name)
                    for sample in caller.samples:
                        info('sample: ' + sample.name)
                        raw_vcf_fpath = sample.find_raw_vcf_by_callername(caller.name)
                        if not raw_vcf_fpath:
                            if sample.phenotype != 'normal':
                                err('Error: raw VCF does not exist: sample ' + sample.name + ', caller "' +
                                    caller.name + '". Phenotype = ' + sample.phenotype + '.')
                        else:
                            filt_vcf_fpath = sample.get_filt_vcf_fpath_by_callername(caller.name, gz=True)
                            if not self.varfilter and sample.phenotype != 'normal' and not verify_file(filt_vcf_fpath, silent=True):
                                err('Error: filtered VCF does not exist: sample ' + sample.name + ', caller "' +
                                    caller.name + '". Phenotype = ' + sample.phenotype + '.' +
                                    ' Note that you need to run VarFilter first, and this step is not in config.')
                            else:
                                self._submit_job(
                                    self.varqc_after, sample.name, caller_suf=caller.name, threads=self.threads_per_sample,
                                    wait_for_steps=([self.varfilter.job_name(caller=caller.name)] if self.varfilter in self.steps else []),
                                    vcf=filt_vcf_fpath, sample=sample.name, caller=caller.name, genome=sample.genome)

            if self.varqc_after_summary in self.steps:
                self._submit_job(
                    self.varqc_after_summary,
                    wait_for_steps=[
                        self.varqc_after.job_name(s.name, v.name)
                        for v in self.bcbio_structure.variant_callers.values()
                        for s in v.samples
                        if self.varqc_after in self.steps])

            # if self.combined_report in self.steps:
            #     wait_for_steps = []
            #     # summaries
            #     wait_for_steps += [self.varqc_summary.job_name()] if self.varqc_summary in self.steps else []
            #     wait_for_steps += [self.varqc_after_summary.job_name()] if self.varqc_after_summary in self.steps else []
            #     wait_for_steps += [self.targqc_summary.job_name()] if self.targqc_summary in self.steps else []
            #     wait_for_steps += [self.fastqc_summary.job_name()] if self.fastqc_summary in self.steps else []
            #     # and individual reports too
            #     wait_for_steps += [self.varqc.job_name(s.name) for s in self.bcbio_structure.samples if self.varqc in self.steps]
            #     wait_for_steps += [self.targetcov.job_name(s.name) for s in self.bcbio_structure.samples if self.targetcov in self.steps]
            #     wait_for_steps += [self.ngscat.job_name(s.name) for s in self.bcbio_structure.samples if self.ngscat in self.steps]
            #     wait_for_steps += [self.qualimap.job_name(s.name) for s in self.bcbio_structure.samples if self.qualimap in self.steps]
            #     self._submit_job(
            #         self.combined_report,
            #         wait_for_steps=wait_for_steps)

            # if self.mongo_loader in self.steps:
            #     for sample in self.bcbio_structure.samples:
            #         for caller in self.bcbio_structure.variant_callers.values():
            #             filt_vcf_fpath = sample.find_filt_vcf_by_callername(caller.name)
            #             self._submit_job(
            #                 self.mongo_loader, sample.name, suf=caller.name, create_dir=False,
            #                 wait_for_steps=([self.varfilter_all.job_name()] if self.varfilter_all in self.steps else []),
            #                 path=filt_vcf_fpath, sample=sample.name, variantCaller=caller.name,
            #                 project=self.bcbio_structure.project_name)

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
                 ', total was: ' + str(len([j for j in self.jobs_running]))
            )

            html_report_fpath = make_project_level_report(self.cnf, bcbio_structure=self.bcbio_structure)

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
            _final_email_notification(html_report_url, self.cnf.jira, self.bcbio_structure)
            if html_report_url:
                info()
                info('HTML report url: ' + html_report_url)

        except:
            raise
        finally:
            del_jobs(self.cnf, self.jobs_running)


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
                        info('Finished with error: ' + j.repr + '. Please, check the log: ' + str(j.log_fpath))
                    if j.is_done:
                        is_waiting = False

            # check flags and wait if not all are done
            if sum(1 for j in self.jobs_running if not j.is_done) <= number_of_jobs_allowed_to_left_running:
                break
            else:
                if not is_waiting:
                    is_waiting = True
                    strs = []
                    for j in self.jobs_running:
                        if not j.is_done:
                            l = sum(1 for j2 in self.jobs_running if not j2.is_done and j2.step.name == j.step.name)
                            strs.append(j.step.name + ' (' + str(l) + ')')
                    info('Waiting for the jobs to be processed on a GRID (monitor with qstat). Jobs running: ' + ', '.join(strs))
                    info('', print_date=True, ending='')
                sleep(20)
                info('.', print_date=False, ending='')


    def _process_vcf(self, sample, bam_fpath, vcf_fpath, caller_name, threads,
                     steps=None, job_names_to_wait=None):
        steps = steps or self.steps

        bam_cmdline = '--bam ' + bam_fpath if bam_fpath else ''

        normal_match_cmdline = ''
        if sample.normal_match:
            normal_match_cmdline = ' --match-normal-sample-name ' + sample.normal_match.name + ' '

        if self.varannotate in steps:
            self._submit_job(
                self.varannotate, sample.name, caller_suf=caller_name, vcf=vcf_fpath, threads=threads,
                bam_cmdline=bam_cmdline, sample=sample.name, caller=caller_name,
                genome=sample.genome, normal_match_cmdline=normal_match_cmdline,
                wait_for_steps=job_names_to_wait)

        if self.varqc in steps:
            self._submit_job(
                self.varqc, sample.name, caller_suf=caller_name, vcf=sample.get_anno_vcf_fpath_by_callername(caller_name, gz=True),
                threads=threads, sample=sample.name, caller=caller_name, genome=sample.genome,
                wait_for_steps=[self.varannotate.job_name(sample.name, caller_name)] if self.varannotate in self.steps else [])

        # anno_dirpath, _ = self.step_output_dir_and_log_paths(self.varannotate, sample_name, caller=caller_name)
        # annotated_vcf_fpath = join(anno_dirpath, basename(add_suffix(vcf_fpath, 'anno')))

    def _symlink_cnv(self):
        cnv_summary_dirpath = join(self.bcbio_structure.date_dirpath, BCBioStructure.cnv_summary_dir)
        try:
            safe_mkdir(cnv_summary_dirpath)
        except OSError:
            pass

        for sample in self.bcbio_structure.samples:
            sample_dirpath = join(self.bcbio_structure.final_dirpath, sample.name)
            sample_cnv_dirpath = join(sample_dirpath, BCBioStructure.cnv_dir)

            for cnv_caller in ['-cn_mops', '-cnvkit', '-lumpy', '-manta', '-wham']:
                for fname in os.listdir(sample_dirpath):
                    if cnv_caller in fname:
                        # Copy to <sample>/cnv
                        safe_mkdir(sample_cnv_dirpath)
                        try:
                            os.rename(join(sample_dirpath, fname), join(sample_cnv_dirpath, fname))
                        except OSError:
                            err(format_exc())
                            info()
                            pass

                        # Symlink to <datestamp>/cnv/<cnvcaller>
                        dst_dirpath = join(cnv_summary_dirpath, cnv_caller[1:])
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


def _final_email_notification(html_report_url, jira_url, bs):
    subj = bs.small_project_path or bs.project_name
    txt = 'Post-processing finished for ' + bs.project_name + '\n'
    txt += '\n'
    txt += 'Path: ' + bs.final_dirpath + '\n'
    txt += 'Report: ' + (html_report_url or bs.project_report_html_fpath) + '\n'
    if jira_url:
        txt += 'Jira: ' + jira_url
    send_email(txt, subj)

