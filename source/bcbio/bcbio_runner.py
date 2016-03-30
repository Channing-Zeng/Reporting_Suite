import getpass
import os
import hashlib
import base64
import shutil
import traceback
from collections import defaultdict
from genericpath import exists
from os.path import join, dirname, abspath, expanduser, pardir, isfile, isdir, islink, getsize
import datetime
from time import sleep
from traceback import format_exc

import source
from scripts.post.qualimap import get_qualimap_max_mem
from source.bcbio.bcbio_filtering import finish_filtering_for_bcbio
from source.bcbio.bcbio_structure import BCBioStructure
from source.calling_process import call
from source.fastqc.summarize_fastqc import write_fastqc_combo_report
from source.file_utils import verify_file, add_suffix, symlink_plus, remove_quotes, verify_dir, adjust_path
from source.bcbio.project_level_report import make_project_level_report, get_run_info
from source.qsub_utils import del_jobs
from source.targetcov.summarize_targetcov import get_bed_targqc_inputs
from source.tools_from_cnf import get_system_path
from source.file_utils import safe_mkdir
from source.logger import info, err, critical, send_email, warn, is_local
from source.targetcov.bam_and_bed_utils import verify_bam, prepare_beds, extract_gene_names_and_filter_exons, verify_bed, \
    check_md5
from source.utils import is_us, md5, is_uk
from source.variants import summarize_qc
from source.variants.filtering import make_vcf2txt_cmdl_params
from source.variants.vcf_processing import verify_vcf
from source.webserver.exposing import sync_with_ngs_server, convert_path_to_url
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

    def add_step(self, step):
        if not self.__contains__(step.name):
            self.append(step)

    def extend(self, iterable):
        for step in iterable:
            self.add_step(step)


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

        user_prid = getpass.getuser()
        timestamp = str(datetime.datetime.now())
        self.run_id = BCBioRunner.__generate_run_id(self.final_dir, bcbio_structure.project_name, user_prid, timestamp)
        info('User PRID: ' + user_prid + ', run_id: ' + self.run_id)

        self.qsub_runner = abspath(expanduser(cnf.qsub_runner))

        self.max_threads = self.cnf.threads
        total_samples_num = len(self.bcbio_structure.samples)
        total_callers_num = total_samples_num * len(self.bcbio_structure.variant_callers)
        self.threads_per_sample = 1  # max(self.max_threads / total_samples_num, 1)

        self._init_steps(cnf, self.run_id)

        if not cnf.steps:
            cnf.steps = []

        self.steps = Steps()
        if 'Variants' in cnf.steps:
            self.steps.extend([
                self.varannotate,
                self.varfilter])
        if Steps.contains(cnf.steps, 'VarAnnotate'):
            self.steps.extend([self.varannotate])
        if Steps.contains(cnf.steps, 'VarFilter'):
            self.steps.extend([self.varfilter])

        if Steps.contains(cnf.steps, 'TargQC'):
            self.steps.extend([self.targetcov, self.targqc_summary, self.abnormal_regions])
        if any(Steps.contains(cnf.steps, name) for name in ['TargetCov', 'TargetSeq']):
            self.steps.extend([self.targetcov, self.targqc_summary])
        # if Steps.contains(cnf.steps, 'AbnormalCovReport'):
        #    self.steps.append(self.abnormal_regions)

        if Steps.contains(cnf.steps, 'Seq2C'):
            self.steps.extend([self.seq2c])

        # if Steps.contains(cnf.steps, 'ClinicalReport') or \
        #         Steps.contains(cnf.steps, 'ClinicalReports') or \
        #         Steps.contains(cnf.steps, source.clinreport_name):
        self.steps.extend([self.clin_report])

        if Steps.contains(cnf.steps, 'Summary'):
            self.steps.extend([self.targqc_summary])

        # fastqc summary and clinical report -- special case (turn on if user uses default steps)
        if set(defaults['steps']) == set(cnf.steps):
            self.steps.extend([self.clin_report])

        from sys import platform as _platform
        if 'linux' in _platform:
            self.steps.append(self.bw_converting)

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
            '-o \'{output_dir}\' -s \'{sample}\' -c {caller} --qc ' +
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

        info()
        info('Checking BED files')
        self.is_wgs = self.cnf.is_wgs or self.bcbio_structure.is_wgs
        target_bed, exons_bed, exons_no_genes_bed, genes_fpath, seq2c_bed = self.prep_bed()

        if target_bed:
            ready_target_bed = join(self.bcbio_structure.date_dirpath, 'target.bed')
            try:
                shutil.copy(target_bed, ready_target_bed)
            except OSError:
                err(traceback.format_exc())
            target_bed = ready_target_bed
            cnf.bed = target_bed

        varfilter_paramline = params_for_one_sample + (' ' +
            '-o {output_dir} --output-file {output_file} -s {sample} -c {caller} --vcf {vcf} {vcf2txt_cmdl} --qc ' +
            '--work-dir ' + join(cnf.work_dir, BCBioStructure.varfilter_name) + '_{sample}_{caller} ')

        self.varfilter = Step(cnf, run_id,
            name=BCBioStructure.varfilter_name, short_name='vf',
            interpreter='python',
            script=join('scripts', 'post', 'varfilter.py'),
            dir_name=BCBioStructure.varfilter_dir,
            log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', BCBioStructure.varfilter_name + '-{caller}.log'),
            paramln=varfilter_paramline,
        )

        self.vcf2txt_single = Step(cnf, run_id,
            name='vcf2txt_single', short_name='vcf2txt_single',
            interpreter='perl',
            script='vcf2txt',
            dir_name=BCBioStructure.var_dir,
            log_fpath_template=join(self.bcbio_structure.log_dirpath, 'vcf2txt-single-{caller}.log'),
            paramln='{paramln}',
        )
        self.vcf2txt_paired = Step(cnf, run_id,
            name='vcf2txt_paired', short_name='vcf2txt_paired',
            interpreter='perl',
            script='vcf2txt',
            dir_name=BCBioStructure.var_dir,
            log_fpath_template=join(self.bcbio_structure.log_dirpath, 'vcf2txt-paired-{caller}.log'),
            paramln='{paramln}',
        )

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

        targetcov_params = params_for_one_sample + ' --bam \'{bam}\' -o \'{output_dir}\' ' \
            '-s \'{sample}\' --work-dir \'' + join(cnf.work_dir, BCBioStructure.targqc_name) + '_{sample}\' '
        if exons_bed:
            targetcov_params += '--exons ' + exons_bed + ' '
        if exons_no_genes_bed:
            targetcov_params += '--exons-no-genes ' + exons_no_genes_bed + ' '
        if target_bed:
            targetcov_params += '--bed ' + target_bed + ' '
        if self.bcbio_structure.bed:
            targetcov_params += '--original-bed ' + self.bcbio_structure.bed + ' '
        if genes_fpath:
            targetcov_params += '--genes ' + genes_fpath + ' '
        if cnf.no_dedup:
            targetcov_params += '--no-dedup '

        targetcov_params += '--no-prep-bed '
        if cnf.steps and 'AbnormalCovReport' in cnf.steps:
            targetcov_params += '--extended '
        self.targetcov = Step(cnf, run_id,
            name=BCBioStructure.targqc_name, short_name='tc',
            interpreter='python',
            script=join('scripts', 'post', 'targetcov.py'),
            dir_name=BCBioStructure.targqc_dir,
            log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', BCBioStructure.targqc_name + '.log'),
            paramln=targetcov_params
        )
        abnormal_regions_cmdl = summaries_cmdline_params + ' --mutations {mutations_fpath} ' + self.final_dir
        if target_bed:
            abnormal_regions_cmdl += ' --bed ' + target_bed
        self.abnormal_regions = Step(cnf, run_id,
            name='AbnormalCovReport', short_name='acr',
            interpreter='python',
            script=join('scripts', 'post', 'abnormal_regions.py'),
            dir_name=BCBioStructure.targqc_dir,
            log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', 'abnormalRegionsReport.log'),
            paramln=abnormal_regions_cmdl
        )
        self.ngscat = Step(cnf, run_id,
            name=BCBioStructure.ngscat_name, short_name='nc',
            interpreter='python',
            script=join('scripts', 'post', 'ngscat.py'),
            dir_name=BCBioStructure.ngscat_dir,
            log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', BCBioStructure.ngscat_name + '.log'),
            paramln=params_for_one_sample + ' --bam \'{bam}\' --bed ' + str(target_bed) + ' -o \'{output_dir}\' -s \'{sample}\' '
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
        # self.varqc_summary = Step(cnf, run_id,
        #     name=BCBioStructure.varqc_name + '_summary', short_name='vqs',
        #     interpreter='python',
        #     script=join('scripts', 'post_bcbio', 'varqc_summary.py'),
        #     dir_name=BCBioStructure.varqc_summary_dir,
        #     log_fpath_template=join(self.bcbio_structure.log_dirpath, BCBioStructure.varqc_name + '_summary.log'),
        #     paramln=summaries_cmdline_params + ' ' + self.final_dir +
        #         ' --varqc-name ' + BCBioStructure.varqc_name +
        #         ' --varqc-dir ' + BCBioStructure.varqc_summary_dir
        # )
        # self.varqc_after_summary = Step(cnf, run_id,
        #     name=BCBioStructure.varqc_after_name + '_summary', short_name='vqas',
        #     interpreter='python',
        #     script=join('scripts', 'post_bcbio', 'varqc_summary.py'),
        #     dir_name=BCBioStructure.varqc_after_summary_dir,
        #     log_fpath_template=join(self.bcbio_structure.log_dirpath, BCBioStructure.varqc_after_name + '_summary.log'),
        #     paramln=summaries_cmdline_params + ' ' + self.final_dir +
        #         ' --varqc-name ' + BCBioStructure.varqc_after_name +
        #         ' --varqc-dir ' + BCBioStructure.varqc_after_summary_dir
        # )

        clinreport_paramline = (params_for_one_sample +
           ' -s {sample}' +
           ' {targqc_cmdl}' +
           ' {mutations_cmdl}' +
           ' {varqc_cmdl}' +
           ' {match_cmdl}' +
           ' {seq2c_cmdl}' +
           ' {sv_cmdl}' +
           ' {sv_vcf_cmdl}' +
           ' {targqc_summary_cmdl}' +
           ' --target-type ' + self.bcbio_structure.target_type +
          (' --bed ' + target_bed if target_bed else '') +
          (' --jira ' + self.cnf.jira if self.cnf.jira else '') +
           ' -o {output_dir} ' +
           ' --project-level-report {project_report_path}')
        self.clin_report = Step(cnf, run_id,
            name=source.clinreport_name, short_name='clin',
            interpreter='python',
            script=join('scripts', 'post', 'clinical_report.py'),
            dir_name=source.clinreport_dir,
            log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', source.clinreport_name  + '.log'),
            paramln=clinreport_paramline + ' --work-dir ' + join(self.bcbio_structure.work_dir, '{sample}_' + source.clinreport_name)
        )
        self.clin_report_perl = Step(cnf, run_id,
            name=source.clinreport_name + '_perl', short_name='clin_perl',
            interpreter='python',
            script=join('scripts', 'post', 'clinical_report.py'),
            dir_name=source.clinreport_dir + '_perl',
            log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', source.clinreport_name  + '_perl.log'),
            paramln=clinreport_paramline + ' --work-dir ' + join(self.bcbio_structure.work_dir, '{sample}_' + source.clinreport_name + '_perl')
        )

        self.mongo_loader = Step(cnf, run_id,
            name='MongoLoader', short_name='ml',
            interpreter='java',
            script='vcf_loader',
            dir_name='mongo_loader',
            log_fpath_template=join(self.bcbio_structure.log_dirpath, 'mongo_loader.log'),
            paramln='-module loader -project {project} -sample {sample} -path {path} -variantCaller {variantCaller}'
        )

        seq2c_cmdline = summaries_cmdline_params + ' ' + self.final_dir + ' --genome {genome}'
        seq2c_cmdline += ' -t ' + str(self.max_threads)
        seq2c_cmdline += ' --bed ' + seq2c_bed + ' --no-prep-bed '
        if self.is_wgs:
            seq2c_cmdline += ' --wgs '
        normal_snames = [b.normal.name for b in self.bcbio_structure.batches.values() if b.normal]
        if normal_snames or cnf.seq2c_controls:
            controls = (normal_snames or []) + (cnf.seq2c_controls.split(':') if cnf.seq2c_controls else [])
            seq2c_cmdline += ' -c ' + ':'.join(controls)
        if cnf.seq2c_opts:
            seq2c_cmdline += ' --seq2c-opts ' + cnf.seq2c_opts
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

        targqc_summary_cmdline = summaries_cmdline_params + ' ' + self.final_dir
        if target_bed:
            targqc_summary_cmdline += ' --bed ' + target_bed
        if exons_bed:
            targqc_summary_cmdline += ' --exons ' + exons_bed

        self.targqc_summary = Step(cnf, run_id,
            name=BCBioStructure.targqc_name + '_summary', short_name='targqc',
            interpreter='python',
            script=join('scripts', 'post_bcbio', 'targqc_summary.py'),
            log_fpath_template=join(self.bcbio_structure.log_dirpath, BCBioStructure.targqc_name + '_summary.log'),
            dir_name=BCBioStructure.targqc_summary_dir,
            paramln=targqc_summary_cmdline
        )
        # self.fastqc_summary = Step(cnf, run_id,
        #     name=BCBioStructure.fastqc_name, short_name='fastqc',
        #     interpreter='python',
        #     script=join('scripts', 'post_bcbio', 'fastqc_summary.py'),
        #     log_fpath_template=join(self.bcbio_structure.log_dirpath, BCBioStructure.fastqc_name + '_summary.log'),
        #     dir_name=BCBioStructure.fastqc_summary_dir,
        #     paramln=summaries_cmdline_params + ' ' + self.final_dir
        # )

        self.bw_converting = Step(cnf, run_id,
            name='bam_to_bigwig', short_name='bamtobw',
            interpreter='python',
            script=join('scripts', 'post', 'bam_to_bigwig.py'),
            log_fpath_template=join(self.bcbio_structure.log_dirpath, '{sample}', BCBioStructure.bigwig_name + '.log'),
            paramln=basic_params + ' --genome {genome}  -s \'{sample}\' --bam \'{bam}\''
               ' --work-dir ' + join(self.bcbio_structure.work_dir, '{sample}_' + BCBioStructure.bigwig_name)
        )

    def prep_bed(self):
        target_bed, exons_bed, genes_fpath = get_bed_targqc_inputs(self.cnf, self.bcbio_structure.bed)
        exons_no_genes_bed = None
        seq2c_bed = self.bcbio_structure.sv_bed

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

        return target_bed, exons_bed, exons_no_genes_bed, genes_fpath, seq2c_bed

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
                    log_out_fpath=None, wait_for_steps=None, threads=1, mem_m=None, **kwargs):
        job_name = step.job_name(sample_name, caller_suf)

        for job_id_to_wait in wait_for_steps or []:
            job_to_wait = next((j for j in self.jobs_running if j.job_id == job_id_to_wait), None)
            if not job_to_wait or job_to_wait.has_errored:
                if job_to_wait:
                    warn('Job ' + job_to_wait.job_id + ' has failed, and it required to run this job ' + job_name)
                    warn()

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
        queue = self.cnf.queue
        runner_script = self.qsub_runner
        bash = get_system_path(self.cnf, 'bash')
        extra_qsub_opts = ''
        if step.run_on_chara and is_us():
            extra_qsub_opts += '-l h="chara|rask" '

        mem_opts = ''
        if mem_m and not is_local():
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
            info(step.name + (' - ' + sample_name if sample_name else '') + (' - ' + caller_suf if caller_suf else ''))
            info(qsub_cmdline)
        else:
            print step.name,

        if isfile(done_marker_fpath): os.remove(done_marker_fpath)
        if isfile(error_marker_fpath): os.remove(error_marker_fpath)
        job = JobRunning(step, job_name, sample_name, caller_suf, log_err_fpath, qsub_cmdline,
                         done_marker_fpath, error_marker_fpath, threads=threads, not_wait=BCBioStructure.bigwig_name in step.name)
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

        try:
            targqc_wait_for_steps = []
            for sample in self.bcbio_structure.samples:
                if not (any(step in self.steps for step in
                            [self.targetcov,
                             # self.varqc,
                             self.varannotate])):
                    continue

                info(sample.name)
                # BAMS
                if self.targetcov in self.steps:
                    if not sample.bam or not verify_bam(sample.bam):
                        err('Cannot run coverage reports (targetcov, qualimap, ngscat) without BAM files.')
                        info('Target coverage for "' + sample.name + '"')
                    else:
                        self._submit_job(
                            self.targetcov, sample.name,
                            bam=sample.bam, sample=sample.name, genome=sample.genome,
                            caller_names='', vcfs='', threads=self.threads_per_sample, wait_for_steps=targqc_wait_for_steps,
                            mem_m=max(3000, get_qualimap_max_mem(sample.bam)) + 4100, run_on_chara=True)  # --java-mem-size (or 3Gb for picard) + 4Gb for JavaVM

                # Processing VCFs: QC, annotation
                for caller in self.bcbio_structure.variant_callers.values():
                    vcf_fpath = sample.vcf_by_callername.get(caller.name)
                    if not vcf_fpath:
                        if sample.phenotype != 'normal':
                            err('VCF does not exist: sample ' + sample.name + ', caller ' + caller.name + '.')
                    else:
                        self._process_vcf(sample, sample.bam, vcf_fpath, caller, threads=self.threads_per_sample)

                info('-' * 70)

            if self.varfilter in self.steps:
                info('Filtering')
                    # info('Cohort mode set, running vcf2txt in cohort mode')
                    # for c in self.bcbio_structure.variant_callers.values():
                    #     if self.varannotate not in self.steps:
                    #         is_err = False
                    #         for anno_vcf_fpath in c.single_anno_vcf_by_sample.values() + c.paired_anno_vcf_by_sample.values():
                    #             if not verify_vcf(anno_vcf_fpath):
                    #                 is_err = True
                    #         if is_err:
                    #             critical('Error: VarAnnotate is not in steps, and some annotated VCFs do not exist.')
                    #
                    #     if c.single_anno_vcf_by_sample:
                    #         self._submit_job(
                    #             self.vcf2txt_single,
                    #             paramln=make_vcf2txt_cmdl_params(self.cnf, c.single_anno_vcf_by_sample, sample.min_af) +
                    #                 ' > ' + c.single_vcf2txt_res_fpath,
                    #             caller_suf=c.name,
                    #             caller=c.name,
                    #             wait_for_steps=[
                    #                 self.varannotate.job_name(s.name, v.name)
                    #                 for v in self.bcbio_structure.variant_callers.values()
                    #                 for s in v.samples
                    #                 if self.varannotate in self.steps])
                    #             # mem_m=sum(getsize(c.single_anno_vcf_by_sample.values())) / 1024 / 1024 * 10 + 500)
                    #
                    #     if c.paired_anno_vcf_by_sample:
                    #         self._submit_job(
                    #             self.vcf2txt_paired,
                    #             paramln=make_vcf2txt_cmdl_params(self.cnf, c.paired_anno_vcf_by_sample, sample.min_af) +
                    #                 ' > ' + c.paired_vcf2txt_res_fpath,
                    #             caller_suf=c.name,
                    #             caller=c.name,
                    #             wait_for_steps=[
                    #                 self.varannotate.job_name(s.name, v.name)
                    #                 for v in self.bcbio_structure.variant_callers.values()
                    #                 for s in v.samples
                    #                 if self.varannotate in self.steps])
                    #             # mem_m=sum([getsize(v) for v in c.paired_anno_vcf_by_sample.values() if v]) / 1024 / 1024 * 10 + 500)

                info('Per-sample variant filtering')
                for sample in self.bcbio_structure.samples:
                    for caller in self.bcbio_structure.variant_callers.values():
                        if sample.vcf_by_callername.get(caller.name):
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

            # TargetSeq reports
            if self.abnormal_regions in self.steps:
                variant_caller = \
                    self.bcbio_structure.variant_callers.get('vardict') or \
                    self.bcbio_structure.variant_callers.get('vardict-java')
                if variant_caller:
                    vardict_txt_fname = source.mut_fname_template.format(caller_name=variant_caller.name)
                    vardict_txt_fpath = join(self.bcbio_structure.date_dirpath, vardict_txt_fname)
                    mutations_fpath = add_suffix(vardict_txt_fpath, source.mut_pass_suffix)

                    wait_for_steps = []
                    wait_for_steps += [self.varannotate.job_name(sample.name, caller=variant_caller.name) for sample in self.bcbio_structure.samples] if self.varannotate in self.steps else []
                    wait_for_steps += [self.varfilter.job_name(sample.name, caller=variant_caller.name) for sample in self.bcbio_structure.samples] if self.varfilter in self.steps else []
                    wait_for_steps += [self.targetcov.job_name(sample.name) for sample in self.bcbio_structure.samples] if self.targetcov in self.steps else []

                    self._submit_job(self.abnormal_regions, wait_for_steps=wait_for_steps, mutations_fpath=mutations_fpath)

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

            # if self.varqc_summary in self.steps:
            #     self._submit_job(
            #         self.varqc_summary,
            #         wait_for_steps=[
            #             self.varannotate.job_name(s.name, v.name)
            #             for v in self.bcbio_structure.variant_callers.values()
            #             for s in v.samples
            #             if self.varannotate in self.steps])

            # if self.varqc_after in self.steps:
            #     info('VarQC_postVarFilter:')
            #     for caller in self.bcbio_structure.variant_callers.values():
            #         for sample in caller.samples:
            #             raw_vcf_fpath = sample.find_raw_vcf_by_callername(caller.name)
            #             if not raw_vcf_fpath:
            #                 if sample.phenotype != 'normal':
            #                     err('Error: raw VCF does not exist: sample ' + sample.name + ', caller "' +
            #                         caller.name + '". Phenotype = ' + sample.phenotype + '.')
            #             else:
            #                 filt_vcf_fpath = sample.get_filt_vcf_fpath_by_callername(caller.name, gz=True)
            #                 if not self.varfilter and sample.phenotype != 'normal' and not verify_file(filt_vcf_fpath, silent=True):
            #                     err('Error: filtered VCF does not exist: sample ' + sample.name + ', caller "' +
            #                         caller.name + '". Phenotype = ' + sample.phenotype + '.' +
            #                         ' Note that you need to run VarFilter first, and this step is not in config.')
            #                 else:
            #                     self._submit_job(
            #                         self.varqc_after, sample.name, caller_suf=caller.name, threads=self.threads_per_sample,
            #                         wait_for_steps=([self.varfilter.job_name(caller=caller.name)] if self.varfilter in self.steps else []),
            #                         vcf=filt_vcf_fpath, sample=sample.name, caller=caller.name, genome=sample.genome)

            # if self.varqc_after_summary in self.steps:
            #     self._submit_job(
            #         self.varqc_after_summary,
            #         wait_for_steps=[
            #             self.varfilter.job_name(s.name, v.name)
            #             for v in self.bcbio_structure.variant_callers.values()
            #             for s in v.samples
            #             if self.varfilter in self.steps])

            if self.clin_report in self.steps:
                clinical_report_caller = \
                    self.bcbio_structure.variant_callers.get('vardict') or \
                    self.bcbio_structure.variant_callers.get('vardict-java')
                if not clinical_report_caller:
                    err('No vardict or vardict-java in the variants callers: ' + ', '.join(self.bcbio_structure.variant_callers.keys()))
                else:
                    seq2c_cmdl = ''
                    if self.seq2c in self.steps or verify_file(self.bcbio_structure.seq2c_fpath, silent=True):
                        seq2c_cmdl = ' --seq2c ' + self.bcbio_structure.seq2c_fpath

                    for sample in self.bcbio_structure.samples:
                        varqc_cmdl = ''
                        mutation_cmdl = ''
                        match_cmdl = ''
                        targqc_cmdl = ''
                        targqc_summary_cmdl = ''
                        sv_cmdl = ''
                        sv_vcf_cmdl = ''

                        wait_for_steps = []
                        wait_for_steps += [self.targetcov.job_name(sample.name)] if self.targetcov in self.steps else []
                        wait_for_steps += [self.seq2c.job_name()] if self.seq2c in self.steps else []

                        if not sample.phenotype or sample.phenotype != 'normal':
                            match_cmdl = ' --match ' + sample.normal_match.name if sample.normal_match else ''
                            varqc_cmdl = ' --varqc ' + sample.get_varqc_fpath_by_callername(clinical_report_caller.name, ext='.json')
                            wait_for_steps += [self.varannotate.job_name(sample.name, caller=clinical_report_caller.name)] if self.varannotate in self.steps else []
                            wait_for_steps += [self.varfilter.job_name(sample.name, caller=clinical_report_caller.name)] if self.varfilter in self.steps else []

                            mut_fpath = add_suffix(sample.get_vcf2txt_by_callername(clinical_report_caller.name), source.mut_pass_suffix)
                            if self.varfilter in self.steps or verify_file(mut_fpath):
                                mutation_cmdl = ' --mutations ' + mut_fpath

                        targqc_dirpath = join(self.final_dir, sample.name, BCBioStructure.targqc_dir)
                        if self.targetcov in self.steps or verify_dir(targqc_dirpath):
                            targqc_cmdl = ' --targqc-dir ' + join(self.final_dir, sample.name, BCBioStructure.targqc_dir)

                            targqc_summary_cmdl = ''
                            if self.targqc_summary in self.steps or verify_file(self.bcbio_structure.targqc_summary_fpath, silent=True):
                                targqc_summary_cmdl += ' --targqc-html ' + self.bcbio_structure.targqc_summary_fpath

                        sv_fpath = sample.find_sv_fpath()
                        if sv_fpath:
                            sv_cmdl = ' --sv ' + sv_fpath

                        sample_dirpath = join(self.bcbio_structure.final_dirpath, sample.name)
                        sample_cnv_dirpath = join(sample_dirpath, BCBioStructure.cnv_dir)

                        if exists(sample_cnv_dirpath):
                            for fname in os.listdir(sample_cnv_dirpath):
                                if '-manta' in fname and fname.endswith('.vcf.gz'):
                                    sv_vcf_cmdl = ' --sv-vcf ' + join(sample_cnv_dirpath, fname)

                        self._submit_job(
                            self.clin_report,
                            sample.name,
                            sample=sample.name, genome=sample.genome,
                            match_cmdl=match_cmdl, mutations_cmdl=mutation_cmdl,
                            varqc_cmdl=varqc_cmdl, targqc_cmdl=targqc_cmdl,
                            seq2c_cmdl=seq2c_cmdl, sv_cmdl=sv_cmdl, sv_vcf_cmdl=sv_vcf_cmdl,
                            targqc_summary_cmdl=targqc_summary_cmdl,
                            project_report_path=self.bcbio_structure.project_report_html_fpath,
                            wait_for_steps=wait_for_steps,
                            threads=self.threads_per_sample)

            if self.bw_converting in self.steps:
                for sample in self.bcbio_structure.samples:
                    self._submit_job(self.bw_converting, sample.name,
                        sample=sample.name, genome=sample.genome, bam=sample.bam,
                        # wait_for_steps=[self.targetcov.job_name(sample.name)] if self.targetcov in self.steps else [],
                        mem_m=getsize(sample.bam) * 1.1 / 1024 / 1024 + 500)
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
            info()

            if any(s.fastqc_html_fpath and isfile(s.fastqc_html_fpath) for s in self.bcbio_structure.samples):
                final_summary_report_fpath = join(self.bcbio_structure.date_dirpath, BCBioStructure.fastqc_summary_dir, source.fastqc_name + '.html')
                safe_mkdir(dirname(final_summary_report_fpath))
                write_fastqc_combo_report(self.cnf, final_summary_report_fpath, self.bcbio_structure.samples)

            if self.varannotate in self.steps:
                info('Making varQC summary reports')
                self._varqc_summary(BCBioStructure.varqc_dir, BCBioStructure.varqc_summary_dir, BCBioStructure.varqc_name)
                info()

            if self.varfilter in self.steps:
                finish_filtering_for_bcbio(self.cnf, self.bcbio_structure, self.bcbio_structure.variant_callers.values())
                info('Making varQC-post-filtering summary reports')
                self._varqc_summary(BCBioStructure.varqc_after_dir, BCBioStructure.varqc_after_summary_dir, BCBioStructure.varqc_after_name)
                info()

            if is_us() or is_uk():
                info('Exposing to jBrowse')
                add_project_files_to_jbrowse(self.cnf, self.bcbio_structure)
                info()

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
        except KeyboardInterrupt:
            info('Interrupted.')
        except SystemExit:
            info('Interrupted.')
        finally:
            info('Deleting running jobs...')
            del_jobs(self.cnf, self.jobs_running)
            info('Changing permissions...')
            if isdir(self.bcbio_structure.final_dirpath):
                change_permissions(self.bcbio_structure.final_dirpath)
            if isdir(self.bcbio_structure.work_dir):
                change_permissions(self.bcbio_structure.work_dir)
            if isdir(join(self.bcbio_structure.work_dir, '..', 'config')):
                change_permissions(join(self.bcbio_structure.work_dir, '..', 'config'))
            info()
            info('Done post-processing.')


    def _varqc_summary(self, sample_qc_path, summary_qc_path, varqc_level_name):
        jsons_by_sample_by_caller = defaultdict(dict)
        htmls_by_sample_by_caller = defaultdict(dict)
        for vc in self.bcbio_structure.variant_callers.values():
            fpath = vc.find_fpaths_by_sample(sample_qc_path, varqc_level_name, 'json', self.bcbio_structure.final_dirpath)
            if fpath:
                jsons_by_sample_by_caller[vc.name] = fpath
            info('Found JSONs: ' + str(', '.join(k + ': ' + str(v) for k, v in jsons_by_sample_by_caller[vc.name].items())))
            fpath = vc.find_fpaths_by_sample(sample_qc_path, varqc_level_name, 'html', self.bcbio_structure.final_dirpath)
            if fpath:
                htmls_by_sample_by_caller[vc.name] = fpath
            info('Found HTMLs: ' + str(', '.join(k + ': ' + str(v) for k, v in htmls_by_sample_by_caller[vc.name].items())))

        info()
        if jsons_by_sample_by_caller and htmls_by_sample_by_caller:
            info('Writing summary reports...')
            summarize_qc.make_summary_reports(self.cnf, 1, join(self.bcbio_structure.date_dirpath, summary_qc_path),
                 self.bcbio_structure.variant_callers.values(), self.bcbio_structure.samples,
                 jsons_by_sample_by_caller, htmls_by_sample_by_caller,
                 varqc_name=BCBioStructure.varqc_name, caption='Variant QC')
        else:
            err('Not JSON and HTML found, cannot generate summary reports.')


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


    def _process_vcf(self, sample, bam_fpath, vcf_fpath, caller, threads,
                     steps=None, job_names_to_wait=None):
        steps = steps or self.steps

        bam_cmdline = '--bam ' + bam_fpath if bam_fpath else ''

        normal_match_cmdline = ''
        if sample.normal_match:
            normal_match_cmdline = ' --match-normal-sample-name ' + sample.normal_match.name + ' '

        if self.varannotate in steps:
            self._submit_job(
                self.varannotate, sample.name, caller_suf=caller.name, vcf=vcf_fpath, threads=threads,
                bam_cmdline=bam_cmdline, sample=sample.name, caller=caller.name,
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
        try:
            safe_mkdir(cnv_summary_dirpath)
        except OSError:
            pass

        for sample in self.bcbio_structure.samples:
            sample_dirpath = join(self.bcbio_structure.final_dirpath, sample.name)
            sample_cnv_dirpath = join(sample_dirpath, BCBioStructure.cnv_dir)

            for cnv_caller in ['-cn_mops', '-cnvkit', '-lumpy', '-manta', '-wham', '-sv-prioritize', '-metasv']:
                for fname in os.listdir(sample_dirpath):
                    if cnv_caller in fname:
                        # Copy to <sample>/cnv
                        safe_mkdir(sample_cnv_dirpath)
                        try:
                            os.rename(join(sample_dirpath, fname), join(sample_cnv_dirpath, fname))
                        except OSError:
                            err(format_exc())
                            info()

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
    txt += 'URL: ' + convert_path_to_url(bs.final_dirpath) + '\n'
    txt += 'Report: ' + (html_report_url or bs.project_report_html_fpath) + '\n'
    if jira_url:
        txt += 'Jira: ' + jira_url
    send_email(txt, subj)


def change_permissions(path):
    try:
        os.system('chmod -R g+w ' + path)
    except:
        pass
