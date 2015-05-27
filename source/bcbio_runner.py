import getpass
import os
import shutil
import sys
import hashlib
import base64
from os.path import join, dirname, abspath, expanduser, basename, pardir, isfile, isdir, exists, islink, relpath
import datetime
from time import sleep
from source.bcbio_structure import BCBioStructure
from source.calling_process import call
from source.file_utils import verify_dir, verify_file, add_suffix, symlink_plus, remove_quotes
from source.tools_from_cnf import get_system_path

from source.file_utils import file_exists, safe_mkdir
from source.logger import info, err, critical, send_email
from source.ngscat.bed_file import verify_bam


class Step:
    def __init__(self, cnf, run_id, name, script, dir_name=None,
                 interpreter=None, short_name=None, paramln='', env_vars=None):
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

    def job_name(self, sample=None, caller=None):
        return self.short_name.upper() + '_' + self.run_id_ + \
               ('_' + sample if sample else '') + \
               ('_' + caller if caller else '')

    def __repr__(self):
        return self.name


class Steps(list):
    def __init__(self):
        super(Steps, self).__init__()

    def add_step(self, step):
        self.append(step)

    def extend(self, iterable):
        for step in iterable:
            self.add_step(step)


class JobRunning:
    def __init__(self, step, sample_name, caller_suf, log_fpath, qsub_cmdline, done_marker):
        self.step = step
        self.sample_name = sample_name
        self.caller_suf = caller_suf
        self.log_fpath = log_fpath
        self.qsub_cmdline = qsub_cmdline
        self.done_marker = done_marker
        self.repr = step.name
        if sample_name:
            self.repr += ' for ' + sample_name
        if caller_suf:
            self.repr += ', ' + caller_suf
        self.is_done = False


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
        self.run_id = self.__generate_run_id(self.final_dir, bcbio_structure.project_name, user_prid, timestamp)
        info('User PRID: ' + user_prid + ', run_id: ' + self.run_id)

        self.qsub_runner = abspath(expanduser(cnf.qsub_runner))

        self.max_threads = self.cnf.threads
        total_samples_num = len(self.bcbio_structure.samples)
        total_callers_num = total_samples_num * len(self.bcbio_structure.variant_callers)
        self.summary_threads = min(self.max_threads, total_samples_num)
        self.threads_per_sample = 1  # max(self.max_threads / total_samples_num, 1)

        self.steps = Steps()
        self._set_up_steps(cnf, self.run_id)

        self.jobs_running = []

        normalize = lambda name: name.lower().replace('_', '').replace('-', '')
        contains = lambda x, xs: normalize(x) in [normalize(y) for y in (xs or [])]

        self.steps.extend([
            self.varqc_summary,
            self.varqc_after_summary,
            self.fastqc_summary,
            self.targqc_summary,
            self.combined_report,
        ])

        self.steps.extend(
            [s for s in [
                self.varqc,
                self.varannotate,
                self.varfilter_all,
                self.mongo_loader,
                self.varqc_after,
                self.targetcov,
                self.abnormal_regions,
                self.seq2c,
                self.ngscat,
                self.qualimap,
            ] if contains(s.name, cnf.steps)])

        if cnf.steps and 'TargetCov' in cnf.steps:
            self.steps.append(self.targetcov)

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

        self._symlink_cnv()

    def __generate_run_id(self, final_dir, project_name, prid='', timestamp=''):
        hasher = hashlib.sha1(final_dir + prid + timestamp)
        path_hash = base64.urlsafe_b64encode(hasher.digest()[0:4])[:-1]
        return path_hash + '_' + project_name

    def _set_up_steps(self, cnf, run_id):
        basic_params = \
            ' --sys-cnf ' + self.cnf.sys_cnf + \
            ' --run-cnf ' + self.cnf.run_cnf

        summaries_cmdline_params = \
            basic_params + \
            ' -t ' + str(self.summary_threads) + \
           (' --reuse ' if self.cnf.reuse_intermediate else '') + \
            ' --log-dir {log_dirpath}'

        # Params for those who doesn't call bcbio_structure
        params_for_one_sample = \
            basic_params + \
            ' -t ' + str(self.threads_per_sample) + \
           (' --reuse ' if self.cnf.reuse_intermediate else '') + \
            ' --log-dir {log_dirpath}' + \
            ' --genome {genome}' + \
            ' --project-name ' + self.bcbio_structure.project_name + ' '

        if cnf.email:
            summaries_cmdline_params += ' --email ' + remove_quotes(self.cnf.email) + ' '
            params_for_one_sample += ' --email ' + remove_quotes(self.cnf.email) + ' '

        if cnf.bed:
            summaries_cmdline_params += ' --bed ' + cnf.bed

        anno_paramline = params_for_one_sample + ('' +
            ' --vcf \'{vcf}\' {bam_cmdline} {normal_match_cmdline} ' +
            '-o \'{output_dir}\' -s \'{sample}\' -c {caller} ' +
            '--work-dir \'' + join(cnf.work_dir, BCBioStructure.varannotate_name) + '_{sample}_{caller}\' ')

        self.varannotate = Step(cnf, run_id,
            name=BCBioStructure.varannotate_name, short_name='va',
            interpreter='python',
            script=join('sub_scripts', 'varannotate.py'),
            dir_name=BCBioStructure.varannotate_dir,
            paramln=anno_paramline,
        )
        self.varqc = Step(cnf, run_id,
            name=BCBioStructure.varqc_name, short_name='vq',
            interpreter='python',
            script=join('sub_scripts', 'varqc.py'),
            dir_name=BCBioStructure.varqc_dir,
            paramln=params_for_one_sample + ' --vcf \'{vcf}\' -o \'{output_dir}\' -s \'{sample}\' -c {caller} '
                    '--work-dir \'' + join(cnf.work_dir, BCBioStructure.varqc_name) + '_{sample}_{caller}\''
        )
        self.varqc_after = Step(cnf, run_id,
            name=BCBioStructure.varqc_after_name, short_name='vqa',
            interpreter='python',
            script=join('sub_scripts', 'varqc.py'),
            dir_name=BCBioStructure.varqc_after_dir,
            paramln=params_for_one_sample + ' --vcf \'{vcf}\' -o \'{output_dir}\' -s \'{sample}\' -c {caller} '
                    '--work-dir \'' + join(cnf.work_dir, BCBioStructure.varqc_after_name) + '_{sample}_{caller}\' ' +
                    '--proc-name ' + BCBioStructure.varqc_after_name
        )
        targetcov_params = params_for_one_sample + ' --bam \'{bam}\' --bed \'{bed}\' -o \'{output_dir}\' ' \
            '-s \'{sample}\' --work-dir \'' + join(cnf.work_dir, BCBioStructure.targetseq_name) + '_{sample}\' '
        if cnf.exons:
            targetcov_params += '--exons {cnf.exons} '
        if cnf.reannotate:
            targetcov_params += '--reannotate '
        if cnf.count_dups:
            targetcov_params += '--count-dups'

        self.targetcov = Step(cnf, run_id,
            name=BCBioStructure.targetseq_name, short_name='tc',
            interpreter='python',
            script=join('sub_scripts', 'targetcov.py'),
            dir_name=BCBioStructure.targetseq_dir,
            paramln=targetcov_params,
        )
        self.abnormal_regions = Step(cnf, run_id,
            name='AbnormalCovReport', short_name='acr',
            interpreter='python',
            script=join('sub_scripts', 'abnormal_regions.py'),
            dir_name=BCBioStructure.targetseq_dir,
            paramln=params_for_one_sample + ' -o \'{output_dir}\' {caller_names} {vcfs} '
                    '-s \'{sample}\' --work-dir \'' + join(cnf.work_dir, BCBioStructure.targetseq_name) + '_{sample}\' '
        )
        self.ngscat = Step(cnf, run_id,
            name=BCBioStructure.ngscat_name, short_name='nc',
            interpreter='python',
            script=join('sub_scripts', 'ngscat.py'),
            dir_name=BCBioStructure.ngscat_dir,
            paramln=params_for_one_sample + ' --bam \'{bam}\' --bed \'{bed}\' -o \'{output_dir}\' -s \'{sample}\' '
                    '--saturation y --work-dir \'' + join(cnf.work_dir, BCBioStructure.ngscat_name) + '_{sample}\''
        )
        self.qualimap = Step(cnf, run_id,
            name=BCBioStructure.qualimap_name, short_name='qm',
            interpreter='python',
            script=join('sub_scripts', 'qualimap.py'),
            dir_name=BCBioStructure.qualimap_dir,
            paramln=params_for_one_sample + ' --bam {bam} {bed} -o {output_dir}',
        )
        #############
        # Summaries #
        self.varqc_summary = Step(cnf, run_id,
            name=BCBioStructure.varqc_name + '_summary', short_name='vqs',
            interpreter='python',
            script=join('sub_scripts', 'varqc_summary.py'),
            dir_name=BCBioStructure.varqc_summary_dir,
            paramln=summaries_cmdline_params + ' ' + self.final_dir
        )
        self.varqc_after_summary = Step(cnf, run_id,
            name=BCBioStructure.varqc_after_name + '_summary', short_name='vqas',
            interpreter='python',
            script=join('sub_scripts', 'varqc_summary.py'),
            dir_name=BCBioStructure.varqc_after_summary_dir,
            paramln=summaries_cmdline_params + ' ' + self.final_dir +
                    ' --name ' + BCBioStructure.varqc_after_name +
                    ' --dir ' + BCBioStructure.varqc_after_dir
        )
        varfilter_paramline = summaries_cmdline_params + ' ' + self.final_dir
        if cnf.datahub_path:
            varfilter_paramline += ' --datahub-path ' + cnf.datahub_path
        if cnf.min_freq is not None:
            varfilter_paramline += ' --freq ' + str(cnf.min_freq)

        self.varfilter_all = Step(cnf, run_id,
            name=BCBioStructure.varfilter_name, short_name='vfs',
            interpreter='python',
            script=join('sub_scripts', 'varfilter.py'),
            dir_name=BCBioStructure.varfilter_dir,
            paramln=varfilter_paramline
        )

        self.mongo_loader = Step(cnf, run_id,
            name='MongoLoader', short_name='ml',
            interpreter='java',
            script='vcf_loader',
            dir_name='mongo_loader',
            paramln='-module loader -project {project} -sample {sample} -path {path} -variantCaller {variantCaller}'
        )
        seq2c_cmdline = summaries_cmdline_params + ' ' + self.final_dir + ' --genome {genome} '
        if cnf.controls:
            seq2c_cmdline += ' -c ' + cnf.controls
        if cnf.seq2c_opts:
            seq2c_cmdline += ' --seq2c_opts ' + cnf.seq2c_opts
        self.seq2c = Step(cnf, run_id,
            name=BCBioStructure.seq2c_name, short_name='seq2c',
            interpreter='python',
            script=join('sub_scripts', 'seq2c.py'),
            dir_name=BCBioStructure.cnv_summary_dir,
            paramln=seq2c_cmdline
        )
        self.targqc_summary = Step(cnf, run_id,
            name=BCBioStructure.targqc_name, short_name='targqc',
            interpreter='python',
            script=join('sub_scripts', 'targqc_summary.py'),
            dir_name=BCBioStructure.targqc_summary_dir,
            paramln=summaries_cmdline_params + ' ' + self.final_dir
        )
        self.fastqc_summary = Step(cnf, run_id,
            name=BCBioStructure.fastqc_name, short_name='fastqc',
            interpreter='python',
            script=join('sub_scripts', 'fastqc_summary.py'),
            dir_name=BCBioStructure.fastqc_summary_dir,
            paramln=summaries_cmdline_params + ' ' + self.final_dir
        )
        project_level_report_cmdline = summaries_cmdline_params + ' ' + self.final_dir
        if cnf.jira:
            project_level_report_cmdline += ' --jira ' + cnf.jira
        self.combined_report = Step(cnf, run_id,
            name='ProjectLevelReport', short_name='cr',
            interpreter='python',
            script=join('sub_scripts', 'combined_report.py'),
            dir_name=self.bcbio_structure.date_dirpath,
            paramln=project_level_report_cmdline
        )


    def step_log_marker_and_output_paths(self, step, sample_name, caller=None):
        if sample_name:
            base_output_dirpath = abspath(join(self.final_dir, sample_name))
        else:
            base_output_dirpath = abspath(self.bcbio_structure.date_dirpath)

        output_dirpath = join(base_output_dirpath, step.dir_name) if step.dir_name else ''

        log_fpath = join(self.bcbio_structure.log_dirpath,
             (step.name + ('_' + sample_name if sample_name else '') +
                          ('_' + caller if caller else '')) + '.log')

        log_dirpath = join(self.bcbio_structure.log_dirpath, step.name)
        safe_mkdir(log_dirpath)

        done_markers_dirpath = join(self.bcbio_structure.work_dir, 'done_markers')
        marker_fpath = join(done_markers_dirpath,
             (step.name + '_' + step.run_id_ +
              ('_' + sample_name if sample_name else '') +
              ('_' + caller if caller else '')) + '.done')
        safe_mkdir(done_markers_dirpath)

        return output_dirpath, log_fpath, log_dirpath, marker_fpath


    def _submit_job(self, step, sample_name='', caller_suf=None, create_dir=True,
                    log_out_fpath=None, wait_for_steps=None, threads=1, **kwargs):

        output_dirpath, log_err_fpath, log_dirpath, marker_fpath = \
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
        params = dict({'output_dir': output_dirpath, 'log_dirpath': log_dirpath}.items() +
                      self.__dict__.items() + kwargs.items())
        cmdline = tool_cmdline + ' ' + step.param_line.format(**params) + ' --done-marker ' + marker_fpath

        hold_jid_line = '-hold_jid ' + ','.join(wait_for_steps or ['_'])
        job_name = step.job_name(sample_name, caller_suf)
        qsub = get_system_path(self.cnf, 'qsub')
        if threads > 1:
            threads += 1
        threads = str(threads)
        queue = self.cnf.queue
        runner_script = self.qsub_runner
        bash = get_system_path(self.cnf, 'bash')
        qsub_cmdline = (
            '{qsub} -pe smp {threads} -S {bash} -q {queue} '
            '-j n -o {log_err_fpath} -e {log_err_fpath} {hold_jid_line} '
            '-N {job_name} {runner_script} {marker_fpath} "{cmdline}"'.format(**locals()))

        if self.cnf.verbose:
            info(step.name)
            info(qsub_cmdline)
        else:
            print step.name,

        if isfile(marker_fpath):
            os.remove(marker_fpath)
        job = JobRunning(step, sample_name, caller_suf, log_err_fpath, qsub_cmdline, marker_fpath)
        self.jobs_running.append(job)
        call(self.cnf, qsub_cmdline, silent=True, env_vars=step.env_vars)

        if self.cnf.verbose: info()
        return output_dirpath


    def _qualimap_bed(self, bed_fpath):
        if self.qualimap in self.steps and bed_fpath:
            qualimap_bed_fpath = join(self.cnf.work_dir, 'tmp_qualimap.bed')

            fix_bed_for_qualimap(bed_fpath, qualimap_bed_fpath)

            return qualimap_bed_fpath
        else:
            return bed_fpath


    def post_jobs(self):
        callers = self.bcbio_structure.variant_callers.values()

        if self.qualimap in self.steps:
            bed_by_sample = dict((s.name, s.bed) for s in self.bcbio_structure.samples if s.bed)
            beds = set(bed_by_sample.values())
            samples_by_bed = dict((b, (s for s in self.bcbio_structure.samples if s.bed and s.bed == b)) for b in beds)
            for bed, samples in samples_by_bed.items():
                qualimap_bed = self._qualimap_bed(bed)
                for s in samples:
                    s.bed = bed
                    s.qualimap_bed = qualimap_bed

        for sample in self.bcbio_structure.samples:
            if not (any(step in self.steps for step in
                        [self.targetcov,
                         self.seq2c,
                         self.qualimap,
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
                   self.qualimap,
                   self.ngscat]):
                if not sample.bam or not verify_bam(sample.bam):
                    err('Cannot run coverage reports (targetcov, qualimap, ngscat) without BAM files.')
                else:
                    # TargetCov reports
                    if self.targetcov in self.steps:
                        info('Target coverage for "' + sample.name + '"')
                        self._submit_job(
                            self.targetcov, sample.name,
                            bam=sample.bam, bed=sample.bed or self.cnf.genomes[sample.genome].exons, sample=sample.name, genome=sample.genome,
                            caller_names='', vcfs='', threads=self.threads_per_sample)

                    # ngsCAT reports
                    if (self.ngscat in self.steps) and (not sample.bed or not verify_file(sample.bed)):
                        err('Warning: no BED file, assuming WGS, thus skipping ngsCAT reports.')
                    else:
                        if self.ngscat in self.steps:
                            self._submit_job(
                                self.ngscat, sample.name, bam=sample.bam, bed=sample.bed or self.cnf.genomes[sample.genome].exons,
                                sample=sample.name, genome=sample.genome, threads=self.threads_per_sample)

                    # Qualimap
                    if self.qualimap in self.steps:
                        qualimap_gff = ''
                        if sample.bed:
                            qualimap_gff = ' --bed ' + sample.qualimap_bed + ' '
                        self._submit_job(
                            self.qualimap, sample.name, bam=sample.bam, sample=sample.name,
                            genome=sample.genome, bed=qualimap_gff, threads=self.threads_per_sample)

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

        # if self.vardict_steps:
        #     self._sumbit_vardict(self.bcbio_structure.batches)

        if self.varqc_summary in self.steps:
            self._submit_job(
                self.varqc_summary,
                wait_for_steps=[
                    self.varqc.job_name(s.name, v.name)
                    for v in self.bcbio_structure.variant_callers.values()
                    for s in v.samples
                    if self.varqc in self.steps])

        if self.varfilter_all in self.steps:
            self._submit_job(
                self.varfilter_all,
                wait_for_steps=[
                    self.varannotate.job_name(s.name, v.name)
                    for v in self.bcbio_structure.variant_callers.values()
                    for s in v.samples
                    if self.varannotate in self.steps],
                create_dir=False,
                threads=self.summary_threads)

        # TargetSeq reports
        if self.abnormal_regions in self.steps:
            for sample in self.bcbio_structure.samples:
                if not self.cnf.verbose:
                    info(ending='')

                if not sample.bed or not verify_file(sample.bed):
                    err('Warning: no BED file, assuming WGS, thus running targetSeq reports '
                        'only to generate Seq2C reports.')
                    continue

                callers_and_filtered_vcfs = [(c, f) for c, f in ((c.name, c.get_filt_vcf_by_sample().get(sample.name)) for c in callers) if f]
                if callers_and_filtered_vcfs:
                    caller_names, filtered_vcfs = zip(*callers_and_filtered_vcfs)
                else:
                    caller_names, filtered_vcfs = [], []

                wait_for_steps = []
                if self.varfilter_all in self.steps:
                    wait_for_steps.extend([self.varfilter_all.job_name()])
                if self.targetcov in self.steps:
                    wait_for_steps.extend([self.targetcov.job_name(sample.name)])

                self._submit_job(
                    self.abnormal_regions, sample.name,
                    wait_for_steps=wait_for_steps,
                    sample=sample, threads=self.threads_per_sample, genome=sample.genome,
                    caller_names='--caller-names ' + ','.join(caller_names) if caller_names else '',
                    vcfs='--vcfs ' + ','.join(filtered_vcfs) if filtered_vcfs else '')

        if self.varqc_after in self.steps:
            info('VarQC_postVarFilter:')
            for caller in self.bcbio_structure.variant_callers.values():
                info('  ' + caller.name)
                for sample in caller.samples:
                    info('    ' + sample.name)
                    raw_vcf_fpath = sample.find_raw_vcf_by_callername(caller.name)
                    if not raw_vcf_fpath:
                        if sample.phenotype != 'normal':
                            err('Error: raw VCF does not exist: sample ' + sample.name + ', caller "' +
                                caller.name + '". Phenotype = ' + sample.phenotype + '.')
                    else:
                        filt_vcf_fpath = sample.get_filt_vcf_fpath_by_callername(caller.name, gz=True)
                        if not self.varfilter_all and sample.phenotype != 'normal' and not verify_file(filt_vcf_fpath):
                            err('Error: filtered VCF does not exist: sample ' + sample.name + ', caller "' +
                                caller.name + '". Phenotype = ' + sample.phenotype + '.' +
                                ' Note that you need to run VarFilter first, and this step is not in config.')
                        else:
                            self._submit_job(
                                self.varqc_after, sample.name, caller_suf=caller.name, threads=self.threads_per_sample,
                                wait_for_steps=([self.varfilter_all.job_name()] if self.varfilter_all in self.steps else []),
                                vcf=filt_vcf_fpath, sample=sample.name, caller=caller.name, genome=sample.genome)

        if self.varqc_after_summary in self.steps:
            self._submit_job(
                self.varqc_after_summary,
                wait_for_steps=[
                    self.varqc_after.job_name(s.name, v.name)
                    for v in self.bcbio_structure.variant_callers.values()
                    for s in v.samples
                    if self.varqc_after in self.steps])

        if self.seq2c in self.steps:
            self._submit_job(
                self.seq2c,
                wait_for_steps=[],
                genome=self.bcbio_structure.samples[0].genome,
                threads=self.summary_threads)

        if self.targqc_summary in self.steps:
            wait_for_steps = []
            wait_for_steps += [self.targetcov.job_name(s.name) for s in self.bcbio_structure.samples if self.targetcov in self.steps]
            wait_for_steps += [self.ngscat.job_name(s.name) for s in self.bcbio_structure.samples if self.ngscat in self.steps]
            wait_for_steps += [self.qualimap.job_name(s.name) for s in self.bcbio_structure.samples if self.qualimap in self.steps]
            self._submit_job(
                self.targqc_summary,
                wait_for_steps=wait_for_steps,
                threads=self.summary_threads)

        if self.fastqc_summary in self.steps:
            self._submit_job(self.fastqc_summary)

        if self.combined_report in self.steps:
            wait_for_steps = []
            # summaries
            wait_for_steps += [self.varqc_summary.job_name()] if self.varqc_summary in self.steps else []
            wait_for_steps += [self.varqc_after_summary.job_name()] if self.varqc_after_summary in self.steps else []
            wait_for_steps += [self.targqc_summary.job_name()] if self.targqc_summary in self.steps else []
            wait_for_steps += [self.fastqc_summary.job_name()] if self.fastqc_summary in self.steps else []
            # and individual reports too
            wait_for_steps += [self.varqc.job_name(s.name) for s in self.bcbio_structure.samples if self.varqc in self.steps]
            wait_for_steps += [self.targetcov.job_name(s.name) for s in self.bcbio_structure.samples if self.targetcov in self.steps]
            wait_for_steps += [self.ngscat.job_name(s.name) for s in self.bcbio_structure.samples if self.ngscat in self.steps]
            wait_for_steps += [self.qualimap.job_name(s.name) for s in self.bcbio_structure.samples if self.qualimap in self.steps]
            self._submit_job(
                self.combined_report,
                wait_for_steps=wait_for_steps)

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

        info()
        waiting = False
        while True:
            for job in self.jobs_running:
                if not job.is_done and isfile(job.done_marker):
                    job.is_done = True
                    info('Done ' + job.repr)
                    waiting = False

            if not all(j.is_done for j in self.jobs_running):
                info()
                if not waiting:
                    waiting = True
                    info('Waiting for the jobs to be proccesed on the GRID (monitor with qstat)...')
                sleep(30)
            else:
                break

        final_report_job = next((j for j in self.jobs_running if j.step == self.combined_report), None)
        if final_report_job.is_done:
            with open(final_report_job.done_marker) as f:
                html_report_url = f.read()
            info('Final report is saved to ' + html_report_url)

        # Waiting for Seq2C if needed
        not_done = [j for j in self.jobs_running if not j.is_done]
        if not_done:
            info('Waiting for:')
            for job in not_done:
                info('  ' + job.repr)

            while True:
                for job in not_done:
                    if not job.is_done and isfile(job.done_marker):
                        job.is_done = True
                        info('Done ' + job.repr)
                        waiting = False

                if not all(j.is_done for j in self.jobs_running):
                    sleep(30)
                    if not waiting:
                        waiting = True
                        info('Waiting for the jobs to be proccesed on the GRID (monitor with qstat)...')
                else:
                    break


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
                wait_for_steps=[self.varannotate.job_name(sample.name, caller_name)]
                                if self.varannotate in self.steps else [])

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
            cnv_dirpath = join(sample_dirpath, BCBioStructure.cnv_dir)

            for fname in os.listdir(sample_dirpath):
                if any(fname.endswith(s) for s in ['-cn_mops.bed', '-ensemble.bed']):
                    if not isdir(cnv_dirpath): safe_mkdir(cnv_dirpath)
                    try:
                        os.rename(join(sample_dirpath, fname), join(cnv_dirpath, fname))
                    except OSError:
                        pass

            if isdir(cnv_dirpath):
                for fname in os.listdir(cnv_dirpath):
                    if not fname.startswith('.'):
                        src_fpath = join(cnv_dirpath, fname)

                        dst_fname = fname
                        if sample.name not in fname:
                            dst_fname = sample.name + '.' + dst_fname

                        dst_fpath = join(cnv_summary_dirpath, dst_fname)
                        try:
                            if islink(dst_fpath):
                                os.unlink(dst_fpath)
                            symlink_plus(src_fpath, dst_fpath)
                        except OSError:
                            pass


def fix_bed_for_qualimap(bed_fpath, qualimap_bed_fpath):
    with open(qualimap_bed_fpath, 'w') as out, open(bed_fpath) as inn:
        for l in inn:
            fields = l.strip().split('\t')

            if len(fields) < 3:
                continue
            try:
                int(fields[1]), int(fields[2])
            except ValueError:
                continue

            if len(fields) < 4:
                fields.append('.')

            if len(fields) < 5:
                fields.append('0')

            if len(fields) < 6:
                fields.append('+')

            out.write('\t'.join(fields) + '\n')