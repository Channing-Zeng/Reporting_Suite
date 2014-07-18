from genericpath import isfile
import os
import re
import sys
from os.path import join, dirname, abspath, expanduser, basename
from source.calling_process import call
from source.file_utils import verify_dir, verify_file, file_transaction, make_tmpfile
from source.tools_from_cnf import get_tool_cmdline

from source.utils_from_bcbio import file_exists, safe_mkdir, add_suffix
from source.logger import info, err
from source.ngscat.bed_file import verify_bam


basic_dirpath = dirname(dirname(abspath(__file__)))


def run_on_bcbio_final_dir(cnf, bcbio_final_dir, bcbio_cnf):
    return Runner(cnf, bcbio_final_dir, bcbio_cnf).run()


class Step():
    def __init__(self, name, cmdline=None):
        self.name = name
        self.cmdline = cmdline

    def job_name(self, sample=None, caller=None):
        return self.name + ('_' + sample if sample else '') + \
               ('_' + caller if caller else '')


class Steps(list):
    def __init__(self, cnf, steps):
        super(Steps, self).__init__(steps)
        self.cnf = cnf

    @staticmethod
    def __normalize(item):
        return item.lower().replace('_', '').replace('-', '')

    def __contains__(self, item):
        return self.__normalize(item) in \
               [self.__normalize(el) for el in self]

    def step(self, name, script=None, interpreter='python', param_line=None):
        if name not in self:
            return Step(name, None)

        cmd = get_tool_cmdline(self.cnf, interpreter, script or name)
        if not cmd:
            sys.exit(1)

        return Step(name, cmd + ' ' + param_line)


class Runner():
    def __init__(self, cnf, bcbio_final_dir, bcbio_cnf):
        self.dir = bcbio_final_dir
        self.cnf = cnf
        self.bcbio_cnf = bcbio_cnf

        self.threads = str(self.cnf.threads)
        self.steps = Steps(cnf, cnf.steps)
        self.qsub_runner = expanduser(cnf.qsub_runner)

        self.varannotate = None
        self.varqc = None
        self.varqc_after = None
        self.varfilter = None
        self.targetcov = None
        self.ngscat = None
        self.qualimap = None
        self.targetcov_summary = None
        self.varqc_summary = None

        self.date_dirpath = join(bcbio_final_dir, bcbio_cnf.fc_date + '_' + bcbio_cnf.fc_name)
        assert verify_dir(self.date_dirpath)

        self.set_up_steps(cnf)


    def set_up_steps(self, cnf):
        cnfs_line = ' --sys-cnf \'' + self.cnf.sys_cnf + '\' --run-cnf \'' + self.cnf.run_cnf + '\' '

        if cnf.overwrite is not None:
            if cnf.overwrite is True:
                overwrite_line = '-w'
            else:
                overwrite_line = '--reuse'
        else:
            overwrite_line = ''

        spec_params = cnfs_line + ' -t ' + self.threads + ' ' + overwrite_line + ' '

        self.varannotate = self.steps.step(
            name='VarAnnotate',
            script='varannotate.py',
            param_line=spec_params + ' --vcf \'{vcf}\' {bam_cmdline} -o \'{output_dir}\' -s \'{sample}\' --work-dir \'' + join(cnf.work_dir, 'varannotate') + '\'')
        self.varqc = self.steps.step(
            name='VarQC',
            script='varqc.py',
            param_line=spec_params + ' --vcf \'{vcf}\' -o \'{output_dir}\' -s \'{sample}\' --work-dir \'' + join(cnf.work_dir, 'varqc') + '\'')
        self.varqc_after = self.steps.step(
            name='VarQC_after',
            script='varqc.py',
            param_line=spec_params + ' --vcf \'{vcf}\' -o \'{output_dir}\' -s \'{sample}\' --work-dir \'' + join(cnf.work_dir, 'varqc_after') + '\'')
        self.varfilter = self.steps.step(
            name='VarFilter',
            script='varfilter.py',
            param_line=spec_params + ' --vcf \'{vcf}\' -o \'{output_dir}\' -s \'{sample}\' --work-dir \'' + join(cnf.work_dir, 'varfilter') + '\'')
        self.targetcov = self.steps.step(
            name='TargetCov',
            script='targetcov.py',
            param_line=spec_params + ' --bam \'{bam}\' --bed \'{bed}\' -o \'{output_dir}\' -s \'{sample}\' --work-dir \'' + join(cnf.work_dir, 'targetcov') + '\'')
        self.ngscat = self.steps.step(
            name='NGScat',
            script='ngscat.py',
            param_line=spec_params + ' --bam \'{bam}\' --bed \'{bed}\' -o \'{output_dir}\' -s \'{sample}\' --saturation y --work-dir \'' + join(cnf.work_dir, 'ngscat') + '\'')
        self.qualimap = self.steps.step(
            name='QualiMap',
            interpreter=None,
            param_line=' bamqc -nt ' + self.threads + ' --java-mem-size=24G -nr 5000 -bam \'{bam}\' -outdir \'{output_dir}\' -gff \'{bed}\' -c -gd HUMAN')

        all_suffixes = set()
        for s_info in self.bcbio_cnf.details:
            all_suffixes |= set(s_info['algorithm'].get('variantcaller')) or set()

        self.varqc_summary = self.steps.step(
            name='VarQC_summary',
            script='varqc_summary.py',
            param_line=cnfs_line + ' -o \'{output_dir}\' -d \'' + self.dir + '\' -s \'{samples}\''
                       ' -n varqc --vcf-suf ' + ','.join(all_suffixes) + ' --work-dir \'' + join(cnf.work_dir, 'varqc_summary') + '\'')
        self.targetcov_summary = self.steps.step(
            name='TargetCov_summary',
            script='targetcov_summary.py',
            param_line=cnfs_line + ' -o \'{output_dir}\' -d \'' + self.dir + '\' -s \'{samples}\''
                       ' -n targetcov --work-dir \'' + join(cnf.work_dir, 'targetcov_summary') + '\'')

    def submit_if_needed(self, step, sample_name='', suf=None, create_dir=True,
                         out_fpath=None, wait_for_steps=list(), **kwargs):

        output_dirpath = self.dir
        if sample_name:
            output_dirpath = join(output_dirpath, sample_name)
        else:
            output_dirpath = self.date_dirpath
            # if output_dirpath != self.dir:
            #     create_dir = False

        log_fpath = join(output_dirpath, step.job_name(sample_name, suf).lower() + '.log')

        if create_dir:
            output_dirpath = join(output_dirpath, step.name.lower())
            log_fpath = join(output_dirpath, 'log' + ('_' + suf if suf else ''))

        # Skipping steps
        if not step or step.name not in self.steps:
            return output_dirpath

        if create_dir:
            safe_mkdir(output_dirpath)

        out_fpath = out_fpath or log_fpath

        if isfile(out_fpath):
            try:
                os.remove(out_fpath)
            except OSError:
                err('Cannot remove log file ' + out_fpath + ', probably permission denied.')

        if isfile(log_fpath):
            try:
                os.remove(log_fpath)
            except OSError:
                err('Cannot remove log file ' + out_fpath + ', probably permission denied.')

        hold_jid_line = '-hold_jid ' + ','.join(wait_for_steps or ['_'])

        job_name = step.job_name(sample_name, suf)

        params = dict({'output_dir': output_dirpath}.items() +
                      self.__dict__.items() + kwargs.items())

        runner_script = self.qsub_runner

        cmdline = step.cmdline.format(**params)

        qsub = get_tool_cmdline(self.cnf, 'qsub')

        threads = str(self.threads)

        queue = self.cnf.queue

        qsub_cmdline = (
            '{qsub} -pe smp {threads} -S /bin/bash -q {queue} '
            '-j n -o {out_fpath} -e {log_fpath} {hold_jid_line} '
            '-N {job_name} {runner_script} "{cmdline}"'.format(**locals()))

        if self.cnf.verbose:
            info(step.name)
            info(qsub_cmdline)
        else:
            print step.name,

        call(self.cnf, qsub_cmdline, silent=True)

        if self.cnf.verbose:
            info()

        return output_dirpath

    def qualimap_bed(self, bed_fpath):
        if 'QualiMap' in self.steps and bed_fpath:
            qualimap_bed_fpath = join(self.cnf.work_dir, 'tmp_qualimap.bed')

            with open(qualimap_bed_fpath, 'w') as out, open(bed_fpath) as inn:
                for l in inn:
                    fields = l.strip().split('\t')

                    if len(fields) < 3:
                        continue
                    try:
                        n = int(fields[1])
                        n = int(fields[2])
                    except ValueError:
                        continue

                    if len(fields) < 4:
                        fields.append('-')

                    if len(fields) < 5:
                        fields.append('0')

                    if len(fields) < 6:
                        fields.append('+')

                    out.write('\t'.join(fields) + '\n')

            bed_fpath = qualimap_bed_fpath
        return bed_fpath

    def run(self):
        for sample_info in self.bcbio_cnf.details:
            sample = sample_info['description']

            bed_fpath = sample_info['algorithm'].get('variant_regions')
            qualimap_bed_fpath = self.qualimap_bed(bed_fpath)

            info(sample)
            if not self.cnf.verbose:
                info(ending='')

            sample_dirpath = join(self.dir, sample)
            if not verify_dir(sample_dirpath):
                sys.exit(1)

            bam_fpath = join(sample_dirpath, sample + '-ready.bam')
            if 'TargetCov' in self.steps or 'NGScat' in self.steps or 'Qualimap' in self.steps:
                if not verify_bam(bam_fpath):
                    sys.exit(1)
            else:
                if not file_exists(bam_fpath):
                    bam_fpath = None

            for variant_caller in sample_info['algorithm'].get('variantcaller') or []:
                vcf_fpath = join(sample_dirpath, sample + '-' + variant_caller + '.vcf')

                if not file_exists(vcf_fpath) and file_exists(vcf_fpath + '.gz'):
                    gz_vcf_fpath = vcf_fpath + '.gz'
                    gunzip = get_tool_cmdline(self.cnf, 'gunzip')
                    cmdline = '{gunzip} -c {gz_vcf_fpath}'.format(**locals())
                    call(self.cnf, cmdline, output_fpath=vcf_fpath)
                    info()

                if not file_exists(vcf_fpath):
                    info('No ' + vcf_fpath + ', skipping')
                    continue

                if not verify_file(vcf_fpath):
                    sys.exit(1)

                self._process_vcf(sample, sample_dirpath, bam_fpath, vcf_fpath, variant_caller)

            self.submit_if_needed(self.targetcov, sample, bam=bam_fpath, bed=bed_fpath, sample=sample)

            self.submit_if_needed(self.ngscat, sample, bam=bam_fpath, bed=bed_fpath, sample=sample)

            self.submit_if_needed(self.qualimap, sample, bam=bam_fpath, bed=qualimap_bed_fpath, sample=sample)

            if self.cnf.verbose:
                info('-' * 70)
            else:
                print ''
                info()

        if not self.cnf.verbose:
            info('', ending='')

        all_variantcallers = set()
        for s_info in self.bcbio_cnf.details:
            all_variantcallers |= set(s_info['algorithm'].get('variantcaller')) or set()

        samples_fpath = abspath(join(self.cnf.work_dir, 'samples.txt'))
        #if not isfile(samples_fpath):
        with open(samples_fpath, 'w') as f:
            for sample_info in self.bcbio_cnf.details:
                sample = sample_info['description']
                f.write(sample + '\n')

        self.submit_if_needed(
            self.varqc_summary,
            wait_for_steps=[
                self.varqc.job_name(d['description'], v)
                for d in self.bcbio_cnf.details
                for v in all_variantcallers],
            samples=samples_fpath)

        self.submit_if_needed(
            self.targetcov_summary,
            wait_for_steps=[
                self.targetcov.job_name(d['description'])
                for d in self.bcbio_cnf.details],
            samples=samples_fpath)

        if not self.cnf.verbose:
            print ''
        if self.cnf.verbose:
            info('Done.')

    def _process_vcf(self, sample, sample_dirpath, bam_fpath, vcf_fpath, suf):
        self.submit_if_needed(self.varqc, sample, suf=suf, vcf=vcf_fpath, sample=sample + '-' + suf)

        bam_cmdline = '--bam ' + bam_fpath if bam_fpath else ''
        anno_dirpath = self.submit_if_needed(self.varannotate, sample, suf=suf, vcf=vcf_fpath, bam_cmdline=bam_cmdline, sample=sample + '-' + suf)

        annotated_vcf_fpath = join(anno_dirpath, basename(add_suffix(vcf_fpath, 'anno')))

        filter_dirpath = self.submit_if_needed(
            self.varfilter, sample, suf=suf,
            wait_for_steps=[self.varannotate.job_name(sample, suf)],
            vcf=annotated_vcf_fpath, sample=sample + '-' + suf)

        filtered_vcf_fpath = join(filter_dirpath, basename(add_suffix(annotated_vcf_fpath, 'filt')))

        self.submit_if_needed(
            self.varqc_after, sample, suf=suf,
            wait_for_steps=[self.varfilter.job_name(sample, suf)] if self.varfilter.name in self.steps else [],
            vcf=filtered_vcf_fpath, sample=sample + '-' + suf)