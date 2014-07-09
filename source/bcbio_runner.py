from genericpath import isfile
import os
import re
import sys
from os.path import join, dirname, abspath, expanduser, basename
from source.calling_process import call
from source.file_utils import verify_dir, verify_file, file_transaction, make_tmpfile
from source.tools_from_cnf import get_tool_cmdline

from source.utils_from_bcbio import file_exists, safe_mkdir, add_suffix
from source.logger import info
from source.ngscat.bed_file import verify_bam


basic_dirpath = dirname(dirname(abspath(__file__)))


def run_on_bcbio_final_dir(cnf, bcbio_final_dir, samples_fpath,
                           bed_fpath, vcf_sufs):
    return Runner(cnf, bcbio_final_dir, samples_fpath,
                  bed_fpath, vcf_sufs).run()


class Step():
    def __init__(self, name, cmdline=None):
        self.name = name
        self.cmdline = cmdline

    def job_name(self, sample=None):
        return self.name + ('_' + sample if sample else '')


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
            return None

        cmd = get_tool_cmdline(self.cnf, interpreter, script or name)
        if not cmd:
            sys.exit(1)

        return Step(name, cmd + ' ' + param_line)


class Runner():
    def __init__(self, cnf, bcbio_final_dir, samples_fpath, bed_fpath, vcf_sufs):
        self.dir = bcbio_final_dir
        self.cnf = cnf
        self.bed = bed_fpath
        self.sufs = vcf_sufs
        self.threads = str(self.cnf.threads)
        self.steps = Steps(cnf, cnf.steps)
        self.qsub_runner = expanduser(cnf.qsub_runner)

        self.varannotate = None
        self.varqc = None
        self.varfilter = None
        self.targetcov = None
        self.ngscat = None
        self.qualimap = None
        self.targetcov_summary = None
        self.varqc_summary = None

        self.samples_fpath = samples_fpath
        with open(samples_fpath) as sample_f:
            self.samples = [s.strip() for s in sample_f.readlines()
                            if s and s.strip() and not s.startswith('#')]

        date_dir_pattern = re.compile(r'\d\d\d\d\-\d\d-\d\d.*')
        self.date_dirpath = next((join(cnf.bcbio_final_dir, dir_name)
                             for dir_name in os.listdir(cnf.bcbio_final_dir)
                             if date_dir_pattern.match(dir_name)),
                            cnf.bcbio_final_dir)

        self.set_up_steps(cnf)


    def set_up_steps(self, cnf):
        cnfs_line = '--sys-cnf \'' + self.cnf.sys_cnf + '\' --run-cnf \'' + self.cnf.run_cnf + '\''
        overwrite_line = '-w' if self.cnf.overwrite else ''
        spec_params = cnfs_line + ' -t ' + self.threads + ' ' + overwrite_line + ' '

        self.varannotate = self.steps.step(
            name='VarAnnotate',
            script='varannotate.py',
            param_line=spec_params + ' --vcf \'{vcf}\' --bam \'{bam}\' -o \'{output_dir}\' -s \'{sample}\' --work-dir \'' + join(cnf.work_dir, 'varannotate') + '\'')
        self.varqc = self.steps.step(
            name='VarQC',
            script='varqc.py',
            param_line=spec_params + ' --vcf \'{vcf}\' -o \'{output_dir}\' -s \'{sample}\' --work-dir \'' + join(cnf.work_dir, 'varqc') + '\'')
        self.varfilter = self.steps.step(
            name='VarFilter',
            script='varfilter.py',
            param_line=spec_params + ' --vcf \'{vcf}\' -o \'{output_dir}\' -s \'{sample}\' --clean --work-dir \'' + join(cnf.work_dir, 'varfilter') + '\'')
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
            param_line=' bamqc -nt ' + self.threads + ' --java-mem-size=24G -nr 5000 -bam '
                       '\'{bam}\' -outdir \'{output_dir}\' -gff \'{bed}\' -c -gd HUMAN')

        if self.varqc:
            self.steps.append('VarQC_summary')
            self.varqc_summary = self.steps.step(
                name='VarQC_summary',
                script='varqc_summary.py',
                param_line=' -o \'{output_dir}\' -d \'' + self.dir + '\' -s \'{samples}\' -n ' +
                           self.varqc.name + ' --vcf-suf ' + ','.join(self.sufs) +
                           ' --work-dir \'' + join(cnf.work_dir, 'varqc_summary') + '\'')
        if self.targetcov:
            self.steps.append('TargetCov_summary')
            self.targetcov_summary = self.steps.step(
                name='TargetCov_summary',
                script='targetcov_summary.py',
                param_line=' -o \'{output_dir}\' -d \'' + self.dir +
                           '\' -s \'{samples}\' -n ' + self.targetcov.name +
                           ' --work-dir \'' + join(cnf.work_dir, 'targetcov_summary') + '\'')

    def submit(self, step, sample_name='', create_dir=False, out_fpath=None,
               wait_for_steps=list(), **kwargs):
        if not step or step.name not in self.steps:
            return None

        output_dirpath = self.dir
        if sample_name:
            output_dirpath = join(output_dirpath, sample_name)
        else:
            output_dirpath = self.date_dirpath
            # if output_dirpath != self.dir:
            #     create_dir = False

        log_fpath = join(output_dirpath, step.name.lower() + '.log')

        if create_dir:
            output_dirpath = join(output_dirpath, step.name.lower())
            safe_mkdir(output_dirpath)
            log_fpath = join(output_dirpath, 'log')

        out_fpath = out_fpath or log_fpath

        if isfile(out_fpath):
            os.remove(out_fpath)
        if isfile(log_fpath):
            os.remove(log_fpath)

        hold_jid_line = '-hold_jid ' + ','.join(wait_for_steps or ['_'])

        job_name = step.job_name(sample_name) if sample_name else step.name

        params = dict({'output_dir': output_dirpath}.items() +
                      self.__dict__.items() + kwargs.items())

        runner_script = self.qsub_runner

        cmdline = step.cmdline.format(**params)

        qsub = get_tool_cmdline(self.cnf, 'qsub')

        threads = str(self.threads)

        qsub_cmdline = (
            '{qsub} -pe smp {threads} -S /bin/bash -q batch.q '
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

    def run(self):
        qualimap_bed_fpath = join(self.cnf.work_dir, 'tmp_qualimap.bed')
        with open(qualimap_bed_fpath, 'w') as out, open(self.bed) as inn:
            for l in inn:
                ts = l.strip().split('\t')
                if len(ts) < 5:
                    ts += ['0']
                if len(ts) < 6:
                    ts += ['+']
                out.write('\t'.join(ts) + '\n')

        for sample in self.samples:
            info(sample)
            if not self.cnf.verbose:
                info(ending='')

            sample_dirpath = join(self.dir, sample)
            if not verify_dir(sample_dirpath):
                sys.exit(1)

            bam_fpath = join(sample_dirpath, sample + '-ready.bam')
            if not verify_bam(bam_fpath):
                sys.exit(1)

            for vcf_suf in self.sufs:
                vcf_fpath = join(sample_dirpath, sample + '-' + vcf_suf + '.vcf')

                if not file_exists(vcf_fpath) and file_exists(vcf_fpath + '.gz'):
                    gz_vcf_fpath = vcf_fpath + '.gz'
                    gunzip = get_tool_cmdline(self.cnf, 'gunzip')
                    cmdline = '{gunzip} -c {gz_vcf_fpath}'.format(**locals())
                    call(self.cnf, cmdline, output_fpath=vcf_fpath)
                    info()

                if 'mutect' in vcf_suf and not file_exists(vcf_fpath):
                    continue

                if not verify_file(vcf_fpath):
                    sys.exit(1)

                self._process_vcf(sample, sample_dirpath, bam_fpath, vcf_fpath, vcf_suf)

            self.submit(self.targetcov, sample, True, bam=bam_fpath, bed=self.bed, sample=sample)

            self.submit(self.ngscat, sample, True, bam=bam_fpath, bed=self.bed, sample=sample)

            self.submit(self.qualimap, sample, True, bam=bam_fpath, bed=qualimap_bed_fpath, sample=sample)

            if self.cnf.verbose:
                info('-' * 70)
            else:
                print ''
                info()

        if not self.cnf.verbose:
            info('', ending='')

        self.submit(
            self.varqc_summary,
            create_dir=True,
            wait_for_steps=[self.varqc.job_name(s) for s in self.samples],
            samples=self.samples_fpath)

        self.submit(
            self.targetcov_summary,
            create_dir=True,
            wait_for_steps=[self.targetcov.job_name(s) for s in self.samples],
            samples=self.samples_fpath)

        if not self.cnf.verbose:
            print ''
        if self.cnf.verbose:
            info('Done.')

    def _process_vcf(self, sample, sample_dirpath, bam_fpath, vcf_fpath, suf):
        if self.varqc:
            self.submit(
                self.varqc, sample, True,
                vcf=vcf_fpath, sample=sample + '-' + suf)

        if self.varannotate:
            anno_dirpath = self.submit(
                self.varannotate, sample, True,
                vcf=vcf_fpath, bam=bam_fpath, sample=sample + '-' + suf)
            annotated_vcf_fpath = join(anno_dirpath, basename(add_suffix(vcf_fpath, 'anno')))

            if self.varfilter:
                filtered_vcf_fpath = join(sample_dirpath,
                    basename(add_suffix(annotated_vcf_fpath, 'filt')))
                if file_exists(filtered_vcf_fpath):
                    os.remove(filtered_vcf_fpath)
                self.submit(self.varfilter, sample, False,
                    wait_for_steps=[self.varannotate.job_name(sample)],
                    vcf=annotated_vcf_fpath,
                    sample=sample + '-' + suf)