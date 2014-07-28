import base64
from collections import defaultdict
from genericpath import isfile, isdir
import hashlib
import os
import re
import sys
from os.path import join, dirname, abspath, expanduser, basename, pardir
from source.calling_process import call
from source.file_utils import verify_dir, verify_file, file_transaction, make_tmpfile
from source.tools_from_cnf import get_tool_cmdline, get_script_cmdline, tool_cmdline

from source.utils_from_bcbio import file_exists, safe_mkdir, add_suffix
from source.logger import info, err, critical
from source.ngscat.bed_file import verify_bam


basic_dirpath = dirname(dirname(abspath(__file__)))


def run_on_bcbio_final_dir(cnf, bcbio_final_dir, bcbio_cnf):
    return Runner(cnf, bcbio_final_dir, bcbio_cnf).run()


def _normalize(name):
    return name.lower().replace('_', '').replace('-', '')


class Step():
    def __init__(self, cnf, run_id, name, script, interpreter=None, short_name=None, paramln=None):
        self.name = name
        self.cnf = cnf
        assert run_id
        self.run_id_ = run_id
        self.short_name = short_name or name
        self.param_line = paramln
        self.run_id = None
        self.script = script
        self.interpreter = interpreter

    def job_name(self, sample=None, caller=None):
        return self.short_name.upper() + '_' + self.run_id_ + \
               ('_' + sample if sample else '') + \
               ('_' + caller if caller else '')

class Steps(list):
    def __init__(self):
        super(Steps, self).__init__()

    def add_step(self, step):
        self.append(step)

    def extend(self, iterable):
        for step in iterable:
            self.add_step(step)


class Runner():
    def __init__(self, cnf, bcbio_final_dir, bcbio_cnf):
        self.dir = bcbio_final_dir
        self.cnf = cnf
        self.bcbio_cnf = bcbio_cnf

        hasher = hashlib.sha1(bcbio_final_dir)
        self.run_id = base64.urlsafe_b64encode(hasher.digest()[0:8])[:-1]

        self.threads = str(self.cnf.threads)
        self.qsub_runner = expanduser(cnf.qsub_runner)

        self.date_dirpath = join(bcbio_final_dir, bcbio_cnf.fc_date + '_' + bcbio_cnf.fc_name)
        if not verify_dir(self.date_dirpath):
            critical('The project directory must have format {fc_date}_{fc_name}, here: ' + self.date_dirpath)

        self.steps = Steps()
        self.vardict_steps = Steps()

        self.set_up_steps(cnf, self.run_id)

        def contains(x, xs):
            return _normalize(x) in [_normalize(y) for y in (xs or [])]

        self.vardict_steps.extend([s for s in [
            self.vardict,
            self.testsomatic,
            self.var_to_vcf_somatic,
            self.varqc,
            self.varannotate,
            self.varfilter,
            self.varqc_after,
            self.varqc_summary] if contains(s.name, cnf.vardict_steps)
        ])

        self.steps.extend([s for s in [
            self.varqc,
            self.varannotate,
            self.varfilter,
            self.varqc_after,
            self.varqc_summary,
            self.targetcov,
            self.ngscat,
            self.qualimap,
            self.targetcov_summary] if contains(s.name, cnf.steps)
        ])

    def set_up_steps(self, cnf, run_id):
        cnfs_line = ' --sys-cnf \'' + self.cnf.sys_cnf + '\' --run-cnf \'' + self.cnf.run_cnf + '\' '

        if cnf.overwrite is not None:
            if cnf.overwrite is True:
                overwrite_line = '-w'
            else:
                overwrite_line = '--reuse'
        else:
            overwrite_line = ''

        spec_params = cnfs_line + ' -t ' + self.threads + ' ' + overwrite_line + ' '

        # VARDICT
        af_thr = str(cnf.variant_filtering.min_freq)
        self.vardict = Step(cnf, run_id,
            name='VarDict',
            interpreter='perl',
            script='vardict_pl',
            paramln=' -G ' + cnf.genome.seq + ' -f ' + af_thr + ' -N {tumor_name} -b \'{tumor_bam}|{normal_bam}\''
                    ' -z -F -c 1 -S 2 -E 3 -g 4 {bed} > {vars_txt}'
        )
        self.testsomatic = Step(cnf, run_id,
            name='TestSomatic',
            script='testsomatic_r',
            paramln=' < {vars_txt} > {somatic_vars_txt}',
        )
        self.var_to_vcf_somatic = Step(cnf, run_id,
            name='Var2Vcf_Somatic',
            interpreter='perl',
            script='var2vcf_somatic_pl',
            paramln=' -N \'{tumor_name}|{normal_name}\' -f ' +
                    af_thr + ' < {somatic_vars_txt} > {vardict_vcf}',
        )
        # END VARDICT

        self.varannotate = Step(cnf, run_id,
            name='VarAnnotate', short_name='va',
            interpreter='python',
            script='varannotate',
            paramln=spec_params + ' --vcf \'{vcf}\' {bam_cmdline} '
                    '-o \'{output_dir}\' -s \'{sample}\' --work-dir \'' + join(cnf.work_dir, 'varannotate') + '\''
        )
        self.varqc = Step(cnf, run_id,
            name='VarQC', short_name='vq',
            interpreter='python',
            script='varqc',
            paramln=spec_params + ' --vcf \'{vcf}\' -o \'{output_dir}\''
                    ' -s \'{sample}\' --work-dir \'' + join(cnf.work_dir, 'varqc') + '\''
        )
        self.varqc_after = Step(cnf, run_id,
            name='VarQC_after', short_name='vqa',
            interpreter='python',
            script='varqc',
            paramln=spec_params + ' --vcf \'{vcf}\' -o '
                    '\'{output_dir}\' -s \'{sample}\' --work-dir \'' + join(cnf.work_dir, 'varqc_after') + '\''
        )
        self.varfilter = Step(cnf, run_id,
            name='VarFilter', short_name='vf',
            interpreter='python',
            script='varfilter',
            paramln=spec_params + ' --vcf \'{vcf}\' -o \'{output_dir}\' '
                    '-s \'{sample}\' --work-dir \'' + join(cnf.work_dir, 'varfilter') + '\''
        )
        self.targetcov = Step(cnf, run_id,
            interpreter='python',
            script='targetcov',
            name='TargetCov', short_name='tc',
            paramln=' --bam \'{bam}\' --bed \'{bed}\' -o \'{output_dir}\' '
                    '-s \'{sample}\' --work-dir \'' + join(cnf.work_dir, 'targetcov') + '\''
        )
        self.ngscat = Step(cnf, run_id,
            interpreter='python',
            script='ngscat',
            name='NGScat', short_name='nc',
            paramln=spec_params + ' --bam \'{bam}\' --bed \'{bed}\' '
                    '-o \'{output_dir}\' -s \'{sample}\' --saturation y --work-dir \'' + join(cnf.work_dir, 'ngscat') + '\''
        )
        self.qualimap = Step(cnf, run_id,
            script='qualimap',
            name='QualiMap', short_name='qm',
            paramln=' bamqc -nt ' + self.threads + ' --java-mem-size=24G -nr 5000 '
                    '-bam \'{bam}\' -outdir \'{output_dir}\' -gff \'{bed}\' -c -gd HUMAN'
        )

        all_suffixes = set()
        for s_info in self.bcbio_cnf.details:
            all_suffixes |= set(s_info['algorithm'].get('variantcaller')) or set()

        self.varqc_summary = Step(cnf, run_id,
            name='VarQC_summary', short_name='vqs',
            interpreter='python',
            script='varqc_summary',
            paramln=cnfs_line + ' -o \'{output_dir}\' -d \''
                    + self.dir + '\' -s \'{samples}\' -n varqc --vcf-suf ' + ','.join(all_suffixes) +
                    ' --work-dir \'' + join(cnf.work_dir, 'varqc_summary') + '\''
        )
        self.targetcov_summary = Step(cnf, run_id,
            name='TargetCov_summary', short_name='tcs',
            interpreter='python',
            script='targetcov_summary',
            paramln=cnfs_line + ' -o \'{output_dir}\' -d \''
                    + self.dir + '\' -s \'{samples}\' -n targetcov --work-dir \'' +
                    join(cnf.work_dir, 'targetcov_summary') + '\''
        )

    def step_output_dir_and_log_paths(self, step, sample_name, suf=None, dir_name=None):
        output_dirpath = self.dir
        if sample_name:
            output_dirpath = join(output_dirpath, sample_name)
        else:
            output_dirpath = self.date_dirpath

        output_dirpath = abspath(output_dirpath)

        log_fpath = join(output_dirpath, (step.name +
                                          ('_' + sample_name if sample_name else '') +
                                          ('_' + suf if suf else '')).lower() + '.log')

        if dir_name is None:
            dir_name = step.name
        dir_name = dir_name.lower()
        if dir_name != '':
            output_dirpath = join(output_dirpath, dir_name)
            log_fpath = join(output_dirpath, 'log' + ('_' + suf if suf else ''))

        return output_dirpath, log_fpath


    def submit(self, step, sample_name='', suf=None, dir_name=None,
               out_fpath=None, wait_for_steps=list(), threads=None, **kwargs):

        output_dirpath, log_fpath = self.step_output_dir_and_log_paths(step, sample_name, suf, dir_name)

        if not isdir(output_dirpath):
            safe_mkdir(output_dirpath)

        out_fpath = out_fpath or log_fpath

        if isfile(out_fpath):
            try:
                os.remove(out_fpath)
            except OSError:
                err('Cannot remove log file ' + out_fpath + ', probably permission denied.')

        if log_fpath and isfile(log_fpath):
            try:
                os.remove(log_fpath)
            except OSError:
                err('Cannot remove log file ' + out_fpath + ', probably permission denied.')

        hold_jid_line = '-hold_jid ' + ','.join(wait_for_steps or ['_'])

        job_name = step.job_name(sample_name, suf)

        params = dict({'output_dir': output_dirpath}.items() +
                      self.__dict__.items() + kwargs.items())

        runner_script = self.qsub_runner

        tool_cmdline = get_tool_cmdline(self.cnf, step.interpreter, step.script)
        if not tool_cmdline:
            sys.exit(1)
        cmdline = tool_cmdline + ' ' + step.param_line.format(**params)

        qsub = get_tool_cmdline(self.cnf, 'qsub')

        threads = str(threads or self.threads)

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

            return qualimap_bed_fpath
        else:
            return bed_fpath

    def run(self):
        batches = defaultdict(dict)

        for sample_info in self.bcbio_cnf.details:
            sample = sample_info['description']

            info('Processing "' + sample + '"')
            if not self.cnf.verbose:
                info(ending='')

            sample_dirpath = join(self.dir, sample)
            if not verify_dir(sample_dirpath):
                sys.exit(1)

            bed_fpath = sample_info['algorithm'].get('variant_regions')
            bam_fpath = join(sample_dirpath, sample + '-ready.bam')

            if any(step in self.steps for step in [self.targetcov, self.qualimap, self.ngscat]) \
                    or self.vardict in self.vardict_steps:

                if not verify_bam(bam_fpath) or not verify_file(bed_fpath):
                    sys.exit(1)
                else:
                    bam_fpath = abspath(expanduser(bam_fpath))
                    bed_fpath = abspath(expanduser(bed_fpath))

                if 'QualiMap' in self.steps:
                    bed_fpath = self.qualimap_bed(bed_fpath)
            else:
                if not file_exists(bam_fpath):
                    bam_fpath = None
                if not file_exists(bed_fpath):
                    bed_fpath = None

            phenotype = None
            if 'metadata' in sample_info:
                phenotype = sample_info['metadata']['phenotype']

                if self.vardict_steps:
                    batch_names = sample_info['metadata']['batch']
                    if isinstance(batch_names, basestring):
                        batch_names = [batch_names]

                    for batch_name in batch_names:
                        batches[batch_name]['bed'] = bed_fpath
                        if phenotype == 'normal':
                            if batches[batch_name].get('normal'):
                                critical('Multiple normal samples for batch ' + batch_name)
                            batches[batch_name]['normal'] = sample, bam_fpath
                        elif phenotype == 'tumor':
                            if 'tumor' not in batches[batch_name]:
                                batches[batch_name]['tumor'] = dict()
                            batches[batch_name]['tumor'][sample] = bam_fpath

            for variant_caller in sample_info['algorithm'].get('variantcaller') or []:
                vcf_fpath = join(sample_dirpath, sample + '-' + variant_caller + '.vcf')

                if not file_exists(vcf_fpath) and file_exists(vcf_fpath + '.gz'):
                    gz_vcf_fpath = vcf_fpath + '.gz'
                    gunzip = get_tool_cmdline(self.cnf, 'gunzip')
                    cmdline = '{gunzip} -c {gz_vcf_fpath}'.format(**locals())
                    call(self.cnf, cmdline, output_fpath=vcf_fpath)
                    info()

                if not file_exists(vcf_fpath):
                    if phenotype != 'normal':
                        err('No ' + vcf_fpath + ', skipping')
                else:
                    if not verify_file(vcf_fpath):
                        sys.exit(1)

                    self._process_vcf(sample, bam_fpath, vcf_fpath, variant_caller)

            for step in [self.targetcov, self.qualimap, self.ngscat]:
                if step in self.steps:
                    self.submit(step, sample, bam=bam_fpath, bed=bed_fpath, sample=sample)

            if self.cnf.verbose:
                info('-' * 70)
            else:
                print ''
                info()

        if not self.cnf.verbose:
            info('', ending='')

        for batch_name, batch in batches.items():
            normal_name, normal_bam_fpath = batch['normal']
            bed_fpath = batch['bed']
            for tumor_name, tumor_bam_fpath in batch['tumor'].items():
                output_dirpath, _ = self.step_output_dir_and_log_paths(self.vardict, tumor_name, dir_name=self.vardict.name)
                vars_txt = join(output_dirpath, 'vardict.txt')

                if self.vardict in self.vardict_steps:
                    self.submit(
                        self.vardict, tumor_name, suf='vardict',
                        tumor_name=tumor_name,
                        normal_name=normal_name,
                        tumor_bam=tumor_bam_fpath,
                        normal_bam=normal_bam_fpath,
                        bed=bed_fpath,
                        vars_txt=vars_txt)

                somatic_vars_txt = join(output_dirpath, 'somatic_variants.txt')
                if self.testsomatic in self.vardict_steps:
                    self.submit(
                        self.testsomatic, tumor_name, suf='testsomatic', dir_name=self.vardict.name,
                        vars_txt=vars_txt,
                        somatic_vars_txt=somatic_vars_txt)

                vardict_vcf = join(output_dirpath, 'somatic_variants-vardict_standalone.vcf')
                if self.var_to_vcf_somatic in self.vardict_steps:
                    self.submit(
                        self.var_to_vcf_somatic, tumor_name, suf='var2vcf', dir_name=self.vardict.name,
                        tumor_name=tumor_name,
                        normal_name=normal_name,
                        somatic_vars_txt=somatic_vars_txt,
                        vardict_vcf=vardict_vcf)

                self._process_vcf(
                    tumor_name, tumor_bam_fpath, vardict_vcf, 'vardict_standalone',
                    dir_name=self.vardict.name, steps=self.vardict_steps)

        all_variantcallers = set()
        for s_info in self.bcbio_cnf.details:
            all_variantcallers |= set(s_info['algorithm'].get('variantcaller')) or set()

        samples_fpath = abspath(join(self.cnf.work_dir, 'samples.txt'))
        #if not isfile(samples_fpath):
        samples_num = 0
        with open(samples_fpath, 'w') as f:
            for sample_info in self.bcbio_cnf.details:
                sample = sample_info['description']
                f.write(sample + '\n')
                samples_num += 1

        if self.varqc_summary in self.steps:
            self.submit(
                self.varqc_summary,
                wait_for_steps=[
                    self.varqc.job_name(d['description'], v)
                    for d in self.bcbio_cnf.details
                    for v in all_variantcallers],
                # threads=samples_num + 1,
                samples=samples_fpath)

        if self.targetcov_summary in self.steps:
            self.submit(
                self.targetcov_summary,
                wait_for_steps=[
                    self.targetcov.job_name(d['description'])
                    for d in self.bcbio_cnf.details],
                # threads=samples_num + 1,
                samples=samples_fpath)

        if not self.cnf.verbose:
            print ''
        if self.cnf.verbose:
            info('Done.')

    def _process_vcf(self, sample, bam_fpath, vcf_fpath, caller, dir_name=None, steps=None):
        steps = steps or self.steps

        if self.varqc in steps:
            self.submit(self.varqc, sample, suf=caller, dir_name=dir_name,
                        vcf=vcf_fpath, sample=sample + '-' + caller)

        bam_cmdline = '--bam ' + bam_fpath if bam_fpath else ''

        if self.varannotate in steps:
            self.submit(
                self.varannotate, sample, suf=caller, dir_name=dir_name, vcf=vcf_fpath,
                bam_cmdline=bam_cmdline, sample=sample + '-' + caller)

        anno_dirpath, _ = self.step_output_dir_and_log_paths(self.varannotate, sample, suf=caller, dir_name=dir_name)
        annotated_vcf_fpath = join(anno_dirpath, basename(add_suffix(vcf_fpath, 'anno')))

        if self.varfilter in steps:
            self.submit(
                self.varfilter, sample, suf=caller, dir_name=dir_name,
                wait_for_steps=[self.varannotate.job_name(sample, caller)],
                vcf=annotated_vcf_fpath, sample=sample + '-' + caller)

        filter_dirpath, _ = self.step_output_dir_and_log_paths(self.varfilter, sample, suf=caller, dir_name=dir_name)
        filtered_vcf_fpath = join(filter_dirpath, basename(add_suffix(annotated_vcf_fpath, 'filt')))

        if self.varqc_after in steps:
            self.submit(
                self.varqc_after, sample, suf=caller, dir_name=dir_name,
                wait_for_steps=[self.varfilter.job_name(sample, caller)] if self.varfilter.name in self.steps else [],
                vcf=filtered_vcf_fpath, sample=sample + '-' + caller)