import hashlib
import os
import sys
import base64
from collections import defaultdict
from os.path import join, dirname, abspath, expanduser, basename, pardir, isfile, isdir, exists
from source.calling_process import call
from source.file_utils import verify_dir, verify_file, add_suffix
from source.tools_from_cnf import get_tool_cmdline

from source.file_utils import file_exists, safe_mkdir
from source.logger import info, err, critical
from source.ngscat.bed_file import verify_bam


# basic_dirpath = join(dirname(abspath(__file__)), pardir)


def run_on_bcbio_final_dir(cnf, bcbio_final_dir, bcbio_cnf):
    return Runner(cnf, bcbio_final_dir, bcbio_cnf).run()


def _normalize(name):
    return name.lower().replace('_', '').replace('-', '')


class Step:
    def __init__(self, cnf, run_id, name, script, dir_name, interpreter=None, short_name=None, paramln=None):
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


# noinspection PyAttributeOutsideInit
class Runner:
    def __init__(self, cnf, bcbio_final_dir, bcbio_cnf):
        self.final_dir = bcbio_final_dir
        self.cnf = cnf
        self.bcbio_cnf = bcbio_cnf

        hasher = hashlib.sha1(bcbio_final_dir)
        self.run_id = base64.urlsafe_b64encode(hasher.digest()[0:8])[:-1]

        self.threads = str(self.cnf.threads)
        self.qsub_runner = abspath(expanduser(cnf.qsub_runner))

        self.date_dirpath = join(bcbio_final_dir, bcbio_cnf.fc_date + '_' + bcbio_cnf.fc_name)
        if not verify_dir(self.date_dirpath):
            err('No project directory of format {fc_date}_{fc_name}, creating ' + self.date_dirpath)
        safe_mkdir(self.date_dirpath)

        self.log_dirpath = join(self.date_dirpath, 'log')
        safe_mkdir(self.log_dirpath)

        self.steps = Steps()
        self.vardict_steps = Steps()

        self.set_up_steps(cnf, self.run_id)

        def contains(x, xs):
            return _normalize(x) in [_normalize(y) for y in (xs or [])]

        self.steps.extend(
            [s for s in [
                self.varqc,
                self.varannotate,
                self.varfilter_all,
                self.varqc_after,
                self.varqc_summary,
                self.targetcov,
                self.ngscat,
                self.qualimap,
                self.targetcov_summary,
                self.ngscat_summary,
                self.qualimap_summary]
             if contains(s.name, cnf.steps)])

        self.vardict_steps.extend(
            [s for s in [
                self.vardict,
                self.testsomatic,
                self.var_to_vcf_somatic,
                self.varqc,
                self.varannotate,
                self.varfilter_all,
                self.varqc_after,
                self.varqc_summary]
             if contains(s.name, cnf.vardict_steps)])


    def set_up_steps(self, cnf, run_id):
        cnfs_line = ' --sys-cnf \'' + self.cnf.sys_cnf + '\' --run-cnf \'' + self.cnf.run_cnf + '\' '
        overwrite_line = {True: '-w', False: '--reuse'}.get(cnf.overwrite, '')
        spec_params = cnfs_line + ' -t ' + str(self.threads) + ' ' + overwrite_line + ' '

        self.varannotate = Step(cnf, run_id,
            name='VarAnnotate', short_name='va',
            interpreter='python',
            script='varannotate',
            dir_name='varAnnotate',
            paramln=spec_params + ' --vcf \'{vcf}\' {bam_cmdline} -o \'{output_dir}\' -s \'{sample}\' '
                                  '--work-dir \'' + join(cnf.work_dir, 'varAnnotate') + '_{sample}\''
        )
        self.varqc = Step(cnf, run_id,
            name='VarQC', short_name='vq',
            interpreter='python',
            script='varqc',
            dir_name='qc/varQC',
            paramln=spec_params + ' --vcf \'{vcf}\' -o \'{output_dir}\''
                    ' -s \'{sample}\' --work-dir \'' + join(cnf.work_dir, 'varQC') + '_{sample}\''
        )
        self.varqc_after = Step(cnf, run_id,
            name='VarQC_postVarFilter', short_name='vqa',
            interpreter='python',
            script='varqc',
            dir_name='qc/varQC_postVarFilter',
            paramln=spec_params + ' --vcf \'{vcf}\' -o \'{output_dir}\' -s \'{sample}\' '
                                  '--work-dir \'' + join(cnf.work_dir, 'varQC_postVarFilter') + '_{sample}\''
        )
        self.targetcov = Step(cnf, run_id,
            name='TargetCov', short_name='tc',
            interpreter='python',
            script='targetcov',
            dir_name='targetSeq',
            paramln=spec_params + ' --bam \'{bam}\' --bed \'{bed}\' -o \'{output_dir}\' '
                    '-s \'{sample}\' --work-dir \'' + join(cnf.work_dir, 'targetSeq') + '_{sample}\''
        )
        self.ngscat = Step(cnf, run_id,
            interpreter='python',
            script='ngscat',
            dir_name='qc/ngscat',
            name='ngsCAT', short_name='nc',
            paramln=spec_params + ' --bam \'{bam}\' --bed \'{bed}\' -o \'{output_dir}\' -s \'{sample}\' '
                                  '--saturation y --work-dir \'' + join(cnf.work_dir, 'ngscat') + '_{sample}\''
        )
        self.qualimap = Step(cnf, run_id,
            script='qualimap',
            dir_name='qc/qualimap',
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
            dir_name='varQC',
            paramln=cnfs_line + ' -o \'{output_dir}\' -d \''
                    + self.final_dir + '\' -s \'{samples}\' -n "' + self.varqc.dir_name + '" --vcf-suf ' + ','.join(all_suffixes) +
                    ' --work-dir \'' + join(cnf.work_dir, 'varQC_summary') + '\''
        )
        self.varfilter_all = Step(cnf, run_id,
            name='VarFilter', short_name='vfs',
            interpreter='python',
            script='varfilter_all',
            dir_name='varFilter',
            paramln=spec_params + ' -d \'' +
                    self.final_dir + '\' -s \'{samples}\' --vcf-suf ' + ','.join(all_suffixes) + ' ' +
                    '--work-dir \'' + join(cnf.work_dir, 'varFilter') + '\' --make_soft_links ' +
                    '--vcf-dir ' + self.varannotate.dir_name
        )
        self.targetcov_summary = Step(cnf, run_id,
            name='TargetCov_summary', short_name='tcs',
            interpreter='python',
            script='targetcov_summary',
            dir_name='targetSeq',
            paramln=cnfs_line + ' -o \'{output_dir}\' -d \''
                    + self.final_dir + '\' -s \'{samples}\' -n "' + self.targetcov.dir_name + '" --work-dir \'' +
                    join(cnf.work_dir, 'targetSeq_summary') + '\''
        )
        self.ngscat_summary = Step(cnf, run_id,
            name='ngsCAT_summary', short_name='ncs',
            interpreter='python',
            script='ngscat_summary',
            dir_name='ngscat',
            paramln=cnfs_line + ' -o \'{output_dir}\' -d \''
                    + self.final_dir + '\' -s \'{samples}\' -n "' + self.ngscat.dir_name + '" --work-dir \'' +
                    join(cnf.work_dir, 'ngscat_summary') + '\''
        )
        self.qualimap_summary = Step(cnf, run_id,
            name='QualiMap_summary', short_name='qms',
            interpreter='python',
            script='qualimap_summary',
            dir_name='qualimap',
            paramln=cnfs_line + ' -o \'{output_dir}\' -d \''
                    + self.final_dir + '\' -s \'{samples}\' -n "' + self.qualimap.dir_name + '" --work-dir \'' +
                    join(cnf.work_dir, 'qualimap_summary') + '\''
        )

        af_thr = str(cnf.variant_filtering.min_freq)
        self.vardict = Step(cnf, run_id,
            name='VarDict',
            interpreter='perl',
            script='vardict_pl',
            dir_name='VarDict',
            paramln=' -G ' + cnf.genome.seq + ' -f ' + af_thr + ' -N {tumor_name} -b \'{tumor_bam}|{normal_bam}\''
                    ' -z -F -c 1 -S 2 -E 3 -g 4 {bed} > {vars_txt}'
        )
        self.testsomatic = Step(cnf, run_id,
            name='TestSomatic',
            script='testsomatic_r',
            dir_name='VarDict',
            paramln=' < {vars_txt} > {somatic_vars_txt}',
        )
        self.var_to_vcf_somatic = Step(cnf, run_id,
            name='Var2Vcf_Somatic',
            interpreter='perl',
            script='var2vcf_somatic_pl',
            dir_name='VarDict',
            paramln=' -N \'{tumor_name}|{normal_name}\' -f ' +
                    af_thr + ' < {somatic_vars_txt} > {vardict_vcf}',
        )

    def step_output_dir_and_log_paths(self, step, sample_name, caller=None):
        if sample_name:
            base_output_dirpath = abspath(join(self.final_dir, sample_name))
        else:
            base_output_dirpath = abspath(self.date_dirpath)

        output_dirpath = join(base_output_dirpath, step.dir_name)

        log_fpath = join(self.log_dirpath,
             (step.dir_name + ('_' + sample_name if sample_name else '') +
                              ('_' + caller if caller else '')) + '.log')

        return output_dirpath, log_fpath


    def submit(self, step, sample_name='', suf=None, create_dir=True,
               out_fpath=None, wait_for_steps=list(), threads=None, **kwargs):

        output_dirpath, log_fpath = self.step_output_dir_and_log_paths(step, sample_name, suf)
        if not isdir(output_dirpath) and create_dir:
            safe_mkdir(output_dirpath)
        out_fpath = out_fpath or log_fpath

        if isfile(out_fpath):
            try:
                os.remove(out_fpath)
            except OSError:
                err('Warning: cannot remove log file ' + out_fpath + ', probably permission denied.')

        if log_fpath and isfile(log_fpath):
            try:
                os.remove(log_fpath)
            except OSError:
                err('Warning: cannot remove log file ' + out_fpath + ', probably permission denied.')

        safe_mkdir(dirname(log_fpath))
        safe_mkdir(dirname(out_fpath))

        tool_cmdline = get_tool_cmdline(self.cnf, step.interpreter, step.script)
        if not tool_cmdline: sys.exit(1)
        params = dict({'output_dir': output_dirpath}.items() + self.__dict__.items() + kwargs.items())
        cmdline = tool_cmdline + ' ' + step.param_line.format(**params)

        hold_jid_line = '-hold_jid ' + ','.join(wait_for_steps or ['_'])
        job_name = step.job_name(sample_name, suf)
        qsub = get_tool_cmdline(self.cnf, 'qsub')
        threads = str(threads or self.threads)
        queue = self.cnf.queue
        runner_script = self.qsub_runner
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

        if self.cnf.verbose: info()
        return output_dirpath

    def _qualimap_bed(self, bed_fpath):
        if 'QualiMap' in self.steps and bed_fpath:
            qualimap_bed_fpath = join(self.cnf.work_dir, 'tmp_qualimap.bed')

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
                        fields.append('-')

                    if len(fields) < 5:
                        fields.append('0')

                    if len(fields) < 6:
                        fields.append('+')

                    out.write('\t'.join(fields) + '\n')

            return qualimap_bed_fpath
        else:
            return bed_fpath

    def _sumbit_vardict(self, batches):
        for batch_name, batch in batches.items():
            normal_name, normal_bam_fpath = batch['normal']
            bed_fpath = batch['bed']
            for tumor_name, tumor_bam_fpath in batch['tumor'].items():
                output_dirpath, _ = self.step_output_dir_and_log_paths(self.vardict, tumor_name)
                vars_txt = join(output_dirpath, 'vardict.txt')

                if not verify_bam(tumor_bam_fpath):
                    sys.exit(1)

                if not file_exists(tumor_bam_fpath + '.bai'):
                    samtools = get_tool_cmdline(self.cnf, 'samtools')
                    cmdline = '{samtools} index {bam}'.format(samtools=samtools, bam=tumor_bam_fpath)
                    call(self.cnf, cmdline)

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
                        self.testsomatic, tumor_name, suf='testsomatic',
                        vars_txt=vars_txt,
                        somatic_vars_txt=somatic_vars_txt,
                        wait_for_steps=[self.vardict.job_name(tumor_name, 'vardict')])

                vardict_vcf = join(output_dirpath, 'somatic_variants-vardict_standalone.vcf')
                if self.var_to_vcf_somatic in self.vardict_steps:
                    self.submit(
                        self.var_to_vcf_somatic, tumor_name, suf='var2vcf',
                        tumor_name=tumor_name,
                        normal_name=normal_name,
                        somatic_vars_txt=somatic_vars_txt,
                        vardict_vcf=vardict_vcf,
                        wait_for_steps=[self.testsomatic.job_name(tumor_name, 'testsomatic')])

                self._process_vcf(
                    tumor_name, tumor_bam_fpath, vardict_vcf, 'vardict_standalone', steps=self.vardict_steps,
                    job_names_to_wait=[self.var_to_vcf_somatic.job_name(tumor_name, 'var2vcf')])

    def run(self):
        batches = defaultdict(dict)

        for sample_info in self.bcbio_cnf.details:
            sample = sample_info['description']

            if not (any(step in self.steps for step in
                        [self.targetcov, self.qualimap, self.ngscat,
                         self.varqc, self.varqc_after, self.varannotate]) or self.vardict_steps):
                continue

            info('Processing "' + sample + '"')
            if not self.cnf.verbose:
                info(ending='')

            sample_dirpath = join(self.final_dir, sample)
            if not verify_dir(sample_dirpath):
                sys.exit(1)

            safe_mkdir(join(sample_dirpath, 'qc'))

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
                    bed_fpath = self._qualimap_bed(bed_fpath)
            else:
                if not file_exists(bam_fpath):
                    bam_fpath = None
                if not file_exists(bed_fpath):
                    bed_fpath = None

            phenotype = None
            if 'metadata' in sample_info:
                phenotype = sample_info['metadata']['phenotype']

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
                vcf_fname = sample + '-' + variant_caller + '.vcf'
                # print 'vcf_fname = ' + vcf_fname
                vcf_fpath = join(sample_dirpath, vcf_fname)
                if not file_exists(vcf_fpath) and file_exists(vcf_fpath + '.gz'):
                    gz_vcf_fpath = vcf_fpath + '.gz'
                    gunzip = get_tool_cmdline(self.cnf, 'gunzip')
                    cmdline = '{gunzip} -c {gz_vcf_fpath}'.format(**locals())
                    call(self.cnf, cmdline, output_fpath=vcf_fpath)
                    info()

                var_dirpath = abspath(join(self.final_dir, sample, 'var'))
                # print 'creating var_dirpath = ' + var_dirpath
                safe_mkdir(var_dirpath)

                for fname in os.listdir(sample_dirpath):
                    # print '  listdir: vcf_fname ' + vcf_fname + ' ' + ('in ' if vcf_fname in fname else ' not in') + ' fname ' + fname
                    if vcf_fname in fname:
                        src_fpath = join(sample_dirpath, fname)
                        dst_fpath = join(var_dirpath, fname)
                        if exists(dst_fpath):
                            os.remove(dst_fpath)
                        safe_mkdir(var_dirpath)
                        info('Moving ' + src_fpath + ' to ' + dst_fpath)
                        os.rename(src_fpath, dst_fpath)

                vcf_fpath = join(var_dirpath, vcf_fname)

                if not file_exists(vcf_fpath):
                    if phenotype != 'normal':
                        err('No ' + vcf_fpath + ', skipping')
                        err()
                    continue

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

        if self.vardict_steps:
            self._sumbit_vardict(batches)

        all_variantcallers = set()
        for s_info in self.bcbio_cnf.details:
            all_variantcallers |= set(s_info['algorithm'].get('variantcaller')) or set()

        samples_fpath = abspath(join(self.cnf.work_dir, 'samples.txt'))
        #if not isfile(samples_fpath):
        with open(samples_fpath, 'w') as f:
            for sample_info in self.bcbio_cnf.details:
                sample = sample_info['description']
                f.write(sample + '\n')

        if self.varqc_summary in self.steps:
            self.submit(
                self.varqc_summary,
                wait_for_steps=[
                    self.varqc.job_name(d['description'], v)
                    for d in self.bcbio_cnf.details
                    for v in all_variantcallers
                    if self.varqc in self.steps],
                # threads=samples_num + 1,
                samples=samples_fpath)

        if self.targetcov_summary in self.steps:
            self.submit(
                self.targetcov_summary,
                wait_for_steps=[
                    self.targetcov.job_name(d['description'])
                    for d in self.bcbio_cnf.details
                    if self.targetcov in self.steps],
                # threads=samples_num + 1,
                samples=samples_fpath)

        if self.ngscat_summary in self.steps:
            self.submit(
                self.ngscat_summary,
                wait_for_steps=[
                    self.ngscat.job_name(d['description'])
                    for d in self.bcbio_cnf.details
                    if self.ngscat in self.steps],
                # threads=samples_num + 1,
                samples=samples_fpath)

        if self.qualimap_summary in self.steps:
            self.submit(
                self.qualimap_summary,
                wait_for_steps=[
                    self.qualimap.job_name(d['description'])
                    for d in self.bcbio_cnf.details
                    if self.qualimap in self.steps],
                # threads=samples_num + 1,
                samples=samples_fpath)

        if self.varfilter_all in self.steps:
            self.submit(
                self.varfilter_all,
                wait_for_steps=[
                    self.varannotate.job_name(d['description'], v)
                    for d in self.bcbio_cnf.details
                    for v in all_variantcallers
                    if self.varannotate in self.steps],
                samples=samples_fpath,
                create_dir=False,
                threads=len(batches))

        if not self.cnf.verbose:
            print ''
        if self.cnf.verbose:
            info('Done.')

    def _process_vcf(self, sample, bam_fpath, vcf_fpath, caller,
                     steps=None, job_names_to_wait=list()):
        steps = steps or self.steps

        if self.varqc in steps:
            self.submit(
                self.varqc, sample, suf=caller, vcf=vcf_fpath,
                sample=sample + '-' + caller, wait_for_steps=job_names_to_wait)

        bam_cmdline = '--bam ' + bam_fpath if bam_fpath else ''

        if self.varannotate in steps:
            self.submit(
                self.varannotate, sample, suf=caller, vcf=vcf_fpath,
                bam_cmdline=bam_cmdline, sample=sample + '-' + caller,
                wait_for_steps=job_names_to_wait)

        anno_dirpath, _ = self.step_output_dir_and_log_paths(self.varannotate, sample, caller=caller)
        annotated_vcf_fpath = join(anno_dirpath, basename(add_suffix(vcf_fpath, 'anno')))

        filter_dirpath = join(dirname(anno_dirpath), self.varfilter_all.dir_name)
        safe_mkdir(filter_dirpath)
        filtered_vcf_fpath = join(filter_dirpath, basename(add_suffix(annotated_vcf_fpath, 'filt')))

        if self.varqc_after in steps:
            self.submit(
                self.varqc_after, sample, suf=caller,
                wait_for_steps=([self.varfilter_all.job_name()]
                                 if self.varfilter_all in steps else []) + job_names_to_wait,
                vcf=filtered_vcf_fpath, sample=sample + '-' + caller)
