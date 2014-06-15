import sys
from os.path import join, dirname, abspath

from source.utils_from_bcbio import file_exists, safe_mkdir
from source.logger import info
from source.ngscat.bed_file import verify_bam
from source.utils import verify_dir, verify_file, get_tool_cmdline, tmpfile, call


basic_dirpath = dirname(dirname(abspath(__file__)))


def run_on_bcbio_final_dir(cnf, bcbio_final_dir, samples_fpath,
                           bed_fpath, vcf_suffix):
    return Runner(cnf, bcbio_final_dir, samples_fpath,
                  bed_fpath, vcf_suffix).run()


class Step():
    def __init__(self, name, cmdline=None):
        self.name = name
        self.cmdline = cmdline


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

    def step(self, name, script=None, interpreter=None, param_line=None):
        if name not in self:
            return None

        cmd = get_tool_cmdline(self.cnf, interpreter, script or name)
        if not cmd:
            sys.exit(1)

        return Step(name, cmd + ' ' + param_line)


class Runner():
    def __init__(self, cnf, bcbio_final_dir, samples_fpath, bed_fpath, vcf_suffix):
        self.dir = bcbio_final_dir
        self.cnf = cnf
        self.bed = bed_fpath
        self.suf = '-' + vcf_suffix if vcf_suffix else ''
        self.threads = str(self.cnf.threads)
        self.steps = Steps(cnf, cnf.steps)
        self.qsub_runner = cnf.qsub_runner

        self.indel_filter = None
        self.varannotate = None
        self.varqc = None
        self.filter_variants = None
        self.targetcov = None
        self.ngscat = None
        self.qualimap = None
        self.targetcov_summary = None
        self.varqc_summary = None

        with open(samples_fpath) as sample_f:
            self.samples = [s.strip() for s in sample_f.readlines()
                            if s and s.strip() and not s.startswith('#')]

        self.set_up_steps(cnf)


    def set_up_steps(self, cnf):
        self.indel_filter = self.steps.step(
            name='InDelFilter',
            interpreter='python',
            param_line='\'{vcf}\'')

        cnfs_line = '--sys-cnf \'' + self.cnf.sys_cnf + '\' --run-cnf \'' + self.cnf.run_cnf + '\''
        overwrite_line = '-w' if self.cnf.overwrite else ''
        spec_params = cnfs_line + ' -t ' + self.threads + ' ' + overwrite_line + ' '

        self.varannotate = self.steps.step(
            name='VarAnnotate',
            script='varannotate.py',
            interpreter='python',
            param_line=spec_params + ' --vcf \'{vcf}\' --bam \'{bam}\' -o \'{output_dir}\'')
        self.varqc = self.steps.step(
            name='VarQC',
            script='varqc.py',
            interpreter='python',
            param_line=spec_params + ' --vcf \'{vcf}\' -o \'{output_dir}\'')
        self.filter_variants = self.steps.step(
            name='FilterVariants',
            script='filter_variants.py',
            interpreter='python',
            param_line=spec_params + ' --vcf \'{vcf}\'')
        self.targetcov = self.steps.step(
            name='TargetCov',
            script='targetcov.py',
            interpreter='python',
            param_line=spec_params + ' --bam \'{bam}\' --bed \'{bed}\' -o \'{output_dir}\'')
        self.ngscat = self.steps.step(
            name='NGScat',
            script='ngscat.py',
            interpreter='python',
            param_line=spec_params + ' --bam \'{bam}\' --bed \'{bed}\' -o \'{output_dir}\' --saturation y')
        self.qualimap = self.steps.step(
            name='QualiMap',
            param_line=' bamqc -nt ' + self.threads + ' --java-mem-size=24G -nr 5000 -bam \'{bam}\' -outdir \'{output_dir}\' -gff \'{bed}\' -c -gd HUMAN')

        if self.varqc:
            self.steps.append('VarQC_summary')
            self.varqc_summary = self.steps.step(
                name='VarQC_summary',
                script='varqc_summary.py',
                interpreter='python',
                param_line=' {dir} {samples} ' + self.varqc.name + ' {vcf_suffix}')
        if self.targetcov:
            self.steps.append('TargetCov_summary')
            self.targetcov_summary = self.steps.step(
                name='TargetCov_summary',
                script='targetcov_summary.py',
                interpreter='python',
                param_line=' {dir} {samples} ' + self.targetcov.name)

    def submit(self, step, sample_name='', create_dir=False, out_fpath=None,
               wait_for_steps=list(), **kwargs):
        if not step or step.name not in self.steps:
            return None

        output_dirpath = self.dir
        if sample_name:
            output_dirpath = join(output_dirpath, sample_name)
        if create_dir:
            output_dirpath = join(output_dirpath, step.name)
            safe_mkdir(output_dirpath)

        log_fpath = join(output_dirpath, step.name + '.log')
        out_fpath = out_fpath or log_fpath

        hold_jid_line = '-hold_jid ' + ','.join(wait_for_steps or ['_'])

        job_name = step.name + '_' + sample_name if sample_name else step.name

        params = dict({'output_dir': output_dirpath}.items() +
                      self.__dict__.items() + kwargs.items())
        runner_script = self.qsub_runner
        cmdline = step.cmdline.format(**params)
        qsub = get_tool_cmdline(self.cnf, 'qsub')
        threads = str(self.threads)
        qsub_cmdline = (
            '{qsub} -pe smp {threads} -S /bin/bash -q batch.q '
            '-j y -o {out_fpath} -e {log_fpath} {hold_jid_line} '
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
        with tmpfile(self.cnf, 'tmp_qualimap.bed') as qualimap_bed_fpath:
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
                if not verify_dir(sample_dirpath): sys.exit(1)

                bam_fpath = join(sample_dirpath, sample + '-ready.bam')
                if not verify_bam(bam_fpath): sys.exit(1)

                vcf_fpath = join(sample_dirpath, sample + self.suf + '.vcf')
                if not file_exists(vcf_fpath) and file_exists(vcf_fpath + '.gz'):
                    gz_vcf_fpath = vcf_fpath + '.gz'
                    gunzip = get_tool_cmdline(self.cnf, 'gunzip')
                    cmdline = '{gunzip} -c {gz_vcf_fpath}'.format(**locals())
                    call(self.cnf, cmdline, output_fpath=vcf_fpath)
                    info()
                if not verify_file(vcf_fpath):
                    sys.exit(1)

                indel_filtered_vcf_fpath = vcf_fpath
                if self.indel_filter:
                    indel_filtered_vcf_fpath = join(sample_dirpath, sample + self.suf + '.filt_indels.vcf')
                    self.submit(self.indel_filter, sample, create_dir=False,
                           out_fpath=indel_filtered_vcf_fpath, vcf=vcf_fpath)

                annotated_vcf_fpath = indel_filtered_vcf_fpath
                if self.varannotate:
                    anno_dirpath = self.submit(self.varannotate, sample, True,
                        wait_for_steps=[sample + '_' + self.indel_filter.name] if self.indel_filter else [],
                        vcf=indel_filtered_vcf_fpath, bam=bam_fpath)
                    annotated_vcf_fpath = join(anno_dirpath, sample + self.suf + '.anno.vcf')

                    self.submit(self.varqc, sample, True, vcf=annotated_vcf_fpath,
                                wait_for_steps=[sample + '_' + self.varannotate.name])

                    if self.filter_variants:
                        filtered_vcf_fpath = join(anno_dirpath, sample + self.suf + '.anno.filt.vcf')
                        self.submit(self.filter_variants, sample, False,
                            wait_for_steps=[sample + '_' + self.varannotate.name],
                            out_fpath=filtered_vcf_fpath, vcf=annotated_vcf_fpath)

                self.submit(self.targetcov, sample, True, bam=bam_fpath, bed=self.bed)
                self.submit(self.ngscat, sample, True, bam=bam_fpath, bed=self.bed)
                self.submit(self.qualimap, sample, True, bam=bam_fpath, bed=qualimap_bed_fpath)

                if self.cnf.verbose:
                    info('-' * 70)
                else:
                    print ''
                    info()

        if not self.cnf.verbose:
            info('', ending='')

        samples = ','.join(self.samples)
        self.submit(
            self.varqc_summary,
            wait_for_steps=[s + '_' + self.varqc.name for s in self.samples],
            vcf_suffix=self.suf + '.anno',
            samples=samples)

        self.submit(
            self.targetcov_summary,
            wait_for_steps=[s + '_' + self.targetcov.name for s in self.samples],
            samples=samples)

        if not self.cnf.verbose:
            print ''
        if self.cnf.verbose:
            info('Done.')
