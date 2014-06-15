#!/usr/bin/env python
import sys
from optparse import OptionParser
from os.path import join, abspath, dirname, pardir

from source.bcbio_utils import file_exists, safe_mkdir
from source.config import Defaults, Config
from source.logger import info
from source.main import check_system_resources, check_inputs, check_keys
from source.ngscat.bed_file import verify_bam
from source.utils import verify_dir, verify_file, get_tool_cmdline, tmpfile, call, tmpdir


if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

cur_dirpath = dirname(abspath(__file__))
basic_runner = join(cur_dirpath, 'sh_runners', 'basic_runner.sh')


class Step():
    def __init__(self, name, cmdline):
        self.name = name
        self.cmdline = cmdline


class Runner():
    def __init__(self, cnf, bcbio_final_dir, samples_fpath, bed_fpath, vcf_suffix):
        self.dir = bcbio_final_dir
        self.cnf = cnf
        self.bed = bed_fpath
        self.suf = '-' + vcf_suffix if vcf_suffix else ''
        self.threads = str(self.cnf.threads)
        self.steps = cnf.steps

        with open(samples_fpath) as sample_f:
            self.samples = [s.strip() for s in sample_f.readlines()
                            if s and s.strip() and not s.startswith('#')]

        # SET UP STEPS
        python = get_tool_cmdline(cnf, 'python')

        self.indel_filter = None
        if 'IndelFilter' in self.steps:
            indel_filter_cmd = get_tool_cmdline(cnf, 'indelfilter')
            if not indel_filter_cmd:
                sys.exit(1)
            indel_filter_cmd = python + ' ' + indel_filter_cmd
            self.indel_filter = Step('IndelFilter',
                indel_filter_cmd + ' "{vcf}" > "{filtered_vcf}"')

        cnfs_line = '--sys-cnf "' + self.cnf.sys_cnf + '" --run-cnf "' + self.cnf.run_cnf + '"'
        overwrite_line = '-w' if self.cnf.overwrite else ''
        spec_params = cnfs_line + ' -t ' + self.threads + ' ' + overwrite_line + ' '
        my_script = python + ' ' + join(cur_dirpath, '%s') + ' ' + spec_params + ' '

        self.varannotate = Step(
            'VarAnnotate',
            my_script % 'varannotate.py --vcf "{vcf}" --bam "{bam}" -o "{output_dir}"')

        self.varqc = Step(
            'VarQC',
            my_script % 'varqc.py --vcf "{vcf}" -o "{output_dir}"')

        self.targetcov = Step(
            'TargetCov',
            my_script % 'targetcov.py --bam "{bam}" --bed "{bed}" -o "{output_dir}"')

        self.ngscat = Step(
            'NGScat',
            my_script % 'ngscat.py --bam "{bam}" --bed "{bed}" -o "{output_dir}" --saturation y')

        self.qualimap = None
        if 'QualiMap' in self.steps:
            qualimap_cmd = get_tool_cmdline(cnf, 'qualimap')
            if not qualimap_cmd:
                sys.exit(1)
            self.qualimap = Step(
                'QualiMap',
                qualimap_cmd + ' bamqc -nt ' + self.threads + ' --java-mem-size=24G -nr 5000 -bam {bam} -outdir ${output_dir} -gff {bed} -c -gd HUMAN')

        self.varqc_summary = Step(
            'VarQC_summary',
            my_script % 'varqc_summary.py {dir} {samples} ' + self.varqc.name + ' {vcf_suffix}')

        self.targetcov_summary = Step(
            'TargetCov_summary',
            my_script % 'targetcov_summary.py {dir} {samples} ' + self.targetcov.name)


    def run(self):
        def submit(step, sample_name='', create_dir=False, wait_for_steps=list(), **kwargs):
            output_dirpath = self.dir
            if create_dir:
                output_dirpath = join(self.dir, sample_name, step.name)
                safe_mkdir(output_dirpath)

            log_fpath = join(output_dirpath, step.name + '.log')
            out_fpath = log_fpath

            hold_jid_line = '-hold_jid ' + ','.join(wait_for_steps) if wait_for_steps else ''

            job_name = sample_name + '_' + step.name

            params = dict({'output_dir': output_dirpath}.items() + self.__dict__.items() + kwargs.items())
            runner_script = basic_runner + ' ' + step.cmdline.format(**params)

            qsub = get_tool_cmdline(self.cnf, 'qsub')
            threads = self.threads
            qsub_cmdline = (
                '{qsub} {hold_jid_line} -pe smp {threads} -S /bin/bash -q batch.q '
                '-b y -o {out_fpath} -e {log_fpath} -N ${job_name} bash ${runner_script}'.format(**locals()))

            call(self.cnf, qsub_cmdline)

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
                sample_dirpath = join(self.dir, sample)
                if not verify_dir(sample_dirpath): sys.exit(1)

                bam_fpath = join(sample_dirpath, sample + '-ready.bam')
                if not verify_bam(bam_fpath): sys.exit(1)

                vcf_fpath = join(sample_dirpath, sample + self.suf + '.vcf')
                if not file_exists(vcf_fpath) and file_exists(vcf_fpath + '.gz'):
                    vcf_fpath += '.gz'
                    # gzip = get_tool_cmdline(cnf, 'gunzip')
                    # call(gunzip -c "${sample}${VCF_SUFFIX}.vcf.gz" > "${sample}${VCF_SUFFIX}.vcf"
                if not verify_file(vcf_fpath):
                    sys.exit(1)

                filtered_vcf_fpath = vcf_fpath
                if self.indel_filter:
                    filtered_vcf_fpath = sample + self.suf + '.filtered.vcf'
                    submit(self.indel_filter, sample, False,
                           vcf=vcf_fpath, filtered_vcf=filtered_vcf_fpath)

                annotated_vcf_fpath = filtered_vcf_fpath
                if self.varannotate.name in self.steps:
                    annotated_vcf_fpath = sample + self.suf + '.anno.vcf'
                    submit(self.varannotate, sample, True,
                           wait_for_steps=[sample + '_' + self.indel_filter.name] if self.indel_filter else [],
                           vcf=filtered_vcf_fpath, bam=bam_fpath)

                if self.varqc.name in self.steps:
                    submit(self.varqc, sample, True,
                           wait_for_steps=[sample + '_' + self.varannotate.name],
                           vcf=annotated_vcf_fpath)

                if self.targetcov.name in self.steps:
                    submit(self.targetcov, sample, True,
                           bam=bam_fpath, bed=self.bed)

                if self.ngscat.name in self.steps:
                    submit(self.ngscat, sample, True,
                           bam=bam_fpath, bed=self.bed)

                if self.qualimap:
                    submit(self.qualimap, sample, True,
                           bam=bam_fpath, bed=qualimap_bed_fpath)

        if self.varqc.name in self.steps:
            submit(self.varqc_summary, wait_for_steps=[s + '_' + self.varqc.name for s in self.samples],
                   vcf_suffix=self.suf + '.anno')

        if self.targetcov.name in self.steps:
            submit(self.targetcov_summary, wait_for_steps=[s + '_' + self.targetcov.name for s in self.samples])


def main(args):
    description = 'This script runs reporting suite on the bcbio final directory.'

    parser = OptionParser(description=description)
    parser.add_option('-d', dest='bcbio_final_dir', default=Defaults.bcbio_final_dir, help='Path to bcbio-nextgen final directory (default is pwd)')
    parser.add_option('-s', dest='samples', help='List of samples (default is samples.txt in bcbio final directory)')
    parser.add_option('-b', dest='bed', help='BED file')
    parser.add_option('--vcf-suffix', dest='vcf_suffix', help='Suffix to choose VCF file s(mutect, ensembl, freebayes, etc)')
    parser.add_option('--qualimap', dest='qualimap', action='store_true', default=Defaults.qualimap, help='Run QualiMap in the end')

    parser.add_option('-v', dest='verbose', action='store_true', default=Defaults.verbose, help='Verbose')
    parser.add_option('-t', dest='threads', type='int', default=Defaults.threads, help='Number of threads for each process')
    parser.add_option('-w', dest='overwrite', action='store_true', default=Defaults.overwrite, help='Overwrite existing results')
    parser.add_option('--sys-cnf', dest='sys_cnf', default=Defaults.sys_cnf, help='system configuration yaml with paths to external tools and genome resources (see default one %s)' % Defaults.sys_cnf)
    parser.add_option('--run-cnf', dest='run_cnf', default=Defaults.run_cnf, help='run configuration yaml (see default one %s)' % Defaults.run_cnf)

    (opts, args) = parser.parse_args()
    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)

    if not cnf.samples:
        cnf.samples = join(cnf.bcbio_final_dir, 'samples.txt')

    info('BCBio "final" dir: ' + cnf.bcbio_final_dir + ' (set with -d)')
    info('Samples: ' + cnf.samples + ' (set with -s)')

    if not check_keys(cnf, ['bcbio_final_dir', 'samples', 'bed', 'vcf_suffix']):
        parser.print_help()
        sys.exit(1)

    if not check_inputs(cnf, ['bcbio']):
        sys.exit(1)

    if not check_inputs(cnf, file_keys=['samples', 'bed'], dir_keys=['bcbio_final_dir']):
        sys.exit(1)

    info('Capture/amplicons BED file: ' + cnf.bed + ' (set with -b)')
    info('Suffix to choose VCF files: ' + cnf.vcf_suffix + ' (set with --vcf-suffix)')

    if opts.qualimap and 'QualiMap' not in cnf.steps:
        cnf.steps.append('QualiMap')

    check_system_resources(cnf, required=['qsub'])

    with tmpdir(cnf):
        runner = Runner(cnf, cnf.bcbio_final_dir, cnf.samples, cnf.bed, cnf.vcf_suffix)
        runner.run()


if __name__ == '__main__':
    main(sys.argv[1:])









