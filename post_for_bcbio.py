#!/usr/bin/env python
from optparse import OptionParser
import os
from os.path import join, expanduser, abspath, dirname

import sys
from source.bcbio_utils import file_exists, safe_mkdir
from source.config import Defaults, Config
from source.logger import critical, info, err
from source.main import check_system_resources
from source.ngscat.bed_file import verify_bed, verify_bam
from source.transaction import file_transaction
from source.utils import verify_dir, verify_file, get_tool_cmdline, tmpfile, make_tmpdir, make_tmpfile, call, tmpdir
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))


basic_runner = join(abspath(dirname(__file__)), 'basic_runner.sh')

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
        cnfs_line = ('--sys-cnf "' + self.cnf.sys_cnf_fpath + '"' +
             '--run-cnf "' + self.cnf.run_cnf_fpath + '"')

        overwrite_line = '-w' if self.cnf.overwrite else ''
        tools_base = 'python /group/ngs/bin/'
        spec_params = cnfs_line + ' -t ' + self.threads + ' ' + overwrite_line + ' '

        self.indel_filter = Step('IndelFilter',
            tools_base + 'InDelFilter.py "{vcf}" > "{filtered_vcf}"')

        self.varannotate = Step('VarAnnotate',
            tools_base + 'varannotate.py ' + spec_params +
            '--vcf "{vcf}" --bam "{bam}" -o "{output_dir}"')

        self.varqc = Step('VarQC',
            tools_base + 'varqc.py ' + spec_params +
            '--vcf "{vcf}" -o "{output_dir}"')

        self.targetcov = Step('TargetCov',
            tools_base + 'targetcov.py ' + spec_params +
            '--bam "{bam}" --bed "{bed}" -o "{output_dir}"')

        self.ngscat = Step('NGScat',
            tools_base + 'ngscat.py ' + spec_params +
            '--bam "{bam}" --bed "{bed}" -o "{output_dir}" --saturation y')

        self.qualimap = Step('QualiMap',
            tools_base + 'qualimap/qualimap ' + spec_params +
            'bamqc -nt ' + self.threads + ' --java-mem-size=24G -nr 5000 -bam {bam} '
            '-outdir ${output_dir} -gff {bed} -c -gd HUMAN')

        self.varqc_summary = Step('VarQC_summary',
            tools_base + 'varqc_summary.py {dir} {samples} ' + self.varqc.name +
            ' {vcf_suffix}'.format(**locals()))

        self.targetcov_summary = Step('TargetCov_summary',
            tools_base + 'targetcov_summary.py {dir} {samples} ' +
            self.targetcov.name)

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

        _, qualimap_bed_fpath = make_tmpfile(self.cnf, 'tmp_qualimap.bed')
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
            if self.indel_filter.name in self.steps:
                filtered_vcf_fpath = sample + self.suf + '.filtered.vcf'
                submit(self.indel_filter, sample, False,
                       vcf=vcf_fpath, filtered_vcf=filtered_vcf_fpath)

            annotated_vcf_fpath = filtered_vcf_fpath
            if self.varannotate.name in self.steps:
                annotated_vcf_fpath = sample + self.suf + '.anno.vcf'
                submit(self.varannotate, sample, True,
                       wait_for_steps=[sample + '_' + self.indel_filter.name],
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

            if self.qualimap.name in self.steps:
                submit(self.qualimap, sample, True,
                       bam=bam_fpath, bed=qualimap_bed_fpath)

        os.remove(qualimap_bed_fpath)

        if self.varqc.name in self.steps:
            submit(self.varqc_summary, wait_for_steps=[s + '_' + self.varqc.name for s in self.samples],
                   vcf_suffix=self.suf + '.anno')

        if self.targetcov.name in self.steps:
            submit(self.targetcov_summary, wait_for_steps=[s + '_' + self.targetcov.name for s in self.samples])


def main(args):
    script_name = __file__
    description = \
'''Usage: {script_name} -s <SAMPLES> -d <BCBIO FINAL DIR> -v <VCF SUFFIX> -b <BED FILE>

This script runs reporting suite on the bcbio final directory.
'''.format(**locals())

    parser = OptionParser(description=description)
    parser.add_option('-d', dest='dir', help='Path to bcbio-nextgen final directory (default is pwd)')
    parser.add_option('-s', dest='samples', help='List of samples (default is samples.txt in bcbio final directory)')
    parser.add_option('-b', dest='bed', help='BED file')
    parser.add_option('--vcf-suffix', dest='vcf_suffix', help='Suffix to file VCF file (mutect, ensembl, freebayes, etc)')
    parser.add_option('--qualimap', dest='qualimap', action='store_true', default=False, help='Run QualiMap in the end')

    parser.add_option('-v', dest='verbose', action='store_true', default=False, help='Verbose')
    parser.add_option('-t', dest='threads', type='int', default=4, help='Number of threads for each process')
    parser.add_option('-w', dest='overwrite', action='store_true', default=False, help='Overwrite existing results')
    parser.add_option('--sys-cnf', dest='sys_cnf', default=Defaults.sys_cnf_path, help='system configuration yaml with paths to external tools and genome resources (see default one %s)' % Defaults.sys_cnf_path)
    parser.add_option('--run-cnf', dest='run_cnf', default=Defaults.run_cnf_path, help='run configuration yaml (see default one %s)' % Defaults.run_cnf_path)

    (opts, args) = parser.parse_args()

    to_exit = False

    bcbio_final_dir = opts.dir
    if not bcbio_final_dir:
        info('The option -d is not specified, using cwd as bcbio final directory.')
        bcbio_final_dir = os.getcwd()
    if not verify_dir(bcbio_final_dir):
        to_exit = True
    else:
        bcbio_final_dir = abspath(expanduser(bcbio_final_dir))

    samples_fpath = opts.samples
    if not samples_fpath:
        info('The option -s is not specified, looking for samples.txt for the list of samples.')
        samples_fpath = join(bcbio_final_dir, 'samples.txt')
    if not verify_file(samples_fpath):
        to_exit = True
    else:
        samples_fpath = abspath(expanduser(samples_fpath))

    bed_fpath = opts.bed
    if not bed_fpath:
        err('BED file is not specified, use the -b option.')
        to_exit = True
    if not verify_bed(bed_fpath, 'BED file'):
        to_exit = True
    else:
        bed_fpath = abspath(expanduser(opts.bed))

    if not opts.vcf_suffix:
        err('VCF suffix is not specified. Please, use --vcf-suffix (i.e. --vcf-suffix mutect)')
        to_exit = True

    if to_exit:
        sys.exit(1)

    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)

    if opts.qualimap and 'QualiMap' not in cnf.steps:
        cnf.steps.append('QualiMap')

    check_system_resources(cnf, required=['qsub'])

    with tmpdir(cnf):
        runner = Runner(cnf, bcbio_final_dir, samples_fpath, bed_fpath, opts.vcf_suffix)

        runner.run()


if __name__ == '__main__':
    main(sys.argv[1:])









