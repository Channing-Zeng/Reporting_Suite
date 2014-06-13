#!/usr/bin/env python
from optparse import OptionParser
import os
from os.path import join, expanduser, abspath, dirname

import sys
from source.bcbio_utils import file_exists, safe_mkdir
from source.logger import critical
from source.main import load_configs, check_system_resources
from source.transaction import file_transaction
from source.utils import verify_bed, verify_dir, verify_file, get_tool_cmdline, iterate_file, verify_bam
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))


default_sys_config_path = 'system_info_Waltham.yaml'
default_run_config_path = 'run_info.yaml'
basic_runner = join(abspath(dirname(__file__)), 'basic_runner.sh')



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
    parser.add_option('--sys-cnf', dest='sys_cnf', default=default_sys_config_path, help='system configuration yaml with paths to external tools and genome resources (see default one %s)' % default_sys_config_path)
    parser.add_option('--run-cnf', dest='run_cnf', default=default_run_config_path, help='run configuration yaml (see default one %s)' % default_run_config_path)

    (opts, args) = parser.parse_args()

    bcbio_final_dir = abspath(expanduser(opts.dir)) or os.getcwd() if verify_dir(opts.dir) else sys.exit(1)
    samples_fpath = abspath(expanduser(opts.samples)) or join(bcbio_final_dir, 'samples.txt') if verify_file(opts.sample) else sys.exit(1)
    bed_fpath = abspath(expanduser(opts.bed)) if verify_bed(opts.bed) else sys.exit(1)

    cnf = load_configs(opts.sys_cnf, opts.run_cnf)

    check_system_resources(cnf, required=['qsub'])

    cnf['tmp_dir'] = cnf.get('tmp_base_dir') or os.getcwd()
    cnf['bcbio_final_dir'] = bcbio_final_dir

    cnfs_line = '--sys-cnf "' + cnf['sys_cnf'] + '" --run-cnf "' + cnf['run_cnf'] + '"'
    overwrite_line = '-w' if opts.overwrite else ''
    tools_base = 'python /group/ngs/bin/'
    spec_params = cnfs_line + ' -t ' + str(opts.threads) + ' ' + overwrite_line + ' '

    class Step():
        def __init__(self, name, cmdline):
            self.name = name
            self.cmdline = cmdline

    indel_filter = Step('IndelFilter',
        tools_base + 'InDelFilter.py "{vcf}" > "{filtered_vcf}"')

    varannotate = Step('VarAnnotate',
        tools_base + 'varannotate.py ' + spec_params +
        '--vcf "{vcf}" --bam "{bam}" -o "{output_dir}"')

    varqc = Step('VarQC',
        tools_base + 'varqc.py ' + spec_params +
        '--vcf "{vcf}" -o "{output_dir}"')

    targetcov = Step('TargetCov',
        tools_base + 'targetcov.py ' + spec_params +
        '--bam "{bam}" --bed "{bed}" -o "{output_dir}"')

    ngscat = Step('NGScat',
        tools_base + 'ngscat.py ' + spec_params +
        '--bam "{bam}" --bed "{bed}" -o "{output_dir}" --saturation y')

    qualimap = Step('QualiMap',
        tools_base + 'qualimap/qualimap ' + spec_params +
        'bamqc -nt ' + str(opts.threads) + ' --java-mem-size=24G -nr 5000 -bam {bam} '
        '-outdir ${output_dir} -gff {bed} -c -gd HUMAN')

    varqc_summary = Step('VarQC_summary',
        tools_base + 'varqc_summary.py {bcbio_final_dir} {samples_fpath} ' + varqc.name +
        ' {vcf_suffix}'.format(**locals()))

    targetcov_summary = Step('TargetCov_summary',
        tools_base + 'targetcov_summary.py {bcbio_final_dir} {samples_fpath} ' +
        targetcov.name)

    qualimap_bed_fpath = join(cnf['tmp_dir'], 'tmp_qualimap.bed')
    with open(qualimap_bed_fpath, 'w') as out, open(bed_fpath) as inn:
        for l in inn:
            ts = l.strip().split('\t')
            if len(ts) < 5:
                ts += ['0']
            if len(ts) < 6:
                ts += ['+']

            out.write('\t'.join(ts) + '\n')

    def submit(step, sample_name='', create_dir=False, wait_for_steps=list(), **kwargs):
        output_dirpath = bcbio_final_dir
        if create_dir:
            output_dirpath = join(bcbio_final_dir, sample_name, step.name)
            safe_mkdir(output_dirpath)

        log_fpath = join(output_dirpath, step.name + '.log')
        out_fpath = log_fpath

        hold_jid_line = '-hold_jid ' + ','.join(wait_for_steps) if wait_for_steps else ''

        job_name = sample_name + '_' + step.name
        threads = str(opts.threads)
        runner_script = basic_runner + ' ' + step.cmdline.format(
            dict({'output_dir': output_dirpath}.items() + kwargs.items()))

        qsub = get_tool_cmdline(cnf, 'qsub')

        qsub_cmdline = (
            '{qsub} {hold_jid_line} -pe smp {threads} -S /bin/bash -q batch.q '
            '-b y -o {out_fpath} -e {log_fpath} -N ${name} bash ${runner_script}'.format(
                **locals()))

    with open(samples_fpath) as sample_f:
        samples = [s.strip() for s in sample_f.readlines() if s and s.strip() and not s.startswith('#')]

    for sample in samples:
        sample_dirpath = join(bcbio_final_dir, sample)
        if not verify_dir(sample_dirpath): sys.exit(1)

        bam_fpath = join(sample_dirpath, sample + '-ready.bam')
        if not verify_bam(bam_fpath): sys.exit(1)

        vcf_fpath = join(sample_dirpath, sample + '-' + opts.vcf_suffix + '.vcf')
        if not file_exists(vcf_fpath) and file_exists(vcf_fpath + '.gz'):
            vcf_fpath += '.gz'
            # gzip = get_tool_cmdline(cnf, 'gunzip')
            # call(gunzip -c "${sample}${VCF_SUFFIX}.vcf.gz" > "${sample}${VCF_SUFFIX}.vcf"
        if not verify_file(vcf_fpath):
            sys.exit(1)


        filtered_vcf_fpath = vcf_fpath
        if indel_filter.name in cnf['steps']:
            filtered_vcf_fpath = sample + '-' + opts.vcf_suffix + '.filtered.vcf'
            submit(indel_filter, sample, False,
                   vcf=vcf_fpath, filtered_vcf=filtered_vcf_fpath)

        annotated_vcf_fpath = filtered_vcf_fpath
        if varannotate.name in cnf['steps']:
            annotated_vcf_fpath = sample + '-' + opts.vcf_suffix + '.anno.vcf'
            submit(varannotate, sample, True,
                   wait_for_steps=[sample + '_' + indel_filter.name],
                   vcf=filtered_vcf_fpath, bam=bam_fpath)

        if varqc.name in cnf['steps']:
            submit(varqc, sample, True,
                   wait_for_steps=[sample + '_' + varannotate.name],
                   vcf=filtered_vcf_fpath)

        if targetcov.name in cnf['steps']:
            submit(targetcov, sample, True,
                   bam=bam_fpath, bed=bed_fpath)

        if ngscat.name in cnf['steps']:
            submit(ngscat, sample, True,
                   bam=bam_fpath, bed=bed_fpath)

        if qualimap.name in cnf['steps']:
            submit(qualimap, sample, True,
                   bam=bam_fpath, bed=qualimap_bed_fpath)

    if varqc.name in cnf['steps']:
        submit(varqc_summary, wait_for_steps=[s + '_' + varqc.name for s in samples],
               vcf_suffix=opts.vcf_suffix + '.anno')

    if targetcov.name in cnf['steps']:
        submit(targetcov_summary, wait_for_steps=[s + '_' + targetcov.name for s in samples])















