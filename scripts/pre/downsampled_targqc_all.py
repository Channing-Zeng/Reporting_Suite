#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import shutil

import bcbio_postproc

import os
import sys
import datetime
from optparse import OptionParser
from os.path import join, isfile, basename, isdir, exists, dirname, splitext, islink
from collections import OrderedDict, namedtuple
import subprocess
import traceback

from ext_modules.joblib import Parallel, delayed
import source
from source.calling_process import call
from source.fastqc.fastq_utils import downsample
from source.targetcov.bam_and_bed_utils import index_bam, markdup_bam, verify_bam
from source.tools_from_cnf import get_system_path, get_script_cmdline
from source.config import Config, CallCnf
from source import logger
from source.logger import info, critical, err, is_local, warn, send_email
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug, check_system_resources, determine_sys_cnf, determine_run_cnf, \
    check_genome_resources, set_up_log
from source.file_utils import safe_mkdir, verify_dir, verify_file, adjust_path, \
    add_suffix, file_transaction, splitext_plus
from source.utils import get_ext_tools_dirname
from targqc import find_fastq_pairs

NGS_WEBSERVER_PREPROC_DIR = '/opt/lampp/htdocs/reports'
if is_local():
    NGS_WEBSERVER_PREPROC_DIR = '/Users/vlad/Sites/reports'


def proc_opts():
    parser = OptionParser()
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser)
    parser.add_option('--expose-only', dest='expose_to_ngs_server_only', action='store_true', default=False, help='Only add project to the webserver')
    parser.add_option('--no-expose', dest='expose', action='store_false', default=True, help='Do not expose the reports')
    parser.add_option('-o', dest='output_dir')
    parser.add_option('--bed', dest='bed', help='BED file to run targetSeq and Seq2C analysis on.')
    parser.add_option('--downsample-to', dest='downsample_to', type='int')

    (opts, args) = parser.parse_args()
    logger.is_debug = opts.debug

    if len(args) < 1:
        critical('Usage: ' + __file__ + ' *.fq.gz -o output_dir')
    # if len(args) < 2:
    #     info('No dataset path specified, assuming it is the current working directory')
    #     dataset_dirpath = adjust_path(os.getcwd())
    #     jira_url = args[0]

    fastq_fpaths = [verify_file(fpath) for fpath in args]
    fastq_fpaths = [fpath for fpath in fastq_fpaths if fpath]
    info(str(len(fastq_fpaths)) + ' fastq files')

    run_cnf = determine_run_cnf(opts)
    cnf = Config(opts.__dict__, determine_sys_cnf(opts), run_cnf)

    cnf.output_dir = adjust_path(cnf.output_dir)
    info('Writing to ' + str(cnf.output_dir))

    cnf.project_name = cnf.project_name or 'preproc'

    if cnf.work_dir:
        cnf.debug = True
    else:
        all_work_dir = join(cnf.output_dir, 'work')
        safe_mkdir(all_work_dir)

        latest_fpath = join(all_work_dir, 'latest')

        if cnf.reuse_intermediate:
            cnf.work_dir = latest_fpath
        else:
            cnf.work_dir = join(all_work_dir, datetime.datetime.now().strftime("%Y-%b-%d_%H-%M"))
            if islink(latest_fpath):
                os.remove(latest_fpath)
            if isdir(latest_fpath):
                shutil.rmtree(latest_fpath)
            if not exists(latest_fpath):
                os.symlink(basename(cnf.work_dir), latest_fpath)

    cnf.work_dir = adjust_path(cnf.work_dir)
    safe_mkdir(cnf.work_dir)
    cnf.log_dir = join(cnf.work_dir, 'log')
    safe_mkdir(cnf.log_dir)
    set_up_log(cnf)
    try:
        subprocess.call(['chmod', '-R', '777', cnf.work_dir])
    except OSError:
        err(traceback.format_exc())
        pass

    if cnf.samplesheet:
        cnf.samplesheet = verify_file(cnf.samplesheet, is_critical=True)

    info(' '.join(sys.argv))
    info()
    info('Created a temporary working directory: ' + cnf.work_dir)

    if cnf.project_name:
        info('Project name: ' + cnf.project_name)

    if cnf.samplesheet:
        info('Using custom sample sheet ' + cnf.samplesheet)

    check_genome_resources(cnf)
    check_system_resources(cnf, optional=['fastq'])

    return cnf, cnf.output_dir, fastq_fpaths


def main():
    cnf, output_dir, fastq_fpaths = proc_opts()

    targqc_dirpath = output_dir

    fastqs_by_sample = find_fastq_pairs(fastq_fpaths)
    samples = []
    for sname, (l, r) in fastqs_by_sample.items():
        s = source.TargQC_Sample(sname, join(cnf.output_dir, sname))
        s.l_fpath = l
        s.r_fpath = r
        samples.append(s)

    threads = len(samples)
    info('Found ' + str(len(samples)) + ' samples.')
    if len(samples) == 0:
        critical('ERROR: No fastq pairs found.')
    info()

    # samples = [source.TargQC_Sample(
    #     s.name,
    #     dirpath=join(targqc_dirpath, s.name),
    #     bed=cnf.bed) for s in fastq_fpaths]

    if cnf.downsample_to == 0:
        lefts = [s.l_fpath for s in samples]
        rights = [s.r_fpath for s in samples]
    else:
        if cnf.downsample_to is None:
            downsample_to = int(5e5)
        else:
            downsample_to = cnf.downsample_to

        info('Downsampling the reads to ' + str(downsample_to))
        lefts, rights = downsample_fastq(cnf, samples, downsample_to)

    bam_by_sample = OrderedDict()
    sambamba = get_system_path(cnf, join(get_ext_tools_dirname(), 'sambamba'), is_critical=True)
    bwa = get_system_path(cnf, 'bwa')
    bammarkduplicates = get_system_path(cnf, 'bammarkduplicates')
    if sambamba and bwa and bammarkduplicates:
        info()
        info('Aligning reads to the reference')
        bam_fpaths = Parallel(n_jobs=threads)(delayed(align)(CallCnf(cnf.__dict__), s, l, r,
            sambamba,
            bwa,
            bammarkduplicates,
            cnf.genome.bwa,
            cnf.is_pcr) for s, l, r in zip(samples, lefts, rights))
        for sample, bam_fpath in zip(samples, bam_fpaths):
            if verify_bam(bam_fpath):
                bam_by_sample[sample.name] = bam_fpath
            else:
                err('Sample ' + sample + ' was not aligned successfully.')
        if not bam_by_sample:
            err('ERROR: No sample was alined.')
        else:
            info()
            cnf.work_dir = join(cnf.work_dir, source.targqc_name)
            safe_mkdir(cnf.work_dir)
            info('Making TargQC reports for BAMs from reads')
            safe_mkdir(targqc_dirpath)
            run_targqc(cnf, bam_by_sample, cnf.bed, targqc_dirpath)
            cnf.work_dir = dirname(cnf.work_dir)
            info('Done TargQC')
    info()
    info('*' * 70)
    # if not cnf.debug and cnf.work_dir:
    #     try:
    #         shutil.rmtree(cnf.work_dir)
    #     except OSError:
    #         err('Can\'t remove work directory ' + cnf.work_dir + ', please, remove it manually.')


def downsample_fastq(cnf, samples, downsample_to=5e5):
    info('Downsampling reads to ' + str(int(downsample_to)) + ' pairs')
    # lefts, rights = [], []
    # for s in samples:
    #     info('Downsampling ' + s.name)
    #     l, r = downsample(cnf, s.l_fpath, s.r_fpath, N=downsample_to, output_dir=cnf.work_dir, suffix='subset')
    #     lefts.append(l)
    #     rights.append(r)

    fastqs = Parallel(n_jobs=len(samples)) \
        (delayed(downsample)(CallCnf(cnf.__dict__), s.name, s.l_fpath, s.r_fpath, N=downsample_to,
                             output_dir=cnf.work_dir, suffix='subset') \
            for s in samples)
    lefts = [l for l, r in fastqs]
    rights = [r for l, r in fastqs]
    return lefts, rights

    # downsample_script = get_script_cmdline(cnf, 'python', join('scripts', 'pre', 'downsample_fastq.py'))
    # cmdl = '{downsample_script} --sys-cnf {cnf.sys_cnf} --run-cnf {cnf.run_cnf} ' \
    #        '--downsample-to {downsample_to} -o {cnf.work_dir} --suffix subset '.format(**locals())
    # if cnf.reuse:
    #     cmdl += ' --reuse'
    # js = []
    # for s in samples:
    #     s_cmdl = cmdl + ' --sample ' + s.name + ' -1 ' + s.l_fpath
    #     l_subset_fpath, r_subset_fpath = None, None
    #     if s.r_fpath:
    #         s_cmdl += ' -2 ' + s.r_fpath
    #         r_subset_fpath = join(cnf.work_dir, add_suffix(basename(s.l_fpath), 'subset'))
    #     l_subset_fpath = join(cnf.work_dir, add_suffix(basename(s.l_fpath), 'subset'))
    #     j = submit_job(cnf, s_cmdl, 'downsample_' + s.name,
    #         l_subset_fpath=l_subset_fpath, r_subset_fpath=r_subset_fpath)
    #     js.append(j)
    # js = wait_for_jobs(cnf, js)
    #
    # lefts = [j.l_subset_fpath for j in js if j.l_subset_fpath]
    # rights = [j.r_subset_fpath for j in js if j.r_subset_fpath]
    # return lefts, rights


def align(cnf, sample, l_fpath, r_fpath, sambamba, bwa, bammarkduplicates, bwa_prefix, is_pcr=False):
    sam_fpath = join(cnf.work_dir, sample.name + '_downsampled.sam')
    bam_fpath = splitext(sam_fpath)[0] + '.bam'
    sorted_bam_fpath = add_suffix(bam_fpath, 'sorted')

    bwa_cmdline = '{bwa} mem {bwa_prefix} {l_fpath} {r_fpath} '.format(**locals())
    res = call(cnf, bwa_cmdline, output_fpath=sam_fpath, exit_on_error=False)
    if not res:
        return None

    cmdline = '{sambamba} view -t {cnf.threads} -S -f bam {sam_fpath}'.format(**locals())
    call(cnf, cmdline, output_fpath=bam_fpath)

    prefix = splitext(sorted_bam_fpath)[0]
    cmdline = '{sambamba} sort -t {cnf.threads} {bam_fpath} -o {sorted_bam_fpath}'.format(**locals())
    call(cnf, cmdline, output_fpath=sorted_bam_fpath, stdout_to_outputfile=False)

    if not is_pcr:
        markdup_bam_fpath = markdup_bam(cnf, sorted_bam_fpath, bammarkduplicates)
        if markdup_bam_fpath:
            sorted_bam_fpath = markdup_bam_fpath

    index_bam(cnf, sorted_bam_fpath, sambamba=sambamba)
    return sorted_bam_fpath


def run_targqc(cnf, bam_by_sample, bed_fpath, output_dirpath):
    info('Running TargQC for downsampled BAMs')

    targqc = get_script_cmdline(cnf, 'python', 'targqc.py', is_critical=True)
    targqc_work_dir = join(cnf.work_dir, 'TargQC')
    targqc_log_dir = join(cnf.log_dir, 'TargQC')
    safe_mkdir(targqc_work_dir)
    safe_mkdir(targqc_log_dir)
    bed_cmdl = ''
    if bed_fpath:
        bed_cmdl = '--bed ' + bed_fpath
    bam_cmdl = ' '.join(bam_fpath + ',' + sname for sname, bam_fpath in bam_by_sample.items())
    cmdl = '{targqc} --sys-cnf {cnf.sys_cnf} {bam_cmdl} {bed_cmdl} ' \
           '--work-dir {targqc_work_dir} --log-dir {targqc_log_dir} --project-name {cnf.project_name} ' \
           '-o {output_dirpath} --genome {cnf.genome.name}'.format(**locals())
    if cnf.reuse_intermediate:
        cmdl += ' --reuse'
    call(cnf, cmdl)

    # samples = [TargQCSample(
    #     s.name,
    #     output_dir=join(targqc_dirpath, s.name),
    #     bed=cnf.bed,
    #     bam=bam_by_sample[s.name])
    #            for s in samples]

    # Parallel(n_jobs=threads)(delayed(make_targetseq_reports(
    #     CallCnf(cnf.__dict__), sample.dirpath, sample,
    #     sample.bam, exons_bed, exons_no_genes_bed, target_bed
    # )(CallCnf(cnf.__dict__), sample) for sample in samples))
    #
    # return summarize_targqc(cnf, 1, targqc_dirpath, samples, bed_fpath, exons_bed)


if __name__ == '__main__':
    main()
