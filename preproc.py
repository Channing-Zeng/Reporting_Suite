#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc

import os
import sys
import datetime
from optparse import OptionParser
from os.path import join, isfile, basename, isdir, exists, dirname, splitext, islink
from collections import OrderedDict, namedtuple
import subprocess
import traceback

from joblib import Parallel, delayed
import source
from source.calling_process import call
from source.fastqc.fastq_utils import downsample
from source.fastqc.summarize_fastqc import write_fastqc_combo_report
from source.jira_utils import retrieve_jira_info
from source.preproc.dataset_structure import DatasetStructure
from source.bcbio.project_level_report import make_project_level_report
from source.qsub_utils import submit_job, wait_for_jobs
from source.targetcov.bam_and_bed_utils import index_bam, markdup_bam
from source.tools_from_cnf import get_system_path, get_script_cmdline
from source.config import Config, CallCnf
from source.logger import info, critical, err, is_local, warn, send_email
from source.utils import is_az
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug, check_system_resources, determine_sys_cnf, determine_run_cnf, \
    check_genome_resources, set_up_log
from source.file_utils import safe_mkdir, verify_dir, verify_file, adjust_path, \
    add_suffix, file_transaction
from source.webserver.exposing import sync_with_ngs_server

NGS_WEBSERVER_PREPROC_DIR = '/opt/lampp/htdocs/reports'
if is_local():
    NGS_WEBSERVER_PREPROC_DIR = '/Users/vlad/Sites/reports'

def proc_opts():
    usage = 'Usage: ' + __file__ + ' <DATASET_DIR_LOCATION> <JIRA_URL> [--bed BED] [--project-name STR] [--expose-only]\n' \
                                   '\tNote that DATASET_DIR_LOCATION can be either the root of a Illumina run ' \
                                   '(in witch case multiple projects may be detected and processed, or ' \
                                   'a specific project in <root>/Unalign/<Project>'

    description = 'This script runs preprocessing. ' + usage

    parser = OptionParser(description=description)
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser)
    parser.add_option('-j', '--jira', dest='jira', help='JIRA case path (goes to the ngs-website)')
    parser.add_option('--bed', dest='bed', help='BED file (used for downsampled TargQC reports)')
    parser.add_option('--samplesheet', dest='samplesheet', help='Sample sheet (default is located under the dataset root as SampleSheet.csv')
    # parser.add_option('--datahub-path', dest='datahub_path', help='DataHub directory path to upload final MAFs and CNV (can be remote).')
    # parser.add_option('--reporter', dest='reporter', help='Reporter name (goes to the ngs-website).')

    # parser.add_option('-e', '--expose', dest='expose_to_ngs_server', action='store_true', default=True, help='Add project to the webserver')
    parser.add_option('--expose-only', dest='expose_to_ngs_server_only', action='store_true', default=False, help='Only add project to the webserver')
    parser.add_option('--no-expose', dest='expose', action='store_false', default=True, help='Do not expose the reports')
    parser.add_option('--targqc', dest='targqc', action='store_true', default=True, help='')
    parser.add_option('--no-targqc', dest='targqc', action='store_false', default=True, help='')
    parser.add_option('--metamapping', dest='metamapping', action='store_true', default=False, help='')
    parser.add_option('--no-metamapping', dest='metamapping', action='store_false', default=False, help='')
    parser.add_option('--fastqc', dest='fastqc', action='store_true', default=True, help='')
    parser.add_option('--no-fastqc', dest='fastqc', action='store_false', default=True, help='')
    parser.add_option('--pcr', dest='is_pcr', action='store_true', default=False, help='')

    (opts, args) = parser.parse_args()
    if len(args) < 1:
        critical(usage)
    # if len(args) < 2:
    #     info('No dataset path specified, assuming it is the current working directory')
    #     dataset_dirpath = adjust_path(os.getcwd())
    #     jira_url = args[0]
    dataset_dirpath = verify_dir(args[0], is_critical=True, description='Dataset directory')  # /ngs/oncology/datasets/hiseq/150521_D00443_0159_AHK2KTADXX

    jira_url = ''
    if len(args) > 1:
        jira_url = args[1]

    run_cnf = determine_run_cnf(opts, is_wgs=not opts.__dict__.get('bed'))
    cnf = Config(opts.__dict__, determine_sys_cnf(opts), run_cnf)

    if cnf.work_dir:
        cnf.debug = True
    else:
        all_work_dir = join(dataset_dirpath, 'work')
        safe_mkdir(all_work_dir)

        latest_fpath = join(all_work_dir, 'latest')

        if cnf.reuse_intermediate:
            cnf.work_dir = latest_fpath
        else:
            cnf.work_dir = join(all_work_dir, datetime.datetime.now().strftime("%Y-%b-%d_%H-%M"))
            if exists(latest_fpath):
                os.remove(latest_fpath)
            os.symlink(cnf.work_dir, latest_fpath)

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

    info(' '.join(sys.argv))
    info()
    info('Created a temporary working directory: ' + cnf.work_dir)

    # check_genome_resources(cnf)
    check_system_resources(cnf, optional=['fastq'])
    check_genome_resources(cnf)

    if cnf.project_name:
        info('Project name: ' + cnf.project_name)

    if cnf.samplesheet:
        cnf.samplesheet = verify_file(cnf.samplesheet, is_critical=True)
        info('Using custom sample sheet ' + cnf.samplesheet)

    return cnf, dataset_dirpath, jira_url


def main():
    cnf, project_dirpath, jira_url = proc_opts()

    jira_case = None
    if is_az() and jira_url:
        info('Getting info from JIRA...')
        jira_case = retrieve_jira_info(jira_url)
        if not cnf.project_name:
            cnf.project_name = jira_case.project_name
            info('Setting project name from JIRA: ' + cnf.project_name)
    elif not cnf.project_name:
        critical('Cannot parse JIRA url ' + str(jira_url) + ', and --project-name is not specified. Please, provide a project name.')
    cnf.project_name = cnf.project_name.replace(' ', '_')
    info('Fixed final project name: ' + cnf.project_name)

    info()
    info('*' * 60)
    ds = DatasetStructure.create(project_dirpath, cnf.project_name, cnf.samplesheet)
    if not ds.project_by_name:
        critical('No projects found')
    for project in ds.project_by_name.values():
        if not project.sample_by_name:
            critical('No samples for project ' + project.name + ' found')

    for project in ds.project_by_name.values():
        samples = project.sample_by_name.values()
        threads = len(samples) if not is_local() else 1
        info('Threads number: ' + str(threads))
        project.concat_fastqs(cnf)

        if cnf.targqc or cnf.metamapping:
            info()
            downsample_to = 1e6
            info('Downsampling the reads to ' + str(downsample_to))
            bam_by_sample = dict()
            fastqs = Parallel(n_jobs=threads) \
                (delayed(downsample_fastq)(CallCnf(cnf.__dict__), sample, reads_num=downsample_to/2) \
                    for sample in samples)
            lefts = [l for l, r in fastqs]
            rights = [r for l, r in fastqs]

            samtools = get_system_path(cnf, 'samtools')
            bwa = get_system_path(cnf, 'bwa')
            seqtk = get_system_path(cnf, 'seqtk')
            bammarkduplicates = get_system_path(cnf, 'bammarkduplicates')
            if samtools and bwa and seqtk and bammarkduplicates:
                info()
                info('Alignming ' + str(downsample_to) + ' random reads to the reference')
                aligned = Parallel(n_jobs=threads)(delayed(align)(CallCnf(cnf.__dict__), s, l, r,
                    samtools,
                    bwa,
                    seqtk,
                    bammarkduplicates,
                    cnf.genome.seq,
                    cnf.is_pcr) for s, l, r in zip(samples, lefts, rights))
                for sample, bam_fpath in zip(samples, aligned):
                    bam_by_sample[sample.name] = bam_fpath
                if any(a is None for a in aligned):
                    pass
                else:
                    if cnf.metamapping:
                        info()
                        info('Metamapping for contamination')
                        safe_mkdir(ds.downsample_metamapping_dirpath)
                        run_metamapping(cnf, samples, bam_by_sample, project.downsample_metamapping_dirpath)

                    if cnf.targqc:
                        info()
                        cnf.work_dir = join(cnf.work_dir, source.targqc_name)
                        safe_mkdir(cnf.work_dir)
                        info('Making TargQC reports for BAMs from ' + str(downsample_to) + ' reads')
                        safe_mkdir(project.downsample_targqc_dirpath)
                        project.downsample_targqc_report_fpath = run_targqc(cnf, project, bam_by_sample)
                        cnf.work_dir = dirname(cnf.work_dir)
            else:
                err('For downsampled targqc and metamappint, bwa, samtools and seqtk are required.')

        if cnf.fastqc:
            if not cnf.expose_to_ngs_server_only:
                info('Making FastQC reports')
                safe_mkdir(project.fastqc_dirpath)
                make_fastqc_reports(cnf, samples, project.fastq_dirpath, project.fastqc_dirpath, project.comb_fastqc_fpath)

        new_project_symlink = join(dirname(project_dirpath), project.az_project_name)
        if not exists(new_project_symlink):
            info()
            info('Creating symlink in Datasets now called as project-name: ' + project_dirpath + ' -> ' + new_project_symlink)
            os.symlink(project_dirpath, new_project_symlink)

        # Creating analysis directory
        __prepare_analysis_directory(cnf.work_dir, project.az_project_name, project_dirpath, samples)

        # Making project-level report
        make_project_level_report(cnf, dataset_structure=ds, dataset_project=project)

        # Exposing
        info()
        info('Synking with the NGS webserver')
        html_report_url = sync_with_ngs_server(cnf,
            jira_url=cnf.jira,
            project_name=project.az_project_name,
            sample_names=[s.name for s in samples],
            dataset_dirpath=project_dirpath,
            summary_report_fpath=project.project_report_html_fpath,
            jira_case=jira_case)

            # FastQC
            # symlink_to_ngs(project_dirpath, NGS_WEBSERVER_PREPROC_DIR)

            # if symlink_to_ngs(comb_fastqc_fpath, ngs_webserver_project_dirpath) is None:
            #     err('Error: cannot connect to the ngs server and make symlinks')
            # else:
            #     # BaseCalls
            #     basecall_stats_dirnames = [fname for fname in os.listdir(basecalls_dirpath) if fname.startswith('Basecall_Stats_')]
            #     if len(basecall_stats_dirnames) > 1:
            #         err('More than 1 Basecall_Stats_* dirs found in unalign_dirpath')
            #     if len(basecall_stats_dirnames) == 0:
            #         err('No Basecall_Stats_* dirs found in unalign_dirpath')
            #     if len(basecall_stats_dirnames) == 1:
            #         basecall_stats_dirpath = join(basecalls_dirpath, basecall_stats_dirnames[0])
            #         fpaths = filter(None, (verify_file(join(basecall_stats_dirpath, html_fname)
            #             for html_fname in ['Demultiplex_Stats.htm', 'All.htm', 'IVC.htm'])))
            #         symlink_to_ngs(fpaths, ngs_webserver_project_dirpath)
            #
            #     # Sample sheet
            #     symlink_to_ngs(sample_sheet_csv_fpath, ngs_webserver_project_dirpath)
            #
            #     # TargQC downsampled
            #     symlink_to_ngs(targqc_html_fpath, ngs_webserver_project_dirpath)

            # jira_case = None
            # if jira_url:
            #     # Add to the NGS list
            #     jira_case = retrieve_jira_info(jira_url)
            #
            # sync_with_ngs_server(cnf, jira_case=jira_case,
            #     project_name=cnf.project_name, sample_names=[s.name for s in samples])

        subj = project.project_report_html_fpath or project.name
        txt = 'Preproc finished for ' + project.name + '\n'
        txt += '\n'
        txt += 'Path: ' + project.dirpath + '\n'
        txt += 'Report: ' + str(html_report_url) + '\n'
        if jira_url:
            txt += 'Jira: ' + jira_url
        send_email(txt, subj)

    info()
    info('*' * 70)
    # if not cnf.debug and cnf.work_dir:
    #     try:
    #         shutil.rmtree(cnf.work_dir)
    #     except OSError:
    #         err('Can\'t remove work directory ' + cnf.work_dir + ', please, remove it manually.')


def downsample_fastq(cnf, sample, reads_num=5e5):
    # downsampled_reads_fpath = join(cnf.work_dir, sample.name + '_' + str(reads_num) + '.fastq')
    info('Downsampling ' + sample.name + ' to ' + str(int(reads_num)))
    l_fpath, r_fpath = downsample(cnf, sample.l_fpath, sample.r_fpath, int(reads_num))
    return l_fpath, r_fpath


def align(cnf, sample, l_fpath, r_fpath, samtools, bwa, seqtk, bammarkduplicates, ref, is_pcr=False):
    sam_fpath = join(cnf.work_dir, sample.name + '_downsampled.sam')
    bam_fpath = splitext(sam_fpath)[0] + '.bam'
    sorted_bam_fpath = add_suffix(bam_fpath, 'sorted')

    bwa_cmdline = '{seqtk} mergepe {l_fpath} {r_fpath} | {bwa} mem {ref} -'.format(**locals())
    res = call(cnf, bwa_cmdline, output_fpath=sam_fpath, exit_on_error=False)
    if not res:
        return None

    cmdline = '{samtools} view -Sb {sam_fpath}'.format(**locals())
    call(cnf, cmdline, output_fpath=bam_fpath)

    prefix = splitext(sorted_bam_fpath)[0]
    cmdline = '{samtools} sort {bam_fpath} {prefix}'.format(**locals())
    call(cnf, cmdline, output_fpath=sorted_bam_fpath, stdout_to_outputfile=False)

    if not is_pcr:
        markdup_bam_fpath = markdup_bam(cnf, sorted_bam_fpath, bammarkduplicates)
        if markdup_bam_fpath:
            sorted_bam_fpath = markdup_bam_fpath

    index_bam(cnf, sorted_bam_fpath, samtools=samtools)
    return sorted_bam_fpath


def run_metamapping(cnf, samples, bam_by_sample, output_dirpath):
    info('Running MetaMapping for downsampled BAMs')


def run_targqc(cnf, project, bam_by_sample):
    info('Running TargQC for downsampled BAMs')

    targqc = get_script_cmdline(cnf, 'python', 'targqc.py', is_critical=True)
    bam_fpaths = ' '.join(bam_by_sample[s.name] + ',' + s.name for s in project.sample_by_name.values())
    targqc_work_dir = join(cnf.work_dir, 'TargQC')
    targqc_log_dir = join(cnf.log_dir, 'TargQC')
    safe_mkdir(targqc_work_dir)
    safe_mkdir(targqc_log_dir)
    cmdl = '{targqc} --sys-cnf {cnf.sys_cnf} {bam_fpaths} --bed {cnf.bed} ' \
           '--work-dir {targqc_work_dir} --log-dir {targqc_log_dir} --project-name {cnf.project_name} ' \
           '-o {project.downsample_targqc_dirpath} --genome {cnf.genome.name}'.format(**locals())
    if cnf.reuse_intermediate:
        cmdl += ' --reuse'
    call(cnf, cmdl)
    info('Waiting for targqc to be done...')
    while True:
        if isfile(project.downsample_targqc_report_fpath):
            break
    verify_file(project.downsample_targqc_report_fpath, is_critical=True)
    return project.downsample_targqc_report_fpath

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


def run_fastqc(cnf, fastq_fpath, output_basename, fastqc_dirpath, need_downsample=True):
    fastqc = get_system_path(cnf, 'fastqc', is_critical=True)
    java = get_system_path(cnf, 'java', is_critical=True)
    tmp_dirpath = join(cnf.work_dir, 'FastQC_' + output_basename + '_tmp')
    safe_mkdir(tmp_dirpath)
    cmdline_l = '{fastqc} --dir {tmp_dirpath} --extract -o {fastqc_dirpath} -f fastq -j {java} {fastq_fpath}'.format(**locals())
    j = submit_job(cnf, cmdline_l, 'FastQC_' + output_basename, run_on_chara=True, stdout_to_outputfile=False,
        output_fpath=join(fastqc_dirpath, output_basename + '_fastqc', 'fastqc_report.html'))
    return j


# def run_fastqc(cnf, sample, fastqc_dirpath, need_downsample=True):
#     # with tx_tmpdir(fastqc_work_dir, fastqc_dirpath) as fastqc_out_tx_dirpath:
#     # cmdline = get_script_cmdline(cnf, 'python', join('scripts', 'pre', 'fastqc.py'))
#     # cmdline += (' --sys-cnf {cnf.sys_cnf} --sample {sample.name} -1 {sample.l_fpath} -2 {sample.r_fpath} -o {fastqc_dirpath}'.format(**locals()))
#     fastqc = get_system_path(cnf, 'fastqc', is_critical=True)
#     java = get_system_path(cnf, 'java', is_critical=True)
#     cmdline_l = '{fastqc} --extract -o {fastqc_dirpath} -f fastq -j {java} {sample.l_fpath}'.format(**locals())
#     cmdline_r = '{fastqc} --extract -o {fastqc_dirpath} -f fastq -j {java} {sample.r_fpath}'.format(**locals())
#     j_l = submit_job(cnf, cmdline_l, 'FastQC_' + sample.l_fastqc_base_name, stdout_to_outputfile=False,
#         output_fpath=join(sample.fastqc_dirpath, sample.l_fastqc_base_name + '_fastqc', 'fastqc_report.html'))
#     j_r = submit_job(cnf, cmdline_r, 'FastQC_' + sample.r_fastqc_base_name, stdout_to_outputfile=False,
#         output_fpath=join(sample.fastqc_dirpath, sample.r_fastqc_base_name + '_fastqc', 'fastqc_report.html'))
#
#     return j_l, j_r

    # parser = FastQCParser(fastqc_out, data["name"][-1])
    # stats = parser.get_fastqc_summary()
    # parser.save_sections_into_file()


def find_fastq_pairs_by_sample_names(fastq_fpaths, sample_names):
    fastq_by_sn = OrderedDict()

    for sn in sample_names:
        sn_fastq_fpaths = sorted([f for f in fastq_fpaths if basename(f).startswith(sn + '_R')])
        if len(sn_fastq_fpaths) == 0:
            err('Error: no fastq found for ' + sn)
            fastq_by_sn[sn] = None
        elif len(sn_fastq_fpaths) > 2:
            critical('Error: more than 2 fastq files starting with ' + sn + '_R: ' + ', '.join(sn_fastq_fpaths))
        elif len(sn_fastq_fpaths) == 1:
            warn('Warning: only single fastq file is found for ' + sn + '. Treating as single reads.')
            fastq_by_sn[sn] = [verify_file(sn_fastq_fpaths[0], description='sn_fastq_fpaths[0] for ' + str(sn)), None]
        else:
            fastq_by_sn[sn] = [verify_file(fpath, description='fpath from sn_fastq_fpaths for ' + str(sn)) for fpath in sn_fastq_fpaths]

    return fastq_by_sn


class FQC_Sample:
    def __init__(self, name, fastq_fpath, sample):
        self.name = name
        self.fastq_fpath = fastq_fpath
        self.sample = sample
        self.fastqc_html_fpath = sample.find_fastqc_html(name)


def make_fastqc_reports(cnf, samples, fastq_dirpath, fastqc_dirpath, comb_fastqc_fpath):
    # if isdir(fastqc_dirpath):
    #     if isdir(fastqc_dirpath + '.bak'):
    #         try:
    #             shutil.rmtree(fastqc_dirpath + '.bak')
    #         except OSError:
    #             pass
    #     if not isdir(fastqc_dirpath + '.bak'):
    #         os.rename(fastqc_dirpath, fastqc_dirpath + '.bak')
    # if isdir(fastqc_dirpath):
    #     err('Could not run and combine fastqc because it already exists and could not be moved to fastqc.bak')
    #     return None

    fastqc = get_system_path(cnf, 'fastqc')
    if not fastqc:
        err('FastQC is not found, cannot make reports')
        return None

    else:
        safe_mkdir(fastqc_dirpath)

        fqc_samples = []
        fastqc_jobs = []
        for s in samples:
            fqc_samples.extend([
                FQC_Sample(name=s.l_fastqc_base_name, fastq_fpath=s.l_fpath, sample=s),
                FQC_Sample(name=s.r_fastqc_base_name, fastq_fpath=s.r_fpath, sample=s)])

        for fqc_s in fqc_samples:
            if cnf.reuse_intermediate and verify_file(fqc_s.fastqc_html_fpath, silent=True):
                info(fqc_s.fastqc_html_fpath + ' exists, reusing')
            else:
                fastqc_jobs.append(run_fastqc(cnf, fqc_s.fastq_fpath, fqc_s.name, fastqc_dirpath))
            info()

        wait_for_jobs(cnf, fastqc_jobs)

        fastqc_jobs = []
        while True:
            for fqc_s in fqc_samples:
                fqc_s.fastqc_html_fpath = fqc_s.sample.find_fastqc_html(fqc_s.name)
            not_done_fqc = [fqc_s for fqc_s in fqc_samples
                if not verify_file(fqc_s.fastqc_html_fpath, description='Not found FastQC html for ' + fqc_s.name)]
            if not not_done_fqc:
                info('')
                info('Every FastQC job is done, moving on.')
                info('-' * 70)
                break
            else:
                info('')
                info('Some FastQC jobs are not done (' + ', '.join(f.name for f in not_done_fqc) + '). Retrying them.')
                info('')
                for fqc_s in not_done_fqc:
                    fastqc_jobs.append(run_fastqc(cnf, fqc_s.fastq_fpath, fqc_s.name, fastqc_dirpath))
                wait_for_jobs(cnf, fastqc_jobs)

        for fqc_s in fqc_samples:
            sample_fastqc_dirpath = join(fastqc_dirpath, fqc_s.name + '_fastqc')
            if isfile(sample_fastqc_dirpath + '.zip'):
                try:
                    os.remove(sample_fastqc_dirpath + '.zip')
                except OSError:
                    pass

        write_fastqc_combo_report(comb_fastqc_fpath, fqc_samples)
        verify_file(comb_fastqc_fpath, is_critical=True)
        info('Combined FastQC saved to ' + comb_fastqc_fpath)
        return comb_fastqc_fpath


def __prepare_analysis_directory(work_dir, project_name, project_dirpath, samples):
    kind = next((kind for pref, kind in source.project_kind_by_prefix.items() if project_name.startswith(pref)), None)
    if kind:
        analysis_proj_dirpath = adjust_path(join(project_dirpath.split('/datasets/')[0], 'analysis', kind, project_name))
        if not exists(analysis_proj_dirpath):
            info('Analysis directory ' + analysis_proj_dirpath + ' does not exist. Creating and preparing...')
            safe_mkdir(analysis_proj_dirpath)

            bcbio_csv_fpath = join(analysis_proj_dirpath, 'bcbio.csv')
            if not isfile(bcbio_csv_fpath):
                with file_transaction(work_dir, bcbio_csv_fpath) as tx:
                    with open(tx, 'w') as f:
                        f.write('samplename,description,batch,phenotype\n')
                        for s in samples:
                            f.write(s.name + ',' + s.name + ',' + s.name + '-batch,tumor\n')

        ds_symlink = join(analysis_proj_dirpath, 'dataset')
        if isdir(analysis_proj_dirpath) and not exists(ds_symlink):
            info('Creating symlink in analysis to datasets: ' + project_dirpath + ' -> ' + ds_symlink)
            os.symlink(project_dirpath, ds_symlink)


if __name__ == '__main__':
    main()


'''
#!/bin/bash/

#Takes 2 arguements, data_loc and project_name, as created in datasets, such as hiseq and Dev_0200_HiSeq_DS
#Usage - upload_seqQC.sh hiseq Dev_0200_HiSeq_DS https://jira.rd.astrazeneca.net/browse/NGSG-313
#Usage - upload_seqQC.sh bioscience Bio_0041_IDT_RR_DS

datasets=/ngs/oncology/datasets
data_loc=$1
project_name=$2

cd /opt/lampp/htdocs/seqQC/
#echo "In /opt/lampp/htdocs/seqQC on NGS Server"
echo " "

mkdir $project_name
cd $project_name
mkdir FastQC

echo "Demultiplex Report linked!"
echo " "

echo "SampleSheet linked!"
echo "DONE!"
echo " "

'''