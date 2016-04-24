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

from ext_modules.joblib import Parallel, delayed
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
from source import logger
from source.logger import info, critical, err, is_local, warn, send_email
from source.utils import is_az
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug, check_system_resources, determine_sys_cnf, determine_run_cnf, \
    check_genome_resources, set_up_log
from source.file_utils import safe_mkdir, verify_dir, verify_file, adjust_path, \
    add_suffix, file_transaction, splitext_plus
from source.webserver.exposing import sync_with_ngs_server

NGS_WEBSERVER_PREPROC_DIR = '/opt/lampp/htdocs/reports'
if is_local():
    NGS_WEBSERVER_PREPROC_DIR = '/Users/vlad/Sites/reports'


def proc_opts():
    parser = OptionParser()
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser)
    parser.add_option('--expose-only', dest='expose_to_ngs_server_only', action='store_true', default=False, help='Only add project to the webserver')
    parser.add_option('--no-expose', dest='expose', action='store_false', default=True, help='Do not expose the reports')
    parser.add_option('-o', dest='output_dir')

    (opts, args) = parser.parse_args()
    logger.is_debug = opts.debug

    if len(args) < 1:
        critical('Usage: ' + __file__ + ' *.fq.gz -o output_dir')
    # if len(args) < 2:
    #     info('No dataset path specified, assuming it is the current working directory')
    #     dataset_dirpath = adjust_path(os.getcwd())
    #     jira_url = args[0]

    fastq_fpaths = [verify_file(fpath) for fpath in args]

    run_cnf = determine_run_cnf(opts)
    cnf = Config(opts.__dict__, determine_sys_cnf(opts), run_cnf)

    cnf.output_dir = adjust_path(cnf.output_dir)
    info('Writing to ' + str(cnf.output_dir))

    cnf.project_name = cnf.project_name or ''

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
            if exists(latest_fpath):
                os.remove(latest_fpath)
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

    # check_genome_resources(cnf)
    check_system_resources(cnf, optional=['fastq'])

    return cnf, cnf.output_dir, fastq_fpaths


def main():
    cnf, output_dir, fastq_fpaths = proc_opts()

    if not cnf.expose_to_ngs_server_only:
        info('Making FastQC reports')
        safe_mkdir(output_dir)
        make_fastqc_reports(cnf, fastq_fpaths, output_dir)


    # Making project-level report
    # make_project_level_report(cnf, dataset_structure=ds, dataset_project=project)

    # # Exposing
    # info()
    # info('Syncing with the NGS webserver')
    # html_report_url = sync_with_ngs_server(cnf,
    #     jira_url=jira_by_subprj.get(project.name, jira_by_subprj.values()[0] if jira_by_subprj.values() else 'unreached'),
    #     project_name=project.az_project_name,
    #     sample_names=[s.name for s in samples],
    #     dataset_dirpath=project_dirpath,
    #     summary_report_fpath=project.project_report_html_fpath,
    #     jira_case=jira_case_by_subprj.get(project.name, jira_case_by_subprj.values()[0] if jira_case_by_subprj.values() else None)
    # )

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

        # subj = project.project_report_html_fpath or project.name
        # txt = 'Preproc finished for ' + project.name + '\n'
        # txt += '\n'
        # txt += 'Path: ' + project.dirpath + '\n'
        # txt += 'Report: ' + str(html_report_url) + '\n'
        # if jira_url:
        #     txt += 'Jira: ' + jira_url
        # send_email(txt, subj)

    info()
    info('*' * 70)
    # if not cnf.debug and cnf.work_dir:
    #     try:
    #         shutil.rmtree(cnf.work_dir)
    #     except OSError:
    #         err('Can\'t remove work directory ' + cnf.work_dir + ', please, remove it manually.')


def run_fastqc(cnf, fastq_fpath, output_basename, fastqc_dirpath, need_downsample=True):
    fastqc = get_system_path(cnf, 'fastqc', is_critical=True)
    java = get_system_path(cnf, 'java', is_critical=True)
    tmp_dirpath = join(cnf.work_dir, 'FastQC_' + output_basename + '_tmp')
    safe_mkdir(tmp_dirpath)
    cmdline_l = '{fastqc} --dir {tmp_dirpath} --extract -o {fastqc_dirpath} -f fastq -j {java} {fastq_fpath}'.format(**locals())
    j = submit_job(cnf, cmdline_l, 'FastQC_' + output_basename, run_on_chara=True, stdout_to_outputfile=False)
        # output_fpath=join(fastqc_dirpath, output_basename + '_fastqc', 'fastqc_report.html'))
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
    def __init__(self, name, fastq_fpath, sample=None):
        self.name = name
        self.fastq_fpath = fastq_fpath
        self.sample = sample


def make_fastqc_reports(cnf, fastq_fpaths, output_dir):
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
        safe_mkdir(output_dir)

        fqc_samples = []
        fastqc_jobs = []
        for fastq_fpath in fastq_fpaths:
            s = FQC_Sample(name=splitext_plus(basename(fastq_fpath))[0], fastq_fpath=fastq_fpath)
            fqc_samples.extend([s])
            info('Added sample ' + s.name)

        for fqc_s in fqc_samples:
            if cnf.reuse_intermediate and verify_file(fqc_s.fastqc_html_fpath, silent=True):
                info(fqc_s.fastqc_html_fpath + ' exists, reusing')
            else:
                fastqc_jobs.append(run_fastqc(cnf, fqc_s.fastq_fpath, fqc_s.name, output_dir))
            info()

        wait_for_jobs(cnf, fastqc_jobs)

        fastqc_jobs = []
        # while True:
        for fqc_s in fqc_samples:
            fqc_s.fastqc_html_fpath = find_fastqc_html(output_dir, fqc_s.name)
        not_done_fqc = [fqc_s for fqc_s in fqc_samples
            if not verify_file(fqc_s.fastqc_html_fpath, description='Not found FastQC html for ' + fqc_s.name)]
        # if not not_done_fqc:
        #     info('')
        #     info('Every FastQC job is done, moving on.')
        #     info('-' * 70)
        #     break
        # else:
        #     info('')
        #     info('Some FastQC jobs are not done (' + ', '.join(f.name for f in not_done_fqc) + '). Retrying them.')
        #     info('')
        #     for fqc_s in not_done_fqc:
        #         fastqc_jobs.append(run_fastqc(cnf, fqc_s.fastq_fpath, fqc_s.name, output_dir))
        #     wait_for_jobs(cnf, fastqc_jobs)

        for fqc_s in fqc_samples:
            sample_fastqc_dirpath = join(output_dir, fqc_s.name + '_fastqc')
            if isfile(sample_fastqc_dirpath + '.zip'):
                try:
                    os.remove(sample_fastqc_dirpath + '.zip')
                except OSError:
                    pass

        comb_fastqc_fpath = join(output_dir, 'fastqc.html')
        write_fastqc_combo_report(cnf, comb_fastqc_fpath, fqc_samples)
        verify_file(comb_fastqc_fpath, is_critical=True)
        info('Combined FastQC saved to ' + comb_fastqc_fpath)
        return comb_fastqc_fpath


def find_fastqc_html(fastqc_dirpath, end_name):
    sample_fastqc_dirpath = join(fastqc_dirpath, end_name + '_fastqc')
    if not isdir(sample_fastqc_dirpath):
        sample_fastqc_dirpath = join(fastqc_dirpath, end_name + '.fq_fastqc')
    if not isdir(sample_fastqc_dirpath):
        sample_fastqc_dirpath = join(fastqc_dirpath, end_name + '.fastq_fastqc')
    if not isdir(sample_fastqc_dirpath):
        return None

    fastqc_html_fpath = join(fastqc_dirpath, end_name + '_fastqc.html')
    if isfile(fastqc_html_fpath):
        return fastqc_html_fpath
    else:
        fastqc_html_fpath = join(sample_fastqc_dirpath, 'fastqc_report.html')
        if isfile(fastqc_html_fpath):
            return fastqc_html_fpath
        else:
            return None


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
