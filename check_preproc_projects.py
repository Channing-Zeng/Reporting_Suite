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
import os.path, time

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


def proc_opts():
    parser = OptionParser(description='')
    (opts, args) = parser.parse_args()
    if len(args) < 1:
        critical('First argument should be a root datasets dir')
    # if len(args) < 2:
    #     info('No dataset path specified, assuming it is the current working directory')
    #     dataset_dirpath = adjust_path(os.getcwd())
    #     jira_url = args[0]
    root_dirpath = verify_dir(args[0], is_critical=True, description='Dataset directory')  # /ngs/oncology/datasets/hiseq/150521_D00443_0159_AHK2KTADXX

    info(' '.join(sys.argv))

    return root_dirpath


def main():
    root_dirpath = proc_opts()
    info('*' * 60)
    info()

    all_issues = []

    info('Iterating over ' + root_dirpath)
    info('-' * 60)
    info()
    for fname in os.listdir(root_dirpath):
        if fname.startswith('.'):
            continue
        info(fname)
        project_dirpath = join(root_dirpath, fname)
        if isdir(project_dirpath) \
                and isfile(join(project_dirpath, 'SampleSheet.csv')) \
                and isdir(join(project_dirpath, 'Unalign')):
            info('Unalign and SampleSheet.csv found')

            ds = DatasetStructure.create(project_dirpath, '')

            issues = []
            if not ds.project_by_name:
                err('No projects found')
            else:
                info('Projects: ' + ', '.join([p.name + ' (' + ', '.join(p.sample_by_name) + ')' for p in ds.project_by_name.values()]))
                for project in ds.project_by_name.values():
                    if not project.sample_by_name:
                        err('No samples for project ' + project.name + ' found')
                    else:
                        for i, s1 in enumerate(project.sample_by_name.values()):
                            for s2 in project.sample_by_name.values()[i + 1:]:
                                if s2.name.startswith(s1.name):
                                    issues.append('   issued samples: ' + s1.name + ' and ' + s2.name + ' from ' + project.name)

            if issues:
                all_issues.append(fname + ' created: %s, last modified: %s' %
                                  (time.ctime(os.path.getctime(project_dirpath)),
                                   time.ctime(os.path.getmtime(project_dirpath))))
                all_issues.extend(issues)
                all_issues.append('')

            info()
            info('-' * 60)

    info()
    info('Failed projects: ')
    for msg in all_issues:
        info(msg)


if __name__ == '__main__':
    main()