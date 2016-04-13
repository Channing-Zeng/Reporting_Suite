#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc

import shutil
from source.clinical_reporting.clinical_parser import ClinicalExperimentInfo, clinical_sample_info_from_cnf

from source.clinical_reporting.clinical_reporting import make_clinical_report
from source.logger import info
from source.prepare_args_and_cnf import check_system_resources, check_genome_resources
from source.main import read_opts_and_cnfs


def main():
    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--targqc-dir'], dict(
                dest='targqc_dirpath',
            )),
            (['--mutations'], dict(
                dest='mutations_fpath',
            )),
            (['--sv'], dict(
                dest='sv_fpath',
            )),
            (['--sv-vcf'], dict(
                dest='sv_vcf_fpath',
            )),
            (['--varqc'], dict(
                dest='varqc_json_fpath',
            )),
            (['--varqc-after'], dict(
                dest='varqc_after_json_fpath',
            )),
            (['--target-type'], dict(
                dest='target_type',
                default='panel',
            )),
            (['--bed'], dict(
                dest='bed_fpath',
            )),
            (['--seq2c'], dict(
                dest='seq2c_tsv_fpath',
            )),
            (['--project-level-report'], dict(
                dest='project_report_path',
            )),
            (['--targqc-html'], dict(
                dest='targqc_report_path',
            )),
            (['--match'], dict(
                dest='match_sample_name',
            )),
            (['--jira'], dict(
                dest='jira_url',
            )),

        ],
        key_for_sample_name=None,
        required_keys=[],
        file_keys=['mutations_fpath',
                   'varqc_json_fpath',
                   'varqc_after_json_fpath',
                   'bed_fpath',
                   'seq2c_tsv_fpath',
                   'sv_fpath',
                   'sv_vcf_fpath',
                   #'project_report_path',  # DO NOT UNCOMMENT! Project level report might not yet exist
                   ],
            # do not check mutations_fpath! could be either of:
            #   vardict.PASS.txt,
            #   vardict-java.PASS.txt,
            #   vardict.single.PASS.txt,
            #   vardict.paired.PASS.txt,
        dir_keys=['targqc_dirpath'],
    )

    check_genome_resources(cnf)
    check_system_resources(
        cnf,
        required=['samtools', 'bedtools'],
        optional=[])

    clin_info = clinical_sample_info_from_cnf(cnf)
    html_fpath = make_clinical_report(cnf, clin_info, clin_info.sample.clinical_html)
    info('Clinical report: ' + html_fpath)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


if __name__ == '__main__':
    main()
