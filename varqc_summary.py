#!/usr/bin/env python

from __future__ import print_function
import sys
from os.path import join
from source.utils_from_bcbio import file_exists
from source.logger import info
from source.utils import critical
from source.summarize import summarize_qc

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))


def main(args):  # dir samples.txt report_basedir
    if len(args) < 2:
        critical('Usage: ' + __file__ + ' <bcbio_final_dir> '
                 '<list_of_sample_names> [report_base_dir (default varqc)] '
                 '[vcf_suffic (defauls is -mutect)]')

    out_dirpath = args[0]
    samples_fpath = args[1]
    report_basedir = args[2] if len(args) >= 3 else 'varqc'
    vcfsuffix = args[3] if len(args) >= 4 else '-mutect'

    report_suffix = '_qc.report'
    summary_report_fpath = join(out_dirpath, 'varqc_summary_report.txt')
    report_fpaths = []

    with open(samples_fpath, 'r') as f:
        for line in f:
            sample_name = line.strip()

            report_fpath = join(out_dirpath, sample_name, report_basedir, sample_name + report_suffix)
            info(report_fpath)

            if file_exists(report_fpath):
                report_fpaths.append(report_fpath)
            else:
                info(report_fpath + ' does not exist, checking another')
                report_fpath = join(out_dirpath, sample_name, report_basedir, sample_name + vcfsuffix + report_suffix)
                info(report_fpath)

                if file_exists(report_fpath):
                    report_fpaths.append(report_fpath)
                else:
                    info(report_fpath + ' does not exist, skipping')
            info('')

    summarize_qc(report_fpaths, summary_report_fpath, report_suffix)


if __name__ == '__main__':
    main(sys.argv[1:])
