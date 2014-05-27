#!/usr/bin/env python

from __future__ import print_function
import sys
import os
from optparse import OptionParser
from os.path import join
from genericpath import isfile
from source.utils import critical
from source.summarize import summarize_cov

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))


def main(argv):  # dir samples.txt report_basedir
    if len(argv) != 4:
        critical('Error: provide bcbio dir, list of samples, and report basedir (e.g. VarQC)')

    out_dirpath = argv[1]
    samples_fname = argv[2]
    report_basedir = argv[3]
    report_suffix='.targetseq.summary.txt'
    summary_report_fpath=join(out_dirpath, 'summary' + report_suffix)
    report_fpaths = []
    with open(samples_fname, 'r') as f:
        for line in f:
            sample_name = line.strip()
            report_fpath = join(out_dirpath, sample_name, report_basedir, sample_name + report_suffix)
            if isfile(report_fpath):
                report_fpaths.append(report_fpath)
    summarize_cov(report_fpaths, summary_report_fpath, report_suffix)


if __name__ == '__main__':
    main(sys.argv)
