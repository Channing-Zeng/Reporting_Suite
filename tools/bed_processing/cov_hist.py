#!/usr/bin/env python

import __check_python_version  # do not remove it: checking for python version and adding site dirs inside

import sys
from os.path import abspath, dirname, realpath, join, basename, splitext
from source.targetcov.coverage_hist import bedcoverage_hist_stats


def main():
    args = sys.argv[1:]

    if len(args) < 2:
        sys.exit('Usage: ' + __file__ + ' bam bed')

    bam, bed = args
    bedcoverage_hist_stats(cnf, splitext(basename(bam))[0], bam, bed)

    # or just script for summarize_bedcoverage_hist_stats?


if __name__ == '__main__':
    main()