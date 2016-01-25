#!/usr/bin/env python
import os
from os.path import basename, join

import bcbio_postproc  # do not remove it: checking for python version and adding site dirs inside

import sys
from logging import info

from source.file_utils import splitext_plus, verify_file
from source.logger import err
from source.targetcov.bam_and_bed_utils import bam_to_bed, bedtools_version, bam_to_bed_nocnf


def main():
    args = sys.argv[1:]

    if len(args) < 2:
        sys.exit('Usage: ' + __file__ + ' work_dir bed bam chr_lengths_fpath bedcov_output_fpath bedtools')

    work_dir, bed, bam, chr_lengths_fpath, bedcov_output_fpath, bedtools = args

    launch_bedcoverage_hist(work_dir, bed, bam, chr_lengths_fpath, bedcov_output_fpath, bedtools)


def launch_bedcoverage_hist(work_dir, bed, bam, chr_lengths_fpath, bedcov_output_fpath=None, bedtools='bedtools'):
    if not bedcov_output_fpath:
        bedcov_output_fpath = join(work_dir,
            splitext_plus(basename(bed))[0] + '__' +
            splitext_plus(basename(bam))[0] + '_bedcov_output.txt')

    if bam.endswith('bam'):
        bam = bam_to_bed_nocnf(bam, bedtools)
    verify_file(bam, is_critical=True, description='bam to bed conversion result')

    v = bedtools_version(bedtools)
    if v and v >= 24:
        cmdline = '{bedtools} coverage -sorted -g {chr_lengths_fpath} -a {bed} -b {bam} -hist'.format(**locals())
    else:
        cmdline = '{bedtools} coverage -a {bam} -b {bed} -hist'.format(**locals())
    cmdline += ' > ' + bedcov_output_fpath
    info(cmdline)
    os.system(cmdline)
    res = verify_file(bedcov_output_fpath)
    if res:
        info('Done, saved to ' + bedcov_output_fpath)
    else:
        err('Error, result is non-existent or empty')


if __name__ == '__main__':
    main()