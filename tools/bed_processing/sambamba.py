#!/usr/bin/env python
import bcbio_postproc  # do not remove it: checking for python version and adding site dirs inside

import sys
import os
import subprocess
from os.path import basename, join, isfile, getctime
from source.file_utils import splitext_plus, verify_file
from source.logger import err, info
from source.targetcov.bam_and_bed_utils import bam_to_bed, bedtools_version, bam_to_bed_nocnf


def main():
    args = sys.argv[1:]

    if len(args) < 2:
        sys.exit('Usage: ' + __file__ + ' bam commandline sambamba')

    bam, cmdline, sambamba = args

    index_bam(bam, sambamba)

    subprocess.call((sambamba + ' ' + cmdline).split())


def index_bam(bam_fpath, sambamba):
    indexed_bam = bam_fpath + '.bai'
    if not isfile(indexed_bam) or getctime(indexed_bam) < getctime(bam_fpath):
        info('Indexing BAM, writing ' + indexed_bam + '...')
        cmdline = '{sambamba} index {bam_fpath}'.format(**locals())
        subprocess.call(cmdline.split())


if __name__ == '__main__':
    main()
