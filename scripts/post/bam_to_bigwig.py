#!/usr/bin/env python
import bcbio_postproc  # do not remove it: checking for python version and adding site dirs inside

'''Convert BAM files to BigWig file format in a specified region.
Usage:
    bam_to_wiggle.py <BAM file> [<YAML config>]
    [--outfile=<output file name>
     --chrom=<chrom>
     --start=<start>
     --end=<end>
     --normalize]
chrom start and end are optional, in which case they default to everything.
The normalize flag adjusts counts to reads per million.
The config file is in YAML format and specifies the location of the wigToBigWig
program from UCSC:
program:
  ucsc_bigwig: wigToBigWig
If not specified, these will be assumed to be present in the system path.
The script requires:
    pysam (http://code.google.com/p/pysam/)
    wigToBigWig from UCSC (http://hgdownload.cse.ucsc.edu/admin/exe/)
If a configuration file is used, then PyYAML is also required (http://pyyaml.org/)
'''

import os
import sys
from os.path import isfile
from os.path import splitext, join


import source
from source.file_utils import file_transaction
from source.targetcov.bam_and_bed_utils import check_md5, get_bedgraph_coverage
from source.utils import get_ext_tools_dirname, get_chr_len_fpath
from source.logger import critical, info, err
from source.main import read_opts_and_cnfs
from source.prepare_args_and_cnf import check_genome_resources
from source.tools_from_cnf import get_system_path
from source.calling_process import call
from tools.add_jbrowse_tracks import create_jbrowse_symlink


def proc_args(argv):
    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--bam'], dict(
                dest='bam',
            )),
        ],
        required_keys=['bam'],
        file_keys=['bam'],
    )

    check_genome_resources(cnf)

    if not cnf.bam:
        critical('No bam file provided to input')
    if not cnf.genome:
        critical('Please, specify the --genome option (e.g. --genome hg19)')

    return cnf


def process_bam(cnf, bam_fpath, chrom='all', start=0, end=None,
         outfile=None, normalize=False, use_tempfile=False):
    if outfile is None:
        outfile = '%s.bigwig' % splitext(bam_fpath)[0]
    if start > 0:
        start = int(start) - 1
    if end is not None:
        end = int(end)

    bigwig_fpath = outfile
    if not (cnf.reuse_intermediate and os.path.exists(bigwig_fpath) and check_md5(cnf.work_dir, bam_fpath, 'bam', silent=True)):
        bedgraph_fpath = get_bedgraph_coverage(cnf, bam_fpath)
        chr_len_fpath = get_chr_len_fpath(cnf)
        convert_to_bigwig(bedgraph_fpath, cnf, chr_len_fpath, bigwig_fpath)
    else:
        info(outfile + ' exists, reusing')
    return bigwig_fpath


def convert_to_bigwig(bedgraph_fpath, cnf, chr_len_fpath, bw_fpath):
    try:
        with file_transaction(cnf.work_dir, bw_fpath) as tx_fpath:
            cmdl = get_system_path(cnf, join(get_ext_tools_dirname(), 'bedGraphToBigWig'), is_critical=True)
            cmdl += ' ' + bedgraph_fpath + ' ' + chr_len_fpath + ' ' + tx_fpath
            call(cnf, cmdl, exit_on_error=True)
    finally:
        os.remove(bedgraph_fpath)
    return bw_fpath


def main():
    cnf = proc_args(sys.argv)
    bigwig_fpath = process_bam(cnf, cnf.bam)
    if isfile(bigwig_fpath) and cnf.project_name and cnf.sample:
        create_jbrowse_symlink(cnf.genome.name, cnf.project_name, cnf.sample, bigwig_fpath)
        info('BAM was successfully converted.')
    elif not isfile(bigwig_fpath):
        err('BAM was not converted to BigWig.')


if __name__ == '__main__':
    main()



