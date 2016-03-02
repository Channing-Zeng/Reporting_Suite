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
from collections import defaultdict


import source
from source.file_utils import file_transaction, add_suffix, adjust_path, intermediate_fname
from source.targetcov.bam_and_bed_utils import check_md5, bam_to_bed, remove_dups, index_bam, verify_bam
from source.utils import get_ext_tools_dirpath, get_chr_len_fpath, get_chr_lengths
from source.logger import critical, info
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
        chr_len_fpath = get_chr_len_fpath(cnf)
        dedup_bam = intermediate_fname(cnf, bam_fpath, source.dedup_bam)
        if not verify_bam(dedup_bam, silent=True):
            info('Deduplicating bam file ' + bam_fpath)
            remove_dups(cnf, bam_fpath, dedup_bam)
        else:
            info(dedup_bam + ' exists')
        index_bam(cnf, dedup_bam)
        bam_bed_fpath = bam_to_bed(cnf, dedup_bam, to_gzip=False)
        sorted_bed_fpath = sort_bed_by_alphabet(cnf, bam_bed_fpath)
        bedgraph_fpath = '%s.bedgraph' % splitext(bam_fpath)[0]
        with file_transaction(cnf.work_dir, bedgraph_fpath) as tx_fpath:
            bedtools = get_system_path(cnf, 'bedtools')
            cmdl = '{bedtools} genomecov -bg -split -g {chr_len_fpath} -i {sorted_bed_fpath}'.format(**locals())
            call(cnf, cmdl, exit_on_error=True, output_fpath=tx_fpath)

        convert_to_bigwig(bedgraph_fpath, cnf, chr_len_fpath, bigwig_fpath)
    else:
        info(outfile + ' exists, reusing')
    return bigwig_fpath


def sort_bed_by_alphabet(cnf, input_bed_fpath, output_bed_fpath=None):
    chr_lengths = get_chr_lengths(cnf)
    chromosomes = set([c for (c, l) in chr_lengths])
    output_bed_fpath = adjust_path(output_bed_fpath) if output_bed_fpath else add_suffix(input_bed_fpath, 'sorted')

    regions = defaultdict(list)

    info('Sorting regions...')
    chunk_size = 10
    chunk_counter = 0
    with open(input_bed_fpath) as f:
        with file_transaction(cnf.work_dir, output_bed_fpath) as tx:
            with open(tx, 'w') as out:
                for l in f:
                    if not l.strip():
                        continue
                    if l.strip().startswith('#'):
                        out.write(l)
                        continue

                    fs = l.strip().split('\t')
                    chrom = fs[0]
                    if chrom not in chromosomes:
                        continue
                    if chunk_counter == chunk_size or not regions[chrom]:
                        chunk_counter = 0
                        regions[chrom].append('')
                    regions[chrom][-1] += l
                    chunk_counter += 1
                for chr in sorted(regions.keys()):
                    for region in regions[chr]:
                        out.write(region)

    return output_bed_fpath


def convert_to_bigwig(bedgraph_fpath, cnf, chr_len_fpath, bw_fpath):
    try:
        with file_transaction(cnf.work_dir, bw_fpath) as tx_fpath:
            cmdl = get_system_path(cnf, join(get_ext_tools_dirpath(), 'bedGraphToBigWig'), is_critical=True)
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


if __name__ == '__main__':
    main()



