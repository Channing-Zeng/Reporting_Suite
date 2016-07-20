#!/usr/bin/env python
import bcbio_postproc
import sys
from os.path import basename
from optparse import OptionParser

from source.config import Config
from source.logger import info, critical
from source.prepare_args_and_cnf import determine_sys_cnf, check_genome_resources
from source.targetcov.bam_and_bed_utils import get_bedgraph_coverage


def main():
    info(' '.join(sys.argv))
    info()
    parser = OptionParser(usage='Usage: ' + basename(__file__) + ' --bed BED_file --bam BAM_file -g hg19 -o Output_BEDGRAPH_file '
                                                                 '--work-dir work_directory')
    parser.add_option('--bam', dest='bam', help='BAM file.')
    parser.add_option('--bed', dest='bed', help='BED file.')
    parser.add_option('-g', '--genome', dest='chr_len_fpath', help='File with chromosomes lengths.')
    parser.add_option('-o', dest='output_fpath', help='Output file.')
    parser.add_option('--work-dir', dest='work_dir', help='Work directory.')
    (opts, args) = parser.parse_args(sys.argv[1:])

    cnf = Config(opts.__dict__, determine_sys_cnf(opts), {})

    if not cnf.output_fpath:
        critical(parser.usage)
    get_bedgraph_coverage(cnf, cnf.bam, chr_len_fpath=opts.chr_len_fpath, output_fpath=cnf.output_fpath, bed_fpath=cnf.bed)
    info('Done.')


if __name__ == '__main__':
    main()
