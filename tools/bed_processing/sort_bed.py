#!/usr/bin/env python
from os.path import basename

import bcbio_postproc  # do not remove it: checking for python version and adding site dirs inside

import sys
from optparse import OptionParser
from source.config import Config
from source.file_utils import adjust_path
from source.logger import critical
from source.prepare_args_and_cnf import determine_sys_cnf, check_genome_resources
from source.targetcov.bam_and_bed_utils import verify_bed, sort_bed


def main():
    parser = OptionParser(usage='Usage: ' + basename(__file__) + ' -o Output_BED_file -g hg19 Input_BED_file')
    parser.add_option('-o', '--output-bed', dest='output_fpath')
    parser.add_option('-g', '--genome', dest='genome')
    (opts, args) = parser.parse_args(sys.argv[1:])

    if len(args) < 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)
    cnf = Config(opts.__dict__, determine_sys_cnf(opts), {})

    check_genome_resources(cnf)

    if not cnf.output_fpath:
        critical(parser.usage)

    sort_bed(cnf, verify_bed(args[0], is_critical=True), adjust_path(cnf.output_fpath))


if __name__ == '__main__':
    main()
