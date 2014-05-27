#!/usr/bin/env python

from __future__ import print_function
import sys
import os
from optparse import OptionParser
from os.path import join

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))


def main(argv):
    pass
    # parser = OptionParser()
    # parser.add_option('-d', '--dir', dest='dir', metavar='DIR')
    # parser.add_option('-s', '--samples', dest='samples')
    # parser.add_option('--vcf-basename', dest='vcf_basename', help='-mutect')
    # parser.add_option('-b', '--bed', dest='bed')
    # parser.add_option('-t', '--nt', dest='threads', default=4)
    # parser.add_option('-w', dest='rewrite', action='store_true', default=False)
    # parser.add_option('--np', '--no-parallel', dest='no_parallel', action='store_true', default=False)
    # (options, args) = parser.parse_args()
    #
    # with open(options.samples) as samples:
    #     for sample in samples:
    #         sample_dir = join(dir, sample)


if __name__ == '__main__':
    main(sys.argv)

