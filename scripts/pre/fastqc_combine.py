#!/usr/bin/env python
import bcbio_postproc

from os.path import join, basename, dirname, splitext
import sys
from source import info, verify_file
from optparse import OptionParser
from source.fastqc.summarize_fastqc import write_fastqc_combo_report
from source.file_utils import adjust_path
from source.logger import critical


class Sample:
    def __init__(self, name, l_fpath=None, r_fpath=None, fastqc_html_fpath=None):
        self.name = name
        self.l_fpath = l_fpath
        self.r_fpath = r_fpath
        self.fastqc_html_fpath = fastqc_html_fpath


def main():
    info(' '.join(sys.argv))
    info()

    description = 'This script runs preprocessing.'

    parser = OptionParser(description=description)
    parser.add_option('-o', dest='output_fpath', help='Output file path')
    (opts, args) = parser.parse_args()
    if len(args) == 0:
        critical('Please, provide paths to fastqc HTML reports in arguments list.')
    if not opts.output_fpath:
        critical('Please, specify output html path with "-o result.html"')

    input_htmls = [verify_file(arg) for arg in args]
    if len(set(input_htmls)) < len(input_htmls):
        info('Some file names are equal; taking base dirnames as sample names')
        samples = [Sample(name=basename(dirname(fp)), fastqc_html_fpath=fp) for fp in input_htmls]
    else:
        info('All file names are different, taking file basenames as sample names')
        samples = [Sample(name=basename(fp).split('.fq_fastqc.html')[0].split('.htm')[0], fastqc_html_fpath=fp) for fp in input_htmls]

    res_fpath = adjust_path(opts.output_fpath)
    write_fastqc_combo_report(res_fpath, samples)

    verify_file(res_fpath, is_critical=True)

    info()
    info('*' * 70)
    info('Fastqc summary:')
    info('  ' + res_fpath)


if __name__ == '__main__':
    main()
