#!/usr/bin/env python
import bcbio_postproc
import subprocess
import sys
from optparse import OptionParser
from source.file_utils import which, add_suffix
from source.logger import info


def main():
    info(' '.join(sys.argv))
    info()
    parser = OptionParser()
    parser.add_option('--bam', dest='bam', help='BAM file.')
    parser.add_option('--bed', dest='bed', help='BED file.')
    parser.add_option('-g', '--genome', dest='genome', help='Genome file.')
    parser.add_option('-o', dest='output_fpath', help='Output file.')
    opts, args = parser.parse_args()

    intermediate_fpath = opts.output_fpath.replace('.txt', '.bam')
    bedtools = which('bedtools')
    sambamba = which('sambamba')

    cmdline = '{sambamba} view -L {opts.bed} {opts.bam} -f bam > {intermediate_fpath}'.format(**locals())
    subprocess.call(cmdline, shell=True)
    cmdline = '{bedtools} genomecov -bg -g {opts.genome} -ibam {intermediate_fpath} > {opts.output_fpath}'.format(**locals())
    subprocess.call(cmdline, shell=True)


if __name__ == '__main__':
    main()
