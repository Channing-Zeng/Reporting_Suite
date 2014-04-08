'''
samtools view -bh your.sorted.bam chr1:100,000-20,000 (region of interest) > regionofinterest.bam
'''
import sys
from os import system
import subprocess


def call(cmd, out):
    print cmd
    res = subprocess.call(cmd, shell=True, stdout=open(out, 'w'))


def main(args):
    bam = args[0]
    bed = args[1]
    call('samtools depth {bam} -b {bed}'.format(**locals()),
         'covs.txt')

    call('samtools bedcov {bed} {bam}'.format(**locals()),
         'bedcovs.txt')

    call('coverageBed -abam {bam} -b {bed}'.format(**locals()),
         'coverageBed.txt')


# samtools view -u -q 30 -f 0x2 aln.bam | \
#    coverageBed -abam stdin -b exons.bed \
#    > exons.bed.proper.q30.coverage


if __name__ == '__main__':
    main(sys.argv[1:])