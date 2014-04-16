#!/usr/bin/env python

import sys
from os.path import join, splitext, basename, realpath, isfile, getsize, dirname, exists
from subprocess import check_output, call


def analyze(self, bam, bed, ref, pad=500):
    reads = int(check_output('samtools view -c -F 4 {bam}'.format(**locals()), shell=True))
    print 'Number of mapped reads\t%d' % reads

    target_reads = int(check_output('samtools view -c -F 4 {bam} -L {bed}'.format(**locals()), shell=True))
    print 'Number of reads on target\t%d' % target_reads

    percent_target_reads = float(target_reads) / reads
    print 'Percent reads on target\t%f' % percent_target_reads
    # self.percent_padded_target_reads = None

    read_bases = int(check_output("samtools depth {bam} | awk '{c+=$3} END {print c}'".format(**locals()), shell=True))
    print 'Percent reads on padded target\t%f' % read_bases

    target_read_bases = int(check_output("samtools depth {bam} -b {bed} | awk '{c+=$3} END {print c}'".format(**locals()), shell=True))
    print 'Total aligned base reads\t%d' % target_read_bases

    percent_target_read_bases = float(target_read_bases) / read_bases
    print 'Percent base reads on target\t%f' % percent_target_read_bases

    genome_size = int(check_output("awk 'BEGIN {gs = 0} {gs \+= $3-$2} END {print gs}' {bed}".format(**locals()), shell=True))
    print 'Genome size\t%d' % genome_size

    covered_bases = int(check_output("samtools depth {bam} | awk '{c+=$3} END {print c}'".format(**locals())))

    pileup = 'starts.pileup'
    call("samtools depth -b {bed} {bam} 2> /dev/null > {pileup}")

    cov_report = check_output("awk src/coverage_analysis.awk -v basereads={read_bases} "
                              "-v genome={genome_size} {pileup}")
    print cov_report


def main(args):
    bam = args[0]
    bed = args[1]
    ref = args[2]
    pad = args[3] if len(args) > 2 else 500

    analyze(bam, bed, ref, pad)

    # with open('report.txt') as rep:
    #     rep.write('Number of mapped reads\t%d\n' % ts._reads())
    #     rep.write('Number of reads on target\t%d\n' % ts._target_reads())
    #     rep.write('Percent reads on target\t%d\n' % ts._())
    #     rep.write('Percent reads on padded target\t%d\n' % ts._number_mapped_reads())
    #     rep.write('Total aligned base reads\t%d\n' % ts._number_mapped_reads())
    #     rep.write('Total base reads on target\t%d\n' % ts._number_mapped_reads())
    #     rep.write('Percent base reads on target\t%d\n' % ts._number_mapped_reads())
    #     rep.write('Bases in targeted reference\t%d\n' % ts._number_mapped_reads())
    #     rep.write('Bases covered (at least 1x)\t%d\n' % ts._number_mapped_reads())
    #     rep.write('Average base coverage depth\t%d\n' % ts._number_mapped_reads())
    #     rep.write('Maximum base read depth\t%d\n' % ts._number_mapped_reads())
    #     rep.write('Average base read depth\t%d\n' % ts._number_mapped_reads())
    #     rep.write('Std.Dev base read depth\t%d\n' % ts._number_mapped_reads())
    #     for x in [1, 5, 10, 20, 50, 100, 500, 1000, 2500, 5000, 10000]:
    #         rep.write('Target coverage at %dx\t%d\n' % (x, ts._target_cov_at(x)))



if __name__ == '__main__':
    main(sys.argv[1:])