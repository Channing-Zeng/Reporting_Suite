#!/usr/bin/env python

import sys
from os.path import join, splitext, basename, realpath, isfile, getsize, dirname, exists


class TargetSeqAnalyzer:
    def __init__(self, bam, bed, ref, pad=500):
        self.bam = bam
        self.bed = bed
        self.ref = ref
        self.pad = pad

        self.reads = None
        self.target_reads = None
        self.percent_target_reads = None
        self.percent_padded_target_reads = None

        self.read_bases = None
        self.target_read_bases = None
        self.percent_target_read_bases = None

        self.reference_target_bases = None
        self.covered_bases




def main(args):
    bam = args[0]
    bed = args[1]
    ref = args[2]
    pad = args[3] if len(args) > 2 else 500

    ts = TargetSeqAnalyzer(bam, bed, ref, pad)

    with open('report.txt') as rep:
        rep.write('Number of mapped reads\t%d\n' % ts._reads())
        rep.write('Number of reads on target\t%d\n' % ts._target_reads())
        rep.write('Percent reads on target\t%d\n' % ts._())
        rep.write('Percent reads on padded target\t%d\n' % ts._number_mapped_reads())
        rep.write('Total aligned base reads\t%d\n' % ts._number_mapped_reads())
        rep.write('Total base reads on target\t%d\n' % ts._number_mapped_reads())
        rep.write('Percent base reads on target\t%d\n' % ts._number_mapped_reads())
        rep.write('Bases in targeted reference\t%d\n' % ts._number_mapped_reads())
        rep.write('Bases covered (at least 1x)\t%d\n' % ts._number_mapped_reads())
        rep.write('Average base coverage depth\t%d\n' % ts._number_mapped_reads())
        rep.write('Maximum base read depth\t%d\n' % ts._number_mapped_reads())
        rep.write('Average base read depth\t%d\n' % ts._number_mapped_reads())
        rep.write('Std.Dev base read depth\t%d\n' % ts._number_mapped_reads())
        for x in [1, 5, 10, 20, 50, 100, 500, 1000, 2500, 5000, 10000]:
            rep.write('Target coverage at %dx\t%d\n' % (x, ts._target_cov_at(x)))



if __name__ == '__main__':
    main(sys.argv[1:])