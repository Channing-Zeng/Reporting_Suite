import os
import sys

warnbasescovered = 90
warnsaturation = 1e-5
warnontarget = 80
warnmeancoverage = 40
warncoverageregion = 100
warncoveragethreshold = 6
warncoveragecorrelation = 0.90
warnstd = 0.3

offtargetoffset = 1000
offtargetthreshold = 15
TMP = '/tmp/'
CHR_LENGTHS = os.path.join(os.path.dirname(sys.argv[0]), 'chr_lengths_hg19.txt')

availablefeatures = ['percbases', 'saturation', 'specificity', 'coveragefreq', 'coveragedistr', 'coveragestd',
                    'gcbias', 'coveragecorr']