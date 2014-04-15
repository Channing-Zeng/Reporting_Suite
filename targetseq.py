#!/usr/bin/env python

import sys
from os.path import join, splitext, basename, realpath, isfile, getsize, dirname, exists

# Number of mapped reads	1,625,036
# Number of reads on target	1,021,939
# Percent reads on target	62.89%
# Percent reads on padded target	89.93%
# Total aligned base reads	61,749,848
# Total base reads on target	33,727,151
# Percent base reads on target	54.62%
# Bases in targeted reference	1,140,710
# Bases covered (at least 1x)	1,052,115
# Average base coverage depth	29.57
# Maximum base read depth	307
# Average base read depth	32.06
# Std.Dev base read depth	27.31
# Target coverage at 1x	92.233%
# Target coverage at 5x	83.661%
# Target coverage at 10x	72.896%
# Target coverage at 20x	54.553%
# Target coverage at 50x	19.890%
# Target coverage at 100x	2.552%
# Target coverage at 500x	0.000%
# Target coverage at 1000x	0.000%
# Target coverage at 2500x	0.000%
# Target coverage at 5000x	0.000%
# Target coverage at 10000x	0.000%


def main(args):
    bam = args[0]
    bed = args[1]
    ref = args[2]
    pad = args[3] if len(args) > 2 else 500


if __name__ == '__main__':
    main(sys.argv[1:])