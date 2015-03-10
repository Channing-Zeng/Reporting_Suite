#!/usr/bin/env python

import sys

with open(sys.argv[1]) as inp, open(sys.argv[2], 'w') as out:
    for l in inp:
        if l and not l.startswith('#') and l.startswith('chr'):
            l = l[3:]
        out.write(l)
