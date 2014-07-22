#!/usr/bin/env python
import sys
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from itertools import repeat, izip, count
import subprocess
import sys

p = subprocess.Popen(['qstat', '-r'], stdout=subprocess.PIPE)

rows = []

cur_tokens = None
for i, l in enumerate(p.stdout):
    if i == 0:
        rows.append(l.split())
    if i == 1:
        continue

    if not l[0] != ' ' and len(l.split()) == 9:
        cur_tokens = l.split()

    elif cur_tokens and l.strip().startswith('Full jobname:'):
        full_name = l.strip().split()[1]
        rows.append(cur_tokens[:2] + [full_name] + cur_tokens[3:])
        cur_tokens = None

col_widths = repeat(0)

for row in rows:
    col_widths = [max(len(v), w) for v, w in izip(row, col_widths)]

for row in rows:
    for i, val, w in izip(count(), row, col_widths):
        line = val + (' ' * (w - len(val) + 2))
        sys.stdout.write(line)
        if i == 0:
            sys.stdout.write('-' * len(line))
    sys.stdout.write('\n')