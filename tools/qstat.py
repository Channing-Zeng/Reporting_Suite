#!/usr/bin/env python
import sys
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from itertools import repeat, izip, count
import subprocess
import sys

f = subprocess.Popen(['qstat', '-r'], stdout=subprocess.PIPE).stdout
# f = open('/Users/vladsaveliev/vagrant/reporting_suite/test/qstat')

rows = []

cur_tokens = None
for i, l in enumerate(f):
    if i == 0:
        rows.append(l.split())
        continue

    if i == 1:
        continue

    if l[0] != ' ' and cur_tokens.strip().split()[3] == 'klpf990':
        cur_tokens = l.split()

    elif cur_tokens and l.strip().startswith('Full jobname:'):
        full_name = l.strip().split()[2]

        if len(cur_tokens) == 8:
            rows.append(cur_tokens[:2] + [full_name] + cur_tokens[3:7] + [''] + cur_tokens[7:])
        else:
            rows.append(cur_tokens[:2] + [full_name] + cur_tokens[3:])
        cur_tokens = None

col_widths = repeat(0)

for row in rows:
    col_widths = [max(len(v), w) for v, w in izip(row, col_widths)]

for i, row in enumerate(rows):
    line = ''
    for val, w in izip(row, col_widths):
        cell = val + (' ' * (w - len(val) + 2))
        sys.stdout.write(cell)
        line += cell
    sys.stdout.write('\n')
    if i == 0:
        sys.stdout.write('-' * len(line[:-2]) + '\n')