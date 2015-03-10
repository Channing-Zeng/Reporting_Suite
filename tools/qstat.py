#!/usr/bin/env python
import sys
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from itertools import repeat, izip, count
import subprocess
import sys


new_example = '''job-ID     prior   name       user         state submit/start at     queue                          jclass                         slots ja-task-ID
------------------------------------------------------------------------------------------------------------------------------------------------
      1157 1.25000 VA_BPmZGg= klpf990      r     03/10/2015 07:43:20 batch.q@orr                                                      21
      1159 1.25000 VA_BPmZGg= klpf990      r     03/10/2015 07:44:05 batch.q@chara                                                    21
      1158 0.00000 VQ_BPmZGg= klpf990      hqw   03/10/2015 07:41:51                                                                  21
      1160 0.00000 VQ_BPmZGg= klpf990      hqw   03/10/2015 07:41:51                                                                  21
      1161 0.00000 VQS_BPmZGg klpf990      hqw   03/10/2015 07:41:51                                                                   1
      1162 0.00000 VFS_BPmZGg klpf990      hqw   03/10/2015 07:41:51                                                                   3
      1163 0.00000 VQA_BPmZGg klpf990      hqw   03/10/2015 07:41:51                                                                  21
      1164 0.00000 VQA_BPmZGg klpf990      hqw   03/10/2015 07:41:51                                                                  21
      1165 0.00000 VQA_BPmZGg klpf990      hqw   03/10/2015 07:41:51                                                                  21
      1166 0.00000 VQA_BPmZGg klpf990      hqw   03/10/2015 07:41:51                                                                  21
      1167 0.00000 VQAS_BPmZG klpf990      hqw   03/10/2015 07:41:51                                                                   1
      1171 0.00000 CR_BPmZGg= klpf990      hqw   03/10/2015 07:41:52                                                                   1
'''


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

    if l[0] != ' ' and l.strip().split()[3] == 'klpf990':
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