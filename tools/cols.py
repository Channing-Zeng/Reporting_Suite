#!/usr/bin/env python


import sys
from itertools import count


def cols(in_f):
    """Visially alignes columns of a tab-delimeted or a comma-delimited file
    :param in_f: input file stream
    """
    rows = []
    max_lens = None
    delimiter = '\t'
    for i, l in enumerate(in_f):
        if i == 0:
            if delimiter not in l and ',' in l:
                delimiter = ','
            row = l.split(delimiter)
            max_lens = [0 for _ in row[:-1]]
        row = l.split(delimiter)
        rows.append(row)
        max_lens = map(max, zip(max_lens, map(len, row[:-1])))
    for row in rows:
        if row:
            for v, max_len in zip(row[:-1], max_lens):
                sys.stdout.write(v + ' '*(max_len - len(v) + 2))
            sys.stdout.write(row[-1])
    sys.stdout.close()

def main():
    cols(open(sys.argv[1]) if len(sys.argv) > 1 else sys.stdin)


if __name__ == '__main__':
    main()
