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
    line_num = 0
    for l in in_f:
        if not l.strip() or l.strip().startswith('##'):
            sys.stdout.write(l)
            continue
        line_num += 1
        l = l.replace('\n', '')  # removing the trailing '\n'
        if line_num == 1:
            if delimiter not in l and ',' in l:
                delimiter = ','
            row = l.split(delimiter)
            max_lens = [0 for _ in row[:-1]]
        else:
            row = l.split(delimiter)
            if len(row) > len(max_lens)+1:
                for i in range(len(max_lens)+1, len(row)):
                    max_lens.append(len(row[i]))
            #     sys.stdout.flush()
            #     sys.stderr.write('Error: line #' + str(line_num) + ' is longer than the first line ' +
            #                      '(has ' + str(len(row)) + ' columns instead of ' + str(len(max_lens) + 1) + '):\n' +
            #                      '  ' + l + '\n')
            #     sys.stderr.flush()
            #     continue
            if len(row) < len(max_lens)+1:  # adding empty columns to a shorter row to make it the same size as the rest
                for _ in range(len(max_lens)+1 - len(row)):
                    row.append('')
        max_lens = map(max, zip(max_lens, map(len, row[:-1])))
        rows.append(row)
    for row in rows:
        if row:
            for v, max_len in zip(row[:-1], max_lens):
                if not sys.stdout.closed:
                    sys.stdout.write(v + ' '*(max_len - len(v) + 2))
            if not sys.stdout.closed:
                sys.stdout.write(row[-1] + '\n')
    sys.stdout.close()

def main():
    try:
        cols(open(sys.argv[1]) if len(sys.argv) > 1 else sys.stdin)
    except IOError:
        pass


if __name__ == '__main__':
    main()
