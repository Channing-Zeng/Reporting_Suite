#!/usr/bin/env python


import sys


def dots_to_empty_cells(in_f):
    """Put dots instead of empty cells in order to view TSV with column -t
    :param in_f: input file stream
    """
    for l in in_f:
        while '\t\t' in l:
            l = l.replace('\t\t', '\t.\t')
        sys.stdout.write(l)
    sys.stdout.close()

def main():
    dots_to_empty_cells(open(sys.argv[1]) if len(sys.argv) > 1 else sys.stdin)


if __name__ == '__main__':
    main()
