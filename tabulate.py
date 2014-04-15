#!/usr/bin/env python


import sys


def dots_to_empty_cells(tsv_fpath):
    """Put dots instead of empty cells in order to view TSV with column -t
    """
    with open(tsv_fpath) as inp:
        for l in inp:
            while '\t\t' in l:
                l = l.replace('\t\t', '\t.\t')
            sys.stdout.write(l)


def main(args):
    dots_to_empty_cells(args[0])


if __name__ == '__main__':
    main(sys.argv[1:])