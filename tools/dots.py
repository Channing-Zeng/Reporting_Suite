#!/usr/bin/env python


import sys


def dots_to_empty_cells():
    """Put dots instead of empty cells in order to view TSV with column -t
    """
    for l in sys.stdin:
        while '\t\t' in l:
            l = l.replace('\t\t', '\t.\t')
        sys.stdout.write(l)
    sys.stdout.close()

def main():
    dots_to_empty_cells()


if __name__ == '__main__':
    main()