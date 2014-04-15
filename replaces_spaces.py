#!/usr/bin/env python


import sys


def spaces(fpath):
    """Put dots instead of empty cells in order to view TSV with column -t
    """
    with open(fpath) as inp:
        for l in inp:
            if not l.startswith('#'):
                l = l.replace(' ', '_')
            sys.stdout.write(l)


def main(args):
    spaces(args[0])


if __name__ == '__main__':
    main(sys.argv[1:])