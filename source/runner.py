#!/usr/bin/env python
from source.file_utils import verify_file
from source.logger import info

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper


def run_one(cnf, process_one_fun, finalize_one_fun=None, *args, **kwargs):
    results_one = process_one_fun(cnf, *args, **kwargs)

    if finalize_one_fun and results_one:
        info('')
        info('*' * 70)
        finalize_one_fun(cnf, *results_one)

    for fpaths in results_one:
        if fpaths:
            ok = True
            info('Checking expected results...')
            if isinstance(fpaths, basestring):
                fpaths = [fpaths]
            for fpath in fpaths:
                if not verify_file(fpath):
                    ok = False
            if ok:
                info('The results are good.')

    return results_one

