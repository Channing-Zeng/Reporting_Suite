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

    for fpath in results_one:
        if fpath:
            try:
                verify_file(fpath)
            except:
                try:
                    for fpath in fpath:
                        verify_file(fpath)
                except:
                    pass

    return results_one

