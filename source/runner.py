#!/usr/bin/env python
import hashlib
import shutil
import os
from os.path import join, isdir, isfile

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper

from source.utils import critical, info, iterate_file, step_greetings, make_tmpdir


def run_all(cnf, cnfs_by_sample, required_inputs, optional_inputs,
            process_one, finalize_one, finalize_all):

    with make_tmpdir(cnf):

        if not cnfs_by_sample or len(cnfs_by_sample) == 1:
            if cnfs_by_sample:
                _, sample_cnf = cnfs_by_sample.items()[0]
            else:
                sample_cnf = cnf
            run_one(sample_cnf, required_inputs, optional_inputs, process_one, finalize_one)

        else:
            results = []
            if cnf.get('parallel'):
                try:
                    from joblib import Parallel, delayed
                except ImportError:
                    from joblib import Parallel, delayed
                    critical(
                        '\nERROR: Joblib not found. You may want samples to be processed '
                        'in parallel, in this case, make sure python joblib intalled. '
                        '(pip install joblib).')
                else:
                    from joblib import Parallel, delayed
                    for sample_name, sample_cnf in cnfs_by_sample.items():
                        sample_cnf['verbose'] = False

                    results = Parallel(n_jobs=len(cnfs_by_sample)) \
                        (delayed(run_one)(sample_cnf, required_inputs, optional_inputs,
                                          process_one, finalize_one,
                                          multiple_samples=True)
                            for sample_name, sample_cnf in cnfs_by_sample.items())
            else:
                results = []
                for sample_name, sample_cnf in cnfs_by_sample.items():
                    results.append(
                        run_one(sample_cnf, required_inputs, optional_inputs,
                                process_one, finalize_one,
                                multiple_samples=True))

            if cnfs_by_sample:
                info('')
                info('*' * 70)
                info('Results for each sample:')
                finalize_all(cnf, cnfs_by_sample, results)

    work_dirpath = cnf['work_dir']
    if not cnf.get('keep_intermediate') and isdir(work_dirpath):
        shutil.rmtree(work_dirpath)


def filter_rejected(cnf, input_fpath):
    step_greetings(cnf, 'Extracting dataset by filename, filtering REJECT line.')
    output_fpath = iterate_file(cnf, input_fpath,
                                (lambda l: l if 'REJECT' not in l else None),
                                cnf['work_dir'])
    info(cnf.get('log'), 'Saved to ' + output_fpath)
    return output_fpath


def run_one(cnf, required_inputs, optional_inputs,
            process_one_fun, finalize_one_fun, multiple_samples=False):
    input_fpaths = []
    for key in required_inputs:
        input_fpaths.append(cnf[key])
    for key in optional_inputs:
        if key in cnf:
            input_fpaths.append(cnf[key])

    if cnf.get('keep_intermediate') and not 'log' in cnf:
        cnf['log'] = join(cnf['work_dir'], cnf['name'] + '_log.txt')
        if isfile(cnf['log']):
            os.remove(cnf['log'])
    else:
        cnf['log'] = None

    if multiple_samples:
        info('')
        info('*' * 70)
        msg = '*' * 3 + ' Sample ' + cnf['name'] + ' '
        info(cnf.get('log'), msg + ('*' * (70 - len(msg)) if len(msg) < 70 else ''))

        info(cnf.get('log'), 'Input:')
        for fpath in input_fpaths:
            info(cnf.get('log'), '  ' + fpath)

    results_one = process_one_fun(cnf, *input_fpaths)

    info(cnf['log'], '')
    info(cnf['log'], '*' * 70)
    info(cnf['log'], cnf['name'])
    finalize_one_fun(cnf, *results_one)

    return results_one

