#!/usr/bin/env python
import hashlib
import shutil
import os
from os.path import join, isdir, isfile, realpath
from source.bcbio_utils import file_exists

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper

from source.utils import critical, info, iterate_file, step_greetings, make_tmpdir


def run_all(cnf, sample_cnfs_by_name, required_inputs, optional_inputs,
            process_one, finalize_one, finalize_all):

    with make_tmpdir(cnf):
        if len(sample_cnfs_by_name) == 1:
            sample_name, sample_cnf = sample_cnfs_by_name.items()[0]
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
                    for sample_name, sample_cnf in sample_cnfs_by_name.items():
                        sample_cnf['verbose'] = False

                    results = Parallel(n_jobs=len(sample_cnfs_by_name)) \
                        (delayed(run_one)(sample_cnf, required_inputs, optional_inputs,
                                          process_one, finalize_one,
                                          multiple_samples=True)
                            for sample_name, sample_cnf in sample_cnfs_by_name.items())
            else:
                results = []
                for sample_name, sample_cnf in sample_cnfs_by_name.items():
                    results.append(
                        run_one(sample_cnf, required_inputs, optional_inputs,
                                process_one, finalize_one,
                                multiple_samples=True))

            if sample_cnfs_by_name:
                info('')
                info('*' * 70)
                info('Results for each sample:')
                finalize_all(cnf, sample_cnfs_by_name, results)

    work_dirpath = cnf['work_dir']
    if not cnf.get('keep_intermediate') and isdir(work_dirpath):
        shutil.rmtree(work_dirpath)


def filter_ensemble(cnf, input_fpath):
    step_greetings(cnf, 'Extracting dataset by filename, filtering ensemble reject line.')
    output_fpath = iterate_file(cnf, input_fpath,
                                (lambda l: 'REJECT' not in l),
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

    if cnf.get('keep_intermediate'):
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

