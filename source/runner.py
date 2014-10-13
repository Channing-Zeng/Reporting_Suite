#!/usr/bin/env python
import hashlib
import shutil
import os
from os.path import join, isdir, isfile
from source.file_utils import verify_file
from source.logger import step_greetings, info
from source.file_utils import file_exists

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper


# def run_all(cnf, cnfs_by_sample, required_inputs, optional_inputs,
#             process_one, finalize_one, finalize_all):
#
#     with tmpdir(cnf):
#
#         if not cnfs_by_sample or len(cnfs_by_sample) == 1:
#             if cnfs_by_sample:
#                 _, sample_cnf = cnfs_by_sample.items()[0]
#             else:
#                 sample_cnf = cnf
#             run_one(sample_cnf, required_inputs, optional_inputs, process_one, finalize_one)
#
#         else:
#             results = []
#             if cnf.get('parallel'):
#                 try:
#                     from joblib import Parallel, delayed
#                 except ImportError:
#                     from joblib import Parallel, delayed
#                     critical(
#                         '\nERROR: Joblib not found. You may want samples to be processed '
#                         'in parallel, in this case, make sure python joblib intalled. '
#                         '(pip install joblib).')
#                 else:
#                     from joblib import Parallel, delayed
#                     for sample_name, sample_cnf in cnfs_by_sample.items():
#                         sample_cnf['verbose'] = False
#
#                     results = Parallel(n_jobs=len(cnfs_by_sample)) \
#                         (delayed(run_one)(sample_cnf, required_inputs, optional_inputs,
#                                           process_one, finalize_one,
#                                           multiple_samples=True)
#                             for sample_name, sample_cnf in cnfs_by_sample.items())
#             else:
#                 results = []
#                 for sample_name, sample_cnf in cnfs_by_sample.items():
#                     results.append(
#                         run_one(sample_cnf, required_inputs, optional_inputs,
#                                 process_one, finalize_one,
#                                 multiple_samples=True))
#
#             if cnfs_by_sample:
#                 info('')
#                 info('*' * 70)
#                 info('Results for each sample:')
#                 finalize_all(cnf, cnfs_by_sample, results)
#
#     if not cnf.get('keep_intermediate') and isdir(cnf['work_dir']):
#         shutil.rmtree(cnf['work_dir'])


def run_one(cnf, process_one_fun, finalize_one_fun=None, multiple_samples=False):
    # input_fpaths = []
    # for key in required_inputs:
    #     input_fpaths.append(cnf[key])
    # for key in optional_inputs:
    #     if key in cnf:
    #         input_fpaths.append(cnf[key])

    if multiple_samples:
        info('')
        info('*' * 70)
        msg = '*' * 3 + ' Sample ' + cnf['name'] + ' '
        info(msg + ('*' * (70 - len(msg)) if len(msg) < 70 else ''))

        # info('Input:')
        # for fpath in input_fpaths:
        #     info('  ' + fpath)

    results_one = process_one_fun(cnf)

    if finalize_one_fun and results_one:
        info('')
        info('*' * 70)
        if multiple_samples:
            info(cnf['name'])

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

