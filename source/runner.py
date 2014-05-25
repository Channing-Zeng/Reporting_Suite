#!/usr/bin/env python
import shutil
import os
from os.path import join, isdir, isfile
from yaml import dump

try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper

from source.utils import critical, info, iterate_file, step_greetings


def run_all(cnf, samples, process_one, finalize_one, finalize_all):
    if len(samples) == 1:
        sample_name, sample_cnf = samples.items()[0]
        run_one(sample_cnf, process_one, finalize_one)
    else:
        results = []
        if cnf.get('parallel'):
            try:
                from joblib import Parallel, delayed
            except ImportError:
                critical(
                    '\nERROR: Joblib not found. You may want samples to be processed '
                    'in parallel, in this case, make sure python joblib intalled. '
                    '(pip install joblib).')
            else:
                for sample_name, sample_cnf in samples.items():
                    sample_cnf['verbose'] = False

                results = Parallel(n_jobs=len(samples)) \
                    (delayed(run_one)(sample_cnf, process_one, finalize_one,
                                      multiple_samples=True)
                        for sample_name, sample_cnf in samples.items())
        else:
            results = []
            for sample_name, sample_cnf in samples.items():
                results.append(
                    run_one(sample_cnf, process_one, finalize_one,
                            multiple_samples=True))

        if samples:
            info('')
            info('*' * 70)
            info('Results for each sample:')
            finalize_all(cnf, samples, results)

    # Cleaning
    for name, data in samples.items():
        work_dirpath = data['work_dir']
        tx_dirpath = join(work_dirpath, 'tx')

        if isdir(tx_dirpath):
            shutil.rmtree(tx_dirpath)

        if not data.get('keep_intermediate') \
                and isdir(work_dirpath):
            shutil.rmtree(work_dirpath)


def _filter_ensemble(cnf, input_fpath):
    step_greetings(cnf, 'Extracting dataset by filename, filtering ensemble reject line.')
    output_fpath = iterate_file(cnf, input_fpath,
                                (lambda l: 'REJECT' not in l),
                                cnf['work_dir'])
    info(cnf.get('log'), 'Saved to ' + output_fpath)
    return output_fpath


def run_one(cnf, process_one_fun, finalize_one_fun, multiple_samples=False):
    if cnf.get('keep_intermediate'):
        cnf['log'] = join(cnf['work_dir'], cnf['name'] + '_log.txt')
        if isfile(cnf['log']):
            os.remove(cnf['log'])
        #with open(join(cnf['work_dir'], cnf['name'] + '_config.yaml'), 'w') as f:
        #    f.write(dump(cnf, Dumper=Dumper))
    else:
        cnf['log'] = None
    # sample['fields'] = []

    if multiple_samples:
        info('')
        info('*' * 70)
        msg = '*' * 3 + ' Sample ' + cnf['name'] + ' '
        info(cnf.get('log'), msg + ('*' * (70 - len(msg)) if len(msg) < 70 else ''))
        info(cnf.get('log'), 'VCF: ' + cnf['vcf'])
        if cnf.get('bam'):
            info(cnf.get('log'), 'BAM: ' + cnf['bam'])

    vcf_fpath = cnf['vcf']
    if cnf.get('ensemble'):
        vcf_fpath = _filter_ensemble(cnf, vcf_fpath)

    results_one = process_one_fun(cnf, vcf_fpath)

    info(cnf['log'], '')
    info(cnf['log'], '*' * 70)
    info(cnf['log'], cnf['name'])
    finalize_one_fun(cnf, *results_one)

    return results_one

