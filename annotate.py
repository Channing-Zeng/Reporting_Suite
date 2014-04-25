#!/usr/bin/env python
import sys
import shutil
import os
from os.path import join, realpath, isdir, isfile, dirname, relpath

from yaml import dump
from src.config import process_config
from src.tsv import make_tsv
try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper

from src.annotation import run_annotators
from src.quality_control import bcftools_qc, gatk_qc, quality_control
from src.my_utils import critical, info


def main(args):
    if len(args) < 1:
        exit('Usage: python ' + __file__ + ' run_info.yaml\n'
             '    or python ' + __file__ + ' system_info.yaml run_info.yaml')

    if (2, 5) > sys.version_info[:2] >= (3, 0):
        exit('Python version 2.5 and higher is supported (you are running ' +
             '.'.join(map(str, sys.version_info[:3])) + ')\n')

    if len(args) == 1:
        run_config_path = args[0]
        system_config_path = join(dirname(realpath(__file__)),
                                  'system_info_rask.yaml')
        sys.stderr.write('Notice: using system_info_rask.yaml as a default'
                         ' tools configutation file.\n\n')
    else:
        system_config_path = args[0]
        run_config_path = args[1]

    if not os.path.isfile(system_config_path):
        exit(system_config_path + ' does not exist or is a directory.\n')
    if not os.path.isfile(run_config_path):
        exit(run_config_path + ' does not exist or is a directory.\n')

    to_exit = False
    if not system_config_path.endswith('.yaml'):
        sys.stderr.write(system_config_path + ' does not end with .yaml,'
                                              ' maybe incorrect parameter?\n')
        to_exit = True
    if not run_config_path.endswith('.yaml'):
        sys.stderr.write(run_config_path + ' does not end with .yaml,'
                                           ' maybe incorrect parameter?\n')
        to_exit = True
    if to_exit:
        exit()

    config, samples = process_config(system_config_path, run_config_path)
    try:
        annotate(samples, config.get('parallel'))
    except KeyboardInterrupt:
        exit()


def annotate(samples, parallel=False):
    if len(samples) == 1:
        sample_name, sample_cnf = samples.items()[0]
        annotate_one(sample_cnf)
    else:
        results = []
        if parallel:
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
                    (delayed(annotate_one)(sample_cnf, multiple_samples=True)
                        for sample_name, sample_cnf in samples.items())
        else:
            for sample_name, sample_cnf in samples.items():
                results.append(
                    annotate_one(sample_cnf, multiple_samples=True))

        info('')
        info('*' * 70)
        info('Results for each sample')
        for (sample_name, cnf), (vcf, tsv, qc) in zip(samples.items(), results):
            info(cnf['log'], sample_name + ':')
            info(cnf['log'], '  ' + vcf)
            info(cnf['log'], '  ' + tsv)
            if qc:
                info(cnf['log'], '  qc: ' + qc)

    for name, data in samples.items():
        work_dirpath = data['work_dir']
        tx_dirpath = join(work_dirpath, 'tx')

        if isdir(tx_dirpath):
            shutil.rmtree(tx_dirpath)

        if not data.get('keep_intermediate') \
                and isdir(work_dirpath):
            shutil.rmtree(work_dirpath)


def annotate_one(cnf, multiple_samples=False):
    if cnf.get('keep_intermediate'):
        cnf['log'] = join(cnf['work_dir'], cnf['name'] + '_log.txt')
        if isfile(cnf['log']):
            os.remove(cnf['log'])

        with open(join(cnf['work_dir'], cnf['name'] + '_config.yaml'), 'w') as f:
            f.write(dump(cnf, Dumper=Dumper))
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

    final_vcf_fpath = run_annotators(cnf, cnf['vcf'])

    final_tsv_fpath = make_tsv(cnf, final_vcf_fpath)

    qc_report_fpath = None
    qc_plots_fpaths = None
    qc_dir = join(cnf['output_dir'], cnf['name'] + '_qc')
    if 'quality_control' in cnf:
        qc_report_fpath, qc_plots_fpaths = quality_control(cnf, qc_dir, final_vcf_fpath)

    final_vcf_fpath = relpath(final_vcf_fpath, os.getcwd())
    final_tsv_fpath = relpath(final_tsv_fpath, os.getcwd())
    qc_report_fpath = relpath(qc_report_fpath, os.getcwd())
    qc_plots_fpaths = [relpath(fpath, os.getcwd()) for fpath in qc_plots_fpaths]

    info(cnf['log'], '')
    info(cnf['log'], '*' * 70)
    info(cnf['log'], cnf['name'])
    info(cnf['log'], 'Saved final VCF to ' + final_vcf_fpath)
    info(cnf['log'], 'Saved final TSV to ' + final_tsv_fpath)
    if qc_report_fpath:
        info(cnf['log'], 'Saved quality control report to ' + qc_report_fpath)
    if qc_plots_fpaths:
        info(cnf['log'], 'Saved quality control plots: ' + ', '.join(qc_plots_fpaths))

    return final_vcf_fpath, final_tsv_fpath, qc_dir


if __name__ == '__main__':
    main(sys.argv[1:])