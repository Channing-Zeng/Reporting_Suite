#!/usr/bin/env python
from os import makedirs
from os.path import dirname, isdir, basename
import shutil

import sys
from source.utils_from_bcbio import file_exists
from source.logger import critical, err
from source.varqc.stats_gatk import gatk_qc
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from source.main import read_opts_and_cnfs, load_genome_resources, check_system_resources, check_inputs
from source.runner import run_one
from source.summarize import summarize_qc
from source.utils import info, verify_module, verify_file
from source.vcf import filter_rejected, extract_sample


def main(args):
    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--var', '--vcf'], dict(
                dest='vcf',
                help='variants to evaluate')
             ),
        ],
        required_keys=['vcf'],
        file_keys=['vcf'],
        key_for_sample_name='vcf')

    check_system_resources(cnf,
        required=['java', 'gatk', 'snpeff', 'bgzip', 'tabix'],
        optional=['bcftools', 'plot_vcfstats'])

    load_genome_resources(cnf,
        required=['seq', 'dbsnp'],
        optional=['cosmic', '1000genomes'])

    check_quality_control_config(cnf)

    info('Using variants ' + cnf['vcf'])

    run_one(cnf, process_one, finalize_one)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


if verify_module('matplotlib'):
    import matplotlib
    matplotlib.use('Agg')  # non-GUI backend
    from source.varqc.distribution_plots import variants_distribution_plot
    from source.varqc.stats_bcftools import bcftools_qc
else:
    info('Warning: matplotlib is not installed, cannot draw plots.')


def check_quality_control_config(cnf):
    qc_cnf = cnf['quality_control']

    to_exit = False
    dbs_dict = {}
    for db in qc_cnf['databases']:
        if not db:
            err('Empty field for quality_control databases')
            to_exit = True
        elif file_exists(db):
            if not verify_file(db, 'VCF'):
                to_exit = True
            dbs_dict[basename(db)] = db
        elif db not in cnf.genome:
            to_exit = True
            err(db + ' for variant qc is not found in genome resources in system config.')
        else:
            dbs_dict[db] = cnf['genome'][db]

    if to_exit:
        exit()

    qc_cnf['database_vcfs'] = dbs_dict

    ## FOR SUMMARIZING ##
    # if 'summary_output' in qc_cnf or 'qc_summary_output' in cnf:
    #     qc_output_fpath = qc_cnf.get('summary_output') or\
    #                       cnf.get('qc_summary_output')
    #     summary_output_dir = dirname(qc_output_fpath)
    #     if not isdir(summary_output_dir):
    #         try:
    #             makedirs(summary_output_dir)
    #         except OSError:
    #             critical('ERROR: cannot create directory for '
    #                      'qc summary report: ' + summary_output_dir)
    #     if not verify_dir(summary_output_dir, 'qc_summary_output'):
    #         exit()


def process_one(cnf):
    vcf_fpath = cnf['vcf']

    if cnf.get('filter_reject'):
        vcf_fpath = filter_rejected(cnf, vcf_fpath)

    if cnf.get('extract_sample'):
        vcf_fpath = extract_sample(cnf, vcf_fpath, cnf['name'])

    qc_report_fpath = gatk_qc(cnf, vcf_fpath)

    if verify_module('matplotlib'):
        qc_plots_fpaths = bcftools_qc(cnf, vcf_fpath)
        qc_var_distr_plot_fpath = variants_distribution_plot(cnf, vcf_fpath)
        return qc_report_fpath, [qc_var_distr_plot_fpath] + qc_plots_fpaths
    else:
        return qc_report_fpath, None


def finalize_one(cnf, qc_report_fpath, qc_plots_fpaths):
    if qc_report_fpath:
        info('Saved QC report to ' + qc_report_fpath)
    if qc_plots_fpaths:
        info('Saved QC plots are in: ' + ', '.join(qc_plots_fpaths))
    elif not verify_module('matplotlib'):
        info('Warning: QC plots were not generated because matplotlib is not installed.')


def finalize_all(cnf, samples, results):
    for (sample_name, cnf), (qc_dir, qc_report, qc_plots) in zip(samples.items(), results):
        if qc_dir:
            info(sample_name + ':')
            info('  ' + qc_report)
            info('  ' + qc_dir)

    qc_cnf = cnf.get('quality_control')
    if qc_cnf and 'summary_output' in qc_cnf or 'qc_summary_output' in cnf:
        qc_output_fpath = cnf.get('qc_summary_output') or qc_cnf.get('summary_output')
        summarize_qc([rep for _, _, _, rep, _ in results], qc_output_fpath)
        info('Variant QC summary:')
        info('  ' + qc_output_fpath)


if __name__ == '__main__':
    main(sys.argv[1:])

