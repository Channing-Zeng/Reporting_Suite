#!/usr/bin/env python

import __check_python_version

import sys
import shutil
from os.path import abspath, dirname, realpath, join, basename, relpath
import source
from source import SingleSample
from source.file_utils import verify_module, verify_file
from source.file_utils import file_exists
from source.logger import err, info, warn, send_email, critical
from source.prepare_args_and_cnf import check_genome_resources, check_system_resources
from source.variants import qc
from source.main import read_opts_and_cnfs
from source.runner import run_one
from source.variants.vcf_processing import remove_rejected
from source.reporting import SampleReport


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
        key_for_sample_name='vcf',
        proc_name=source.varqc_name,
    )

    check_system_resources(cnf)

    check_genome_resources(cnf)

    info('Using variants ' + cnf['vcf'])

    run_one(cnf, process_one, finalize_one)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


if verify_module('matplotlib'):
    import matplotlib
    matplotlib.use('Agg')  # non-GUI backend
    from source.variants.qc_plots import draw_plots
else:
    warn('Warning: matplotlib is not installed, cannot draw plots.')


def process_one(cnf):
    vcf_fpath = cnf['vcf']
    if not verify_file(vcf_fpath):
        critical('Annotated VCF ' + vcf_fpath + ' does not exist, thus cannot run VarQC')

    sample = SingleSample(cnf.sample, cnf.output_dir, vcf=cnf.vcf, bam=cnf.bam, genome=cnf.genome)

    if cnf.get('filter_reject'):
        vcf_fpath = remove_rejected(cnf, vcf_fpath)

    report = qc.make_report(cnf, vcf_fpath, sample)

    if verify_module('matplotlib'):
        qc_plots_fpaths = draw_plots(cnf, vcf_fpath)
    else:
        qc_plots_fpaths = []

    qc_plots_for_html_report_fpaths = qc_plots_fpaths
    report.plots = [relpath(plot_fpath, cnf.output_dir) for plot_fpath in qc_plots_for_html_report_fpaths]

    summary_report_html_fpath = report.save_html(
        cnf.output_dir, cnf.sample + '-' + cnf.caller + '.' + cnf.proc_name,
        caption='Variant QC for ' + cnf.sample + ' (caller: ' + cnf.caller + ')')

    return summary_report_html_fpath, qc_plots_fpaths


def finalize_one(cnf, qc_report_fpath, qc_plots_fpaths):
    if qc_report_fpath:
        info('Saved QC report to ' + qc_report_fpath)
    if qc_plots_fpaths:
        info('Saved QC plots are in: ' + ', '.join(qc_plots_fpaths))
    elif not verify_module('matplotlib'):
        warn('Warning: QC plots were not generated because matplotlib is not installed.')


def finalize_all(cnf, samples, results):
    for (sample_name, cnf), (qc_dir, qc_report, qc_plots) \
            in zip(samples.items(), results):
        if qc_dir:
            info(sample_name + ':')
            info('  ' + qc_report)
            info('  ' + qc_dir)

    qc_cnf = cnf.get('quality_control')
    if qc_cnf and 'summary_output' in qc_cnf or 'qc_summary_output' in cnf:
        qc_output_fpath = cnf.get('qc_summary_output') or qc_cnf.get('summary_output')
        # summarize_qc([rep for _, _, _, rep, _ in results], qc_output_fpath)
        info('Variant QC summary:')
        info('  ' + qc_output_fpath)


if __name__ == '__main__':
    main(sys.argv[1:])

