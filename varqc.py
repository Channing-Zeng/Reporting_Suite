#!/usr/bin/env python

import sys
from source.logger import err, critical
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from source.main import read_opts_and_cnfs, load_genome_resources, check_system_resources
from source.runner import run_all, run_one
from source.summarize import summarize_qc
from source.varqc import qc
from source.utils import info, verify_module


def main(args):
    required_keys = ['vcf']
    optional_keys = []

    cnf = read_opts_and_cnfs(
        extra_opts=[(['--var', '--vcf'], 'variants.vcf', {
            'dest': 'vcf',
            'help': 'variants to evaluate'}),
        ],
        required_keys=required_keys,
        optional_keys=optional_keys)

    check_system_resources(
        cnf,
        required=['java', 'gatk', 'snpeff', 'bgzip', 'tabix'],
        optional=['bcftools', 'plot_vcfstats'])
    load_genome_resources(
        cnf,
        required=['seq', 'dbsnp'],
        optional=['cosmic', '1000genomes'])

    if 'quality_control' not in cnf:
        critical('No quality_control section in the report, cannot run quality control.')

    qc.check_quality_control_config(cnf)

    run_one(cnf, required_keys, optional_keys, process_one, finalize_one)


def process_one(cnf, vcf_fpath):
    if 'quality_control' in cnf:
        return qc.run_qc(cnf, cnf['output_dir'], vcf_fpath)
    else:
        return None, None


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

