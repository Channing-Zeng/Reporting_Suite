#!/usr/bin/env python

import sys
from source.vcf_read import read_samples_info_and_split
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from source.main import common_main, load_genome_resources, check_system_resources
from source.runner import run_all
from source.summarize import summarize_qc
from source.varqc import qc
from source.utils import info, verify_module


def main(args):
    required = ['vcf']
    optional = []

    config, options = common_main(
        'varqc',
        extra_opts=[(['--var', '--vcf'], 'variants.vcf', {
            'dest': 'vcf',
            'help': 'variants to evaluate'}),
        ],
        required=required)

    check_system_resources(config, ['java', 'gatk', 'snpeff', 'bcftools', 'plot_vcfstats', 'bgzip', 'tabix'])
    load_genome_resources(config, ['seq', 'dbsnp'])

    if 'quality_control' in config:
        qc.check_quality_control_config(config)

    sample_cnfs_by_name = read_samples_info_and_split(config, options, required + optional)

    try:
        run_all(config, sample_cnfs_by_name, required, optional,
                process_one, finalize_one, finalize_all)
    except KeyboardInterrupt:
        exit()


def process_one(cnf, vcf_fpath):
    if 'quality_control' in cnf:
        return qc.run_qc(cnf, cnf['output_dir'], vcf_fpath)
    else:
        return None, None


def finalize_one(cnf, qc_report_fpath, qc_plots_fpaths):
    if qc_report_fpath:
        info(cnf['log'], 'Saved QC report to ' + qc_report_fpath)
    if qc_plots_fpaths:
        info(cnf['log'], 'Saved QC plots are in: ' + ', '.join(qc_plots_fpaths))
    elif not verify_module('matplotlib'):
        info('Warning: QC plots were not generated because matplotlib is not installed.')


def finalize_all(cnf, samples, results):
    for (sample_name, cnf), (qc_dir, qc_report, qc_plots) in zip(samples.items(), results):
        if qc_dir:
            info(cnf['log'], sample_name + ':')
            info(cnf['log'], '  ' + qc_report)
            info(cnf['log'], '  ' + qc_dir)

    qc_cnf = cnf.get('quality_control')
    if qc_cnf and 'summary_output' in qc_cnf or 'qc_summary_output' in cnf:
        qc_output_fpath = cnf.get('qc_summary_output') or qc_cnf.get('summary_output')
        summarize_qc([rep for _, _, _, rep, _ in results], qc_output_fpath)
        info(cnf['log'], 'Variant QC summary:')
        info(cnf['log'], '  ' + qc_output_fpath)


if __name__ == '__main__':
    main(sys.argv[1:])

