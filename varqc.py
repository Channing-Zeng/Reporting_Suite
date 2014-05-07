#!/usr/bin/env python
from os.path import join
import sys
from src.main import common_main, read_samples_info_and_split
from src.runner import run_all
from src.summarize import summarize_qc
try:
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Dumper

from src.quality_control import quality_control, check_quality_control_config
from src.my_utils import info


def main(args):
    config, options = common_main(
        'varqc',
        opts=[(['--var', '--vcf'], 'variants.vcf', {
            'dest': 'vcf',
            'help': 'variants to evaluate'}),
        ])

    if 'quality_control' in config:
        check_quality_control_config(config)

    samples = read_samples_info_and_split(config, options)

    try:
        run_all(config, samples, process_one, finalize_one, finalize_all)
    except KeyboardInterrupt:
        exit()


def process_one(cnf, vcf_fpath):
    qc_dir = join(cnf['output_dir'], 'qc')
    if 'quality_control' in cnf:
        return quality_control(cnf, qc_dir, vcf_fpath)
    else:
        return None, None


def finalize_one(cnf, qc_report_fpath, qc_plots_fpaths):
    if qc_report_fpath:
        info(cnf['log'], 'Saved QC report to ' + qc_report_fpath)
    if qc_plots_fpaths:
        info(cnf['log'], 'Saved QC plots are in: ' + ', '.join(qc_plots_fpaths))


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

