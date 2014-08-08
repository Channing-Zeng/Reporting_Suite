from collections import OrderedDict
import json
from source.logger import info
from source.reporting import summarize, write_summary_reports, get_per_sample_fpaths_for_bcbio_final_dir, Metric, \
    Record, parse_value
from source.utils import OrderedDefaultDict
from source.variants.qc_gatk import gatk_metrics, varqc_json_ending


class VariantCaller:
    def __init__(self, suf):
        self.name = suf
        self.suf = suf
        self.single_qc_rep_fpaths = []
        self.summary_qc_report = None
        self.summary_qc_rep_fpaths = []


def make_summary_reports(cnf, sample_names):
    varqc_dir = cnf['base_name']

    vcf_sufs = cnf['vcf_suf'].split(',')
    callers = [VariantCaller(suf) for suf in vcf_sufs]

    for caller in callers:
        fpaths, sample_names = get_per_sample_fpaths_for_bcbio_final_dir(
            cnf['bcbio_final_dir'], sample_names, varqc_dir,
            '-' + caller.suf + varqc_json_ending)
        if fpaths:
            caller.single_qc_rep_fpaths = fpaths

    if len(callers) > 1:
        _make_for_multiple_variant_callers(callers, cnf, sample_names)

    else:
        _make_for_single_variant_caller(callers, cnf, sample_names)


def _make_for_single_variant_caller(callers, cnf, sample_names):
    reports = summarize(sample_names, callers[0].single_qc_rep_fpaths, _parse_qc_sample_report)
    # reports = [records=[], records=[], records=[]...]

    full_summary_fpaths = write_summary_reports(
        cnf['output_dir'], cnf['work_dir'], reports, 'varQC', 'Variant QC')

    info()
    info('*' * 70)
    for fpath in full_summary_fpaths:
        info(fpath)


def _make_for_multiple_variant_callers(callers, cnf, sample_names):
    for caller in callers:
        caller.summary_qc_report = summarize(
            sample_names, caller.single_qc_rep_fpaths, _parse_qc_sample_report)

        caller.summary_qc_rep_fpaths = write_summary_reports(
            cnf['output_dir'], cnf['work_dir'], caller.summary_qc_report,
            caller.suf + '.varQC', 'Variant QC for ' + caller.name)

    all_single_reports = [r for c in callers for r in c.single_qc_rep_fpaths]
    all_sample_names = [sample_name + '-' + c.suf for sample_name in sample_names for c in callers]

    full_summary_report = summarize(all_sample_names, all_single_reports, _parse_qc_sample_report)

    full_summary_fpaths = write_summary_reports(
        cnf['output_dir'], cnf['work_dir'], full_summary_report, 'varQC', 'Variant QC')

    info()
    info('*' * 70)

    for caller in callers:
        info(caller.name)
        for fpath in caller.summary_qc_rep_fpaths:
            info('  ' + fpath)
        info()

    info('Total')
    for fpath in full_summary_fpaths:
        info('  ' + fpath)


def _parse_qc_sample_report(json_fpath):
    return Record.load_records(json_fpath)

