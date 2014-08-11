from collections import defaultdict
from source.bcbio_structure import BCBioStructure, VariantCaller
from source.logger import info
from source.reporting import summarize, write_summary_reports, get_per_sample_fpaths_for_bcbio_final_dir, Record
from source.variants.qc_gatk import varqc_json_ending


def make_summary_reports(cnf, bcbio_structure):
    if len(bcbio_structure.variant_callers) > 1:
        _make_for_multiple_variant_callers(cnf, bcbio_structure.variant_callers.values())

    else:
        _make_for_single_variant_caller(cnf, bcbio_structure.variant_callers.values()[0])


def _make_for_single_variant_caller(cnf, caller):
    reports = summarize(caller.get_qc_reports_by_samples(), _parse_qc_sample_report)

    full_summary_fpaths = write_summary_reports(
        cnf.output_dir,
        cnf.work_dir,
        reports,
        'varQC',
        'Variant QC')

    info()
    info('*' * 70)
    for fpath in full_summary_fpaths:
        info(fpath)


def _make_for_multiple_variant_callers(cnf, callers):
    for caller in callers:
        caller.summary_qc_report = summarize(
            caller.get_qc_reports_by_samples(),
            _parse_qc_sample_report)

        caller.summary_qc_rep_fpaths = write_summary_reports(
            cnf.output_dir,
            cnf.work_dir,
            caller.summary_qc_report,  # TODO outputdir - read from structure too?
            caller.suf + '.varQC',
            'Variant QC for ' + caller.name)

    all_single_reports = dict((sname + '-' + c.suf, rep)
                              for c in callers
                              for sname, rep in c.get_qc_reports_by_samples().items())

    full_summary_report = summarize(
        all_single_reports,
        _parse_qc_sample_report)

    full_summary_fpaths = write_summary_reports(
        cnf.output_dir,
        cnf.work_dir,
        full_summary_report,
        'varQC',
        'Variant QC')

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

