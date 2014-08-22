from collections import defaultdict, OrderedDict
from os.path import join
from source.bcbio_structure import Sample
from source.logger import info
from source.reporting import summarize, write_summary_reports, Record


def make_summary_reports(cnf, bcbio_structure):
    if len(bcbio_structure.variant_callers) > 1:
        _make_for_multiple_variant_callers(cnf, bcbio_structure.variant_callers.values())

    else:
        _make_for_single_variant_caller(cnf, bcbio_structure.variant_callers.values()[0])


def _make_for_single_variant_caller(cnf, caller):
    reports = summarize(cnf, caller.get_report_fpaths_by_sample(cnf), _parse_qc_sample_report, '')

    full_summary_fpaths = write_summary_reports(
        cnf.output_dir,
        cnf.work_dir,
        reports,
        base_fname=join(cnf.output_dir, cnf.name),
        caption='Variant QC')

    info()
    info('*' * 70)
    for fpath in full_summary_fpaths:
        info(fpath)


def _make_for_multiple_variant_callers(cnf, callers):
    for caller in callers:
        caller.summary_qc_report = summarize(
            cnf,
            caller.get_report_fpaths_by_sample(cnf),
            _parse_qc_sample_report,
            'Variant QC')

        caller.summary_qc_rep_fpaths = write_summary_reports(
            cnf.output_dir,
            cnf.work_dir,
            caller.summary_qc_report,  # TODO outputdir - read from structure too?
            join(cnf.output_dir, caller.suf + '.' + cnf.name),
            'Variant QC for ' + caller.name)

    all_single_reports = OrderedDict(
        ((Sample(s + '-' + c), rep)
         for (s, c, rep) in sorted(
            (s.name, c.name, rep)
            for c in callers
            for s, rep in c.get_report_fpaths_by_sample(cnf).items()
        )))

    full_summary_report = summarize(
        cnf,
        all_single_reports,
        _parse_qc_sample_report,
        'Variant QC')

    full_summary_fpaths = write_summary_reports(
        cnf.output_dir,
        cnf.work_dir,
        full_summary_report,
        join(cnf.output_dir, cnf.name),
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


def _parse_qc_sample_report(cnf, json_fpath):
    return Record.load_records(json_fpath)

