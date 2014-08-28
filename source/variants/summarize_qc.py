from collections import defaultdict, OrderedDict
from os.path import join
from source.bcbio_structure import Sample
from source.logger import info
from source.reporting import write_summary_report, Record, FullReport, SampleReport


def make_summary_reports(cnf, bcbio_structure):
    callers = bcbio_structure.variant_callers.values()

    if len(callers) == 1:
        report = _full_report_for_caller(cnf, callers[0])

        full_summary_fpaths = write_summary_report(
            cnf.output_dir, cnf.work_dir, report,
            base_fname=cnf.name, caption='Variant QC')

        info()
        info('*' * 70)
        for fpath in full_summary_fpaths:
            info(fpath)

    else:
        for caller in callers:
            caller.summary_qc_report = _full_report_for_caller(cnf, caller)

            caller.summary_qc_report_fpaths = write_summary_report(
                cnf.output_dir, cnf.work_dir, caller.summary_qc_report,
                base_fname=caller.suf + '.' + cnf.name, caption='Variant QC for ' + caller.name)

        # Combining
        combined_full_report = FullReport('', [
            s_report.set_display_name(s_report.sample.name + ' ' + c_name)
            for (_, c_name, s_report) in sorted(
                (s_report.sample.key_to_sort(), c.name, s_report)
                 for c in callers
                 for s_report in c.summary_qc_report.sample_reports)
        ])

        full_summary_fpaths = write_summary_report(
            cnf.output_dir, cnf.work_dir, combined_full_report,
            base_fname=cnf.name, caption='Variant QC')

        info()
        info('*' * 70)

        for caller in callers:
            info(caller.name)
            for fpath in caller.summary_qc_report_fpaths:
                info('  ' + fpath)
            info()

        info('Total')
        for fpath in full_summary_fpaths:
            info('  ' + fpath)


def _full_report_for_caller(cnf, caller):
    jsons_by_sample = caller.get_fpaths_by_sample(cnf.dir, cnf.name, 'json')
    htmls_by_sample = caller.get_fpaths_by_sample(cnf.dir, cnf.name, 'html')

    return FullReport('', [
        SampleReport(sample,
                     records=_parse_qc_sample_report(jsons_by_sample[sample]),
                     html_fpath=htmls_by_sample[sample])
            for sample in caller.samples
            if sample in jsons_by_sample and sample in htmls_by_sample])


def _parse_qc_sample_report(json_fpath):
    return Record.load_records(json_fpath)

