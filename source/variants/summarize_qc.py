from ext_modules.simplejson import load
from source.logger import info
from source.reporting import FullReport, SampleReport


def make_summary_reports(cnf, bcbio_structure):
    callers = bcbio_structure.variant_callers.values()

    if len(callers) == 1:
        report = _full_report_for_caller(cnf, callers[0])

        full_summary_fpaths = report.save_into_files(
            cnf.output_dir, base_fname=cnf.name, caption='Variant QC')

        info()
        info('*' * 70)
        for fpath in full_summary_fpaths:
            if fpath:
                info(fpath)

    else:
        for caller in callers:
            caller.summary_qc_report = _full_report_for_caller(cnf, caller)

            caller.summary_qc_report_fpaths = caller.summary_qc_report.save_into_files(
                cnf.output_dir, base_fname=caller.suf + '.' + cnf.name,
                caption='Variant QC for ' + caller.name)

        # Combining
        combined_full_report = FullReport('', [
            s_report.set_display_name(s_report.sample.name + ' ' + c_name)
            for (_, c_name, s_report) in sorted(
                (s_report.sample.key_to_sort(), c.name, s_report)
                 for c in callers
                 for s_report in c.summary_qc_report.sample_reports)
        ])

        full_summary_fpaths = combined_full_report.save_into_files(
            cnf.output_dir, base_fname=cnf.name, caption='Variant QC')

        info()
        info('*' * 70)

        for caller in callers:
            info(caller.name)
            for fpath in caller.summary_qc_report_fpaths:
                if fpath:
                    info('  ' + fpath)
            info()

        info('Total')
        for fpath in full_summary_fpaths:
            if fpath:
                info('  ' + fpath)


def _full_report_for_caller(cnf, caller):
    jsons_by_sample = caller.find_fpaths_by_sample(cnf.dir, cnf.name, 'json')
    htmls_by_sample = caller.find_fpaths_by_sample(cnf.dir, cnf.name, 'html')

    return FullReport.construct_from_sample_report_jsons(
        caller.samples, caller.bcbio_structure, jsons_by_sample,
        htmls_by_sample, cnf.output_dir)


def _load_sample_report(json_fpath, sample, bcbio_structure):
    with open(json_fpath) as f:
        json = load(f)
    return SampleReport.load(json, sample, bcbio_structure)

