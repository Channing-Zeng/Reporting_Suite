from ext_modules.simplejson import load
from source.logger import info
from source.reporting import FullReport, SampleReport


def make_summary_reports(cnf, threads, output_dir, callers, samples, jsons_by_sample_by_caller,
                         htmls_by_sample_by_caller):
    if len(jsons_by_sample_by_caller) == 1:
        report = _full_report_for_caller(cnf, samples, jsons_by_sample_by_caller.values()[0],
            htmls_by_sample_by_caller.values()[0])

        full_summary_fpaths = report.save_into_files(
            output_dir, base_fname=cnf.name, caption='Variant QC')

        info()
        info('*' * 70)
        for fpath in full_summary_fpaths:
            if fpath:
                info(fpath)

    else:
        for caller in callers:
            caller.summary_qc_report = _full_report_for_caller(cnf, samples, output_dir,
                jsons_by_sample_by_caller[caller.name], htmls_by_sample_by_caller[caller.name])

            caller.summary_qc_report_fpaths = caller.summary_qc_report.save_into_files(
                output_dir, base_fname=caller.suf + '.' + cnf.name,
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
            output_dir, base_fname=cnf.name, caption='Variant QC')

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


def _full_report_for_caller(cnf, samples, output_dir, jsons_by_sample, htmls_by_sample):
    return FullReport.construct_from_sample_report_jsons(
        samples, None, output_dir, jsons_by_sample, htmls_by_sample)


def _load_sample_report(json_fpath, sample, bcbio_structure):
    with open(json_fpath) as f:
        json = load(f)
    return SampleReport.load(json, sample, bcbio_structure)

