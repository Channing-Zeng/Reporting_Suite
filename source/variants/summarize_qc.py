from ext_modules.simplejson import load
import source
from source.logger import info
from source.reporting import FullReport, SampleReport


def make_summary_reports(cnf, threads, output_dir, callers, samples,
        jsons_by_sample_by_caller, htmls_by_sample_by_caller, tag_by_sample=None,
        varqc_name=source.varqc_name, caption='Variant QC'):

    if len(jsons_by_sample_by_caller) == 1:
        report = _full_report_for_caller(cnf, samples, output_dir,
            jsons_by_sample_by_caller.values()[0], htmls_by_sample_by_caller.values()[0])

        full_summary_fpaths = report.save_into_files(
            output_dir, base_fname=varqc_name, caption=caption)

        info()
        info('*' * 70)
        for fpath in full_summary_fpaths:
            if fpath:
                info(fpath)

    else:
        for caller in callers:
            caller.summary_qc_report = _full_report_for_caller(cnf, samples, output_dir,
                jsons_by_sample_by_caller[caller.name], htmls_by_sample_by_caller[caller.name])

            if tag_by_sample:
                for s_report in caller.summary_qc_report.sample_reports:
                    s_report.set_project_tag(tag_by_sample[s_report.sample.name])

            caller.summary_qc_report_fpaths = caller.summary_qc_report.save_into_files(
                output_dir, base_fname=caller.suf + '.' + varqc_name,
                caption=caption + ' for ' + caller.name)

        # Combining
        sample_reports = []
        for caller in callers:
            for s_report in caller.summary_qc_report.sample_reports:
                s_report.set_caller_tag(caller.name)
                sample_reports.append(s_report)

        sample_reports = [sr for _, _, sr in sorted((sr.sample.key_to_sort(), sr.caller_tag, sr) for sr in sample_reports)]

        combined_full_report = FullReport('', sample_reports)
        full_summary_fpaths = combined_full_report.save_into_files(
            output_dir, base_fname=varqc_name, caption=caption)

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

