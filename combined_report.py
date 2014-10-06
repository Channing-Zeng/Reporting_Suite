#!/usr/bin/env python

import sys

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from os.path import abspath, dirname, realpath, join
from site import addsitedir
from collections import OrderedDict
source_dir = abspath(dirname(realpath(__file__)))
addsitedir(join(source_dir, 'ext_modules'))

from source.prepare_args_and_cnf import summary_script_proc_params
from source.bcbio_structure import BCBioStructure
from source.logger import info, step_greetings
from source.file_utils import verify_file
from source.reporting import Metric, Record, MetricStorage, ReportSection, SampleReport, FullReport


def main():
    info(' '.join(sys.argv))
    info()

    cnf, bcbio_structure = summary_script_proc_params(BCBioStructure.combined_report_name,
                                                      BCBioStructure.combined_report_dir)

    make_combined_report(cnf, bcbio_structure)


def make_combined_report(cnf, bcbio_structure):
    step_greetings('Combined report')

    # summary reports
    general_section = ReportSection('general_section', '', [])
    general_records = []
    for (step_name, step_summary_dir) in [(bcbio_structure.targetseq_name, bcbio_structure.targetseq_summary_dir),
                                          (bcbio_structure.varqc_name, bcbio_structure.varqc_summary_dir),
                                          (bcbio_structure.ngscat_name, bcbio_structure.ngscat_summary_dir),
                                          (bcbio_structure.qualimap_name, bcbio_structure.qualimap_summary_dir)]:
        summary_report_fpath = join(bcbio_structure.date_dirpath, step_summary_dir, step_name + '.html')
        if verify_file(summary_report_fpath):
            cur_metric = Metric(step_name + ' summary', common=True)
            general_section.add_metric(cur_metric)
            general_records.append(Record(metric=cur_metric,
                                          value=cur_metric.name,
                                          html_fpath=summary_report_fpath))

    # individual reports
    individual_reports_section = ReportSection('individual_reports', '', [])

    # varQC reports -- special case
    callers = bcbio_structure.variant_callers.values()
    if len(callers) == 0:
        vq_htmls_by_sample = None
    elif len(callers) == 1:
        vq_htmls_by_sample = callers[0].find_fpaths_by_sample(bcbio_structure.varqc_dir,
                                                              bcbio_structure.varqc_name, 'html')
    else:
        vq_htmls_by_sample = OrderedDict()
        for sample in bcbio_structure.samples:
            vq_htmls_by_sample[sample.name] = OrderedDict()
        for caller in callers:
            for sample, fpath in caller.find_fpaths_by_sample(bcbio_structure.varqc_dir,
                                                              bcbio_structure.varqc_name, 'html').items():
                vq_htmls_by_sample[sample][caller.name] = fpath
    # other reports
    tc_htmls_by_sample = bcbio_structure.get_targetcov_report_fpaths_by_sample('html')
    nc_htmls_by_sample = bcbio_structure.get_ngscat_report_fpaths_by_sample()
    qm_htmls_by_sample = bcbio_structure.get_qualimap_report_fpaths_by_sample()

    sample_reports_records = dict()
    for sample in bcbio_structure.samples:
        sample_reports_records[sample.name] = list(general_records)

    for (step_name, htmls_by_sample) in [(bcbio_structure.targetseq_name, tc_htmls_by_sample),
                                         (bcbio_structure.varqc_name, vq_htmls_by_sample),
                                         (bcbio_structure.ngscat_name, nc_htmls_by_sample),
                                         (bcbio_structure.qualimap_name, qm_htmls_by_sample)]:
        if htmls_by_sample:
            cur_metric = Metric(step_name)
            individual_reports_section.add_metric(cur_metric)
            for sample in bcbio_structure.samples:
                if sample.name in htmls_by_sample:
                    sample_reports_records[sample.name].append(Record(metric=cur_metric,
                                                                      value=cur_metric.name,
                                                                      html_fpath=htmls_by_sample[sample.name]))
                else:
                    sample_reports_records[sample.name].append(Record(metric=cur_metric,
                                                                      value=None,
                                                                      html_fpath=None))

    metric_storage = MetricStorage(general_section=general_section, sections=[individual_reports_section])
    sample_reports = []
    for sample in bcbio_structure.samples:
        sample_reports.append(SampleReport(sample,
                                           records=sample_reports_records[sample.name],
                                           html_fpath=None,
                                           metric_storage=metric_storage))
    full_report = FullReport(cnf.name, sample_reports, metric_storage=metric_storage)
    final_summary_report_fpath = full_report.save_html(
        cnf.output_dir, bcbio_structure.combined_report_name, 'Combined report for ' +
                                                              bcbio_structure.project_name)

    info()
    info('*' * 70)
    info('Combined report saved in: ')
    info('  ' + final_summary_report_fpath)


if __name__ == '__main__':
    main()