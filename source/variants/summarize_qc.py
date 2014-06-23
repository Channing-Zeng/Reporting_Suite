from collections import OrderedDict
from os.path import join, basename

from source.reporting import write_summary_reports, summarize
from source.utils import verify_file
from source.logger import info, critical


database = 'cosmic'
novelty = 'all'
metrics_header = 'Metric'
novelty_header = 'Novelty'
sample_header = 'Sample name:'


def write_qc_summart_reports(bcbio_final_dir, samples_fpath, output_dirpath,
                             report_basedir, vcf_suf=''):
    sample_report_suffix = '.varqc.txt'
    sample_report_fpaths = []
    sample_names = []

    with open(samples_fpath) as f:
        for line in f:
            sample_name = line.strip()
            sample_names.append(sample_name)

            sample_report_fpath = join(
                bcbio_final_dir, sample_name, report_basedir,
                sample_name + '-' + vcf_suf + sample_report_suffix)

            info(basename(sample_report_fpath))

            if not verify_file(sample_report_fpath):
                critical(sample_report_fpath + ' does not exist.')

            sample_report_fpaths.append(sample_report_fpath)

    report = summarize(sample_names, sample_report_fpaths, _parse_qc_report)

    return write_summary_reports(
        report, sample_names, output_dirpath, vcf_suf + '.varqc.summary')


def _parse_qc_report(report_fpath):
    pairs = []

    with open(report_fpath, 'r') as f:
        # parsing Sample name and Database columns
        database_col_id = None
        novelty_col_id = None
        for line in f:
            if line.startswith(sample_header):
                sample_name = line[len(sample_header):].strip()
            elif line.startswith(metrics_header):
                if database in line:
                    database_col_id = line.split().index(database)
                if novelty_header in line:
                    novelty_col_id = line.split().index(novelty_header)
                break

        if database_col_id:
            # parsing rest of the report
            for line in f:
                if novelty_col_id and line.split()[novelty_col_id] != novelty:
                    continue
                cur_metric_name = line.split()[0]
                cur_value = line.split()[database_col_id]
                pairs.append((cur_metric_name, cur_value))

    return pairs