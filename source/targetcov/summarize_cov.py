from os.path import join, basename

from source.reporting import parse_tsv, summarize, write_summary_reports
from source.targetcov.copy_number import run_copy_number
from source.logger import info, critical, step_greetings
from source.utils import verify_file


def write_cov_summary_reports(out_dirpath, samples_fpath, report_basedir):
    sample_report_suffix = '.targetseq.summary.txt'
    sample_report_fpaths = []
    sample_names = []

    with open(samples_fpath, 'r') as f:
        for line in f:
            sample_name = line.strip()
            sample_names.append(sample_name)

            sample_report_fpath = join(
                out_dirpath, sample_name,
                report_basedir, sample_name + sample_report_suffix)

            info(basename(sample_report_fpath))

            if not verify_file(sample_report_fpath):
                critical(sample_report_fpath + ' does not exist.')

            sample_report_fpaths.append(sample_report_fpaths)

    report = summarize(sample_names, sample_report_fpaths,
                       _parse_targetseq_report)

    return write_summary_reports(
        report, sample_names, out_dirpath, '.targetseq.summary')


def _parse_targetseq_report(report_fpath):
    with open(report_fpath) as f:
        return [tuple(l.split('\t')) for l in report_fpath]


def write_cov_gene_summary_reports(out_dirpath, samples_fname, report_basedir):
    report_details_suffix = '.targetseq.details.gene.txt'
    report_summary_suffix = '.targetseq.summary.txt'

    copy_number_summary_fpath = join(out_dirpath, 'targetcov_cnv.txt')

    report_fpaths = []
    report_summary_fpaths = []

    with open(samples_fname) as f:
        for line in f:
            sample_name = line.strip()

            report_details_fpath = join(out_dirpath, sample_name, report_basedir, sample_name + report_details_suffix)
            summary_report_fpath = join(out_dirpath, sample_name, report_basedir, sample_name + report_summary_suffix)
            info(sample_name + ': ' + report_details_fpath)
            info(sample_name + ': ' + summary_report_fpath)
            if not verify_file(report_details_fpath):
                critical(report_details_fpath + ' does not exist.')
            if not verify_file(summary_report_fpath):
                critical(summary_report_fpath + ' does not exist.')

            report_fpaths.append(report_details_fpath)
            report_summary_fpaths.append(summary_report_fpath)

    return summarize_copy_number(report_fpaths, report_summary_fpaths,
        report_summary_suffix, copy_number_summary_fpath)


def summarize_copy_number(sample_names, report_details_fpaths,
                          report_summary_fpaths):
    """ Parsing gene coverage and sample summary report as an input to copy number report
        "Gene-Amplicon" row's used from gene coverage and "Mapped reads" form summary
    """
    gene_summary_lines = []
    cov_by_sample = dict()

    for sample_name, report_details_fpath, report_summary_fpath in \
            zip(sample_names, report_details_fpaths, report_summary_fpaths):

        gene_summary_lines += _get_lines_by_region_type(report_details_fpath, 'Gene-Amplicon')
        report_lines = parse_tsv(report_summary_fpath)
        cov_by_sample[sample_name] = int(report_lines.get('Mapped reads').replace(',', ''))

    report_data = run_copy_number(cov_by_sample, gene_summary_lines)

    return report_data


def _get_lines_by_region_type(report_fpath, region_type):
    gene_summary_lines = []

    with open(report_fpath, 'r') as f:
        for line in f:
            if region_type in line:
                gene_summary_lines.append(line.split()[:8])

    if not gene_summary_lines:
        critical('Regions of type ' + region_type +
                 ' not found in ' + gene_summary_lines)

    return gene_summary_lines
