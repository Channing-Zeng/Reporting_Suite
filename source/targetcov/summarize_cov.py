from collections import OrderedDict
from source.reporting import parse_tsv, get_sample_report_fpaths_for_bcbio_final_dir, \
    summarize, write_summary_reports, write_tsv, Metric, Record
from source.targetcov.copy_number import run_copy_number
from source.logger import critical, step_greetings, info, err
from source.utils import OrderedDefaultDict


def summary_reports(cnf, sample_names):
    step_greetings('Coverage statistics for all samples')

    sample_sum_reports, sample_names = get_sample_report_fpaths_for_bcbio_final_dir(
        cnf['bcbio_final_dir'], sample_names, cnf['base_name'], '.targetSeq.json')

    sum_report = summarize(sample_names, sample_sum_reports, _parse_targetseq_sample_report)

    sum_report_fpaths = write_summary_reports(
        cnf['output_dir'], cnf['work_dir'], sum_report, 'targetSeq', 'Target coverage statistics')

    return sample_sum_reports, sum_report_fpaths


def cnv_reports(cnf, sample_names, sample_sum_reports):
    step_greetings('Coverage statistics for each gene for all samples')

    info('Collecting sample reports...')
    sample_gene_reports, sample_names = get_sample_report_fpaths_for_bcbio_final_dir(
        cnf['bcbio_final_dir'], sample_names, cnf['base_name'], '.targetseq.details.gene.txt')

    if not sample_gene_reports:
        err('No gene reports, cannot call copy numbers.')
        return None

    info('Calculating normalized coverages for CNV...')
    cnv_rows = _summarize_copy_number(sample_names, sample_gene_reports, sample_sum_reports)

    cnv_report_fpath = write_tsv(cnv_rows, cnf['output_dir'], 'Seq2C')

    return cnv_report_fpath


def _parse_targetseq_sample_report(json_fpath):
    with open(json_fpath) as f:
        return Record.load_records(f)


def _summarize_copy_number(sample_names, report_details_fpaths, report_summary_fpaths):
    gene_summary_lines = []
    cov_by_sample = dict()

    for sample_name, report_details_fpath, report_summary_fpath in \
            zip(sample_names, report_details_fpaths, report_summary_fpaths):

        gene_summary_lines += _get_lines_by_region_type(report_details_fpath, 'Gene-Amplicon')

        report_lines = OrderedDict(parse_tsv(report_summary_fpath))

        cov_by_sample[sample_name] = int(report_lines.get('Mapped reads').replace(',', ''))

    return run_copy_number(cov_by_sample, gene_summary_lines)


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
