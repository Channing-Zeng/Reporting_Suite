from source.reporting import parse_tsv
from source.targetcov.copy_number import run_copy_number
from source.logger import critical


def parse_targetseq_sample_report(report_fpath):
    """ returns row_per_sample =
            dict(metricName=None, value=None,
            isMain=True, quality='More is better')
    """
    row_per_sample = []

    with open(report_fpath) as f:
        rows = [l.split('\t') for l in f]

    for row in rows:
        row_per_sample.append(dict(
            metricName=row[0], value=row[1],
            isMain=True, quality='More is better'))

    return row_per_sample



def summarize_copy_number(sample_names, report_details_fpaths, report_summary_fpaths):
    """ Parsing gene coverage and sample summary report as an input to copy number report
        "Gene-Amplicon" row's used from gene coverage and "Mapped reads" form summary
    """
    gene_summary_lines = []
    cov_by_sample = dict()

    for sample_name, report_details_fpath, report_summary_fpath in \
            zip(sample_names, report_details_fpaths, report_summary_fpaths):

        gene_summary_lines += _get_lines_by_region_type(report_details_fpath, 'Gene-Amplicon')
        report_lines = dict(parse_tsv(report_summary_fpath))
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
