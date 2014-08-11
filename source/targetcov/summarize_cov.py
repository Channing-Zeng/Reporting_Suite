from source.reporting import parse_tsv, get_per_sample_fpaths_for_bcbio_final_dir, \
    summarize, write_summary_reports, write_tsv_rows, Metric, Record
from source.targetcov.copy_number import run_copy_number
from source.logger import critical, step_greetings, info, err
from source.targetcov.cov import detail_gene_report_ending, cov_json_ending


def summary_reports(cnf, sample_names):
    step_greetings('Coverage statistics for all samples')

    sample_json_fpaths, sample_names = get_per_sample_fpaths_for_bcbio_final_dir(
        cnf['bcbio_final_dir'], sample_names, cnf['base_name'], cov_json_ending)

    sum_report = summarize(sample_names, sample_json_fpaths, _parse_targetseq_sample_report)

    final_summary_report_fpaths = write_summary_reports(
        cnf['output_dir'], cnf['work_dir'], sum_report, 'targetSeq', 'Target coverage statistics')

    return sample_json_fpaths, final_summary_report_fpaths


def cnv_reports(cnf, sample_names, sample_json_reports):
    step_greetings('Coverage statistics for each gene for all samples')

    info('Collecting sample reports...')
    sample_gene_reports, sample_names = get_per_sample_fpaths_for_bcbio_final_dir(
        cnf['bcbio_final_dir'], sample_names, cnf['base_name'], detail_gene_report_ending)

    if not sample_gene_reports:
        err('No gene reports, cannot call copy numbers.')
        return None

    info('Calculating normalized coverages for CNV...')
    cnv_rows = _summarize_copy_number(sample_names, sample_gene_reports, sample_json_reports)

    cnv_report_fpath = write_tsv_rows(cnv_rows, cnf['output_dir'], 'Seq2C')

    return cnv_report_fpath


def _parse_targetseq_sample_report(json_fpath):
    return Record.load_records(json_fpath)


# def save_results_separate_for_samples(results):
#     header = results[0]
#
#     results_per_sample = OrderedDict()
#
#     for fields in results[1:]:
#         sample_name = fields[0]
#         if sample_name not in results_per_sample:
#             results_per_sample[sample_name] = [header]
#
#         results_per_sample[sample_name].append(fields)
#
#     for sample_name, fields in results_per_sample.items():
#


def _summarize_copy_number(sample_names, report_details_fpaths, sample_json_fpaths):
    gene_summary_lines = []
    cov_by_sample = dict()

    for sample_name, report_details_fpath, report_json_fpath in \
            zip(sample_names, report_details_fpaths, sample_json_fpaths):

        gene_summary_lines += _get_lines_by_region_type(report_details_fpath, 'Gene-Amplicon')

        records = Record.load_records(report_json_fpath)

        cov_by_sample[sample_name] = int(next(rec.value for rec in records if rec.metric.name == 'Mapped reads'))

    results = run_copy_number(cov_by_sample, gene_summary_lines)

    # save_results_separate_for_samples(results)

    return results


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
