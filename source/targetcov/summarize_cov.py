from os.path import join
from source.bcbio_structure import BCBioStructure
from source.reporting import summarize, write_summary_reports, Record
from source.logger import step_greetings, info


def summary_reports(cnf, bcbio_structure):
    step_greetings('Coverage statistics for all samples')

    json_by_sample = bcbio_structure.get_targetcov_json_by_sample()
    sum_report = summarize(json_by_sample, _parse_targetseq_sample_report)

    final_summary_report_fpaths = write_summary_reports(
        cnf.output_dir,
        cnf.work_dir,
        sum_report,
        join(cnf.output_dir, BCBioStructure.targetseq_name),
        'Target coverage statistics')

    info()
    info('*' * 70)
    info('Summary:')
    for fpath in final_summary_report_fpaths:
        if fpath:
            info('  ' + fpath)


def _parse_targetseq_sample_report(json_fpath):
    return Record.load_records(json_fpath)


