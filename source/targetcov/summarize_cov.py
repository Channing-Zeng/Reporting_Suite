from os.path import join
from source.bcbio_structure import BCBioStructure
from source.reporting import summarize, write_summary_reports, write_tsv_rows, Record
from source.targetcov.copy_number import run_copy_number
from source.logger import critical, step_greetings, info, err


def summary_reports(cnf, bcbio_structure):
    step_greetings('Coverage statistics for all samples')

    json_by_sample = bcbio_structure.get_targetcov_json_by_sample()
    sum_report = summarize(json_by_sample, _parse_targetseq_sample_report)

    final_summary_report_fpaths = write_summary_reports(
        cnf.output_dir,
        cnf.work_dir,
        sum_report,
        join(bcbio_structure.date_dirpath, BCBioStructure.targetseq_summary_dir),
        'Target coverage statistics')

    return json_by_sample.values(), final_summary_report_fpaths


def _parse_targetseq_sample_report(json_fpath):
    return Record.load_records(json_fpath)


