from os.path import join
from source.logger import step_greetings, info, send_email
import source

from source.fastqc.html_template_fastqc import write_fastqc_combo_report


def summary_reports(cnf, bcbio_structure):
    step_greetings('FastQC summary for all samples')

    final_summary_report_fpath = join(cnf.output_dir, source.fastqc_name + '.html')

    write_fastqc_combo_report(final_summary_report_fpath, bcbio_structure.samples)

    info()
    info('*' * 70)
    info('Fastqc summary:')
    info('  ' + final_summary_report_fpath)
    # send_email('Fastqc summary: ' + final_summary_report_fpath)

    return final_summary_report_fpath


