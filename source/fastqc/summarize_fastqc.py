from os.path import join
from source.logger import step_greetings, info
from source.bcbio_structure import BCBioStructure

from source.fastqc.html_template_fastqc import print_html


def summary_reports(cnf, bcbio_structure):
    step_greetings('FastQC summary for all samples')

    htmls_by_sample = bcbio_structure.get_fastqc_report_fpaths_by_sample()

    if not htmls_by_sample:
        return None

    final_summary_report_fpath = join(cnf.output_dir, BCBioStructure.fastqc_name + '.html')

    info()
    info('*' * 70)
    info('Summary:')
    info('  ' + final_summary_report_fpath )
    print_html(final_summary_report_fpath, htmls_by_sample)
    return final_summary_report_fpath


