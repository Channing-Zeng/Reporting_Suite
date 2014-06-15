import shutil
import sys
from os import mkdir
from os.path import basename, join

from source.utils_from_bcbio import file_exists
from source.logger import info, step_greetings
from source.utils import call_subprocess, get_tool_cmdline, bgzip_and_tabix_vcf


def bcftools_qc(cnf, vcf_fpath):
    step_greetings('Quality control basic plots')

    bcftools = get_tool_cmdline(cnf, 'bcftools')
    plot_vcfstats = get_tool_cmdline(cnf, 'plot_vcfstats')
    if not bcftools or not bcftools or not plot_vcfstats:
        sys.exit(1)

    gzipped_fpath, tbi_fpath = bgzip_and_tabix_vcf(cnf, vcf_fpath)

    text_report_fpath = join(cnf.work_dir, cnf.name + '_bcftools.report')
    cmdline = '{bcftools} stats {gzipped_fpath}'.format(**locals())
    call_subprocess(cnf, cmdline, None, text_report_fpath)

    viz_report_dir = join(cnf.output_dir, cnf.name + '_qc_plots/')
    if file_exists(viz_report_dir):
        shutil.rmtree(viz_report_dir)
    mkdir(viz_report_dir)
    cmdline = '{plot_vcfstats} -s {text_report_fpath} -p {viz_report_dir} ' \
              '--no-PDF'.format(**locals())
    call_subprocess(cnf, cmdline, text_report_fpath, None, output_is_dir=False)

    results = _get_plots_from_bcftools(cnf, viz_report_dir)
    info('Done.')
    return results


def _get_plots_from_bcftools(cnf, bcftools_report_dir):
    original_plots_names = ['indels.0.png', 'substitutions.0.png']
    final_plots_names = [cnf['name'] + '_indels.png', cnf['name'] + '_substitution.png']
    for i, original_plot in enumerate(original_plots_names):
        plot_src_filename = join(bcftools_report_dir, original_plot)
        plot_dst_filename = join(cnf.output_dir, final_plots_names[i])
        if file_exists(plot_src_filename):
            shutil.copyfile(plot_src_filename, plot_dst_filename)
    return [join(cnf.output_dir, fname) for fname in final_plots_names]

