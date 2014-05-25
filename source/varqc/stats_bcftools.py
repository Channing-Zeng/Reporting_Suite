import shutil
import textwrap
from os import mkdir, makedirs
from os.path import basename, join, isdir, dirname, expanduser

from source.bcbio_utils import file_exists
from source.utils import info, err, verify_file, verify_dir, critical, call, get_tool_cmdline, step_greetings


def bcftools_qc(cnf, qc_dir, vcf_fpath):
    step_greetings(cnf, 'Quality control basic plots')

    work_dir = cnf['work_dir']

    bgzip = get_tool_cmdline(cnf, 'bgzisp')
    tabix = get_tool_cmdline(cnf, 'tabix')
    bcftools = get_tool_cmdline(cnf, 'bcftools')
    plot_vcfstats = get_tool_cmdline(cnf, 'plot_vcfstats')
    if not bcftools or not tabix or not bcftools or not plot_vcfstats:
        exit()

    gzipped_fpath = join(work_dir, basename(vcf_fpath) + '.gz')
    cmdline = '{bgzip} -c {vcf_fpath}'.format(**locals())
    call(cnf, cmdline, None, gzipped_fpath)

    tbi_fpath = gzipped_fpath + '.tbi'
    cmdline = '{tabix} -f -p vcf {gzipped_fpath}'.format(**locals())
    call(cnf, cmdline, None, tbi_fpath)

    text_report_fpath = join(work_dir, cnf['name'] + '_bcftools.report')
    cmdline = '{bcftools} stats {gzipped_fpath}'.format(**locals())
    call(cnf, cmdline, None, text_report_fpath, to_remove=[gzipped_fpath, tbi_fpath])

    viz_report_dir = join(work_dir, cnf['name'] + '_qc_plots/')
    if file_exists(viz_report_dir):
        shutil.rmtree(viz_report_dir)
    mkdir(viz_report_dir)
    cmdline = '{plot_vcfstats} -s {text_report_fpath} -p {viz_report_dir} ' \
              '--no-PDF'.format(**locals())
    call(cnf, cmdline, text_report_fpath, None, output_is_file=False)
    return _get_plots_from_bcftools(cnf, viz_report_dir, qc_dir)


def _get_plots_from_bcftools(cnf, bcftools_report_dir, output_dir):
    original_plots_names = ['indels.0.png', 'substitutions.0.png']
    final_plots_names = [cnf['name'] + '_indels.png', cnf['name'] + '_substitution.png']
    for i, original_plot in enumerate(original_plots_names):
        plot_src_filename = join(bcftools_report_dir, original_plot)
        plot_dst_filename = join(output_dir, final_plots_names[i])
        if file_exists(plot_src_filename):
            shutil.copyfile(plot_src_filename, plot_dst_filename)
    return [join(output_dir, fname) for fname in final_plots_names]

