import shutil
import textwrap
from os import mkdir, makedirs
from os.path import basename, join, isdir, dirname, expanduser

from source.bcbio_utils import file_exists
from source.utils import info, err, verify_file, step_greetings, \
    get_tool_cmdline, get_java_tool_cmdline, call, verify_dir, \
    critical, human_sorted

import matplotlib
matplotlib.use('Agg')  # non-GUI backen
import matplotlib.pyplot


def check_quality_control_config(cnf):
    qc_cnf = cnf.get('quality_control')
    if not qc_cnf:
        return

    if 'databases' not in qc_cnf:
        qc_cnf['databases'] = ['dbsnp']
        info('Warning: not databases for quality control, using [dbsnp]')

    if 'novelty' not in qc_cnf:
        qc_cnf['novelty'] = ['all', 'known', 'novel']
        info('Warning: no novelty specified for quality control, '
             'using default ' + ', '.join(qc_cnf['novelty']))

    if 'metrics' not in qc_cnf:
        qc_cnf['metircs'] = [
           'nEvalVariants', 'nSNPs', 'nInsertions', 'nDeletions',
           'nVariantsAtComp', 'compRate', 'nConcordant', 'concordantRate',
           'variantRate', 'variantRatePerBp', 'hetHomRatio', 'tiTvRatio']
        info('Warning: no metrics for quality control, using '
             'default ' + ', '.join(qc_cnf['metircs']))

    if 'variants_distribution_scale' not in qc_cnf:
        qc_cnf['variants_distribution_scale'] = 1000
        info('Warning: no variants distribution scale specified for quality control, '
             'using default ' + str(qc_cnf['variants_distribution_scale']))

    to_exit = False
    dbs_dict = {}
    for db in qc_cnf['databases']:
        if not db:
            err('Empty field for quality_control databases')
            to_exit = True
        elif file_exists(db):
            if not verify_file(db, 'Vcf'):
                to_exit = True
            dbs_dict[basename(db)] = db
        elif db not in cnf['genome']:
            to_exit = True
            err(cnf.get('log'), db + ' for variant qc is not found '
                                     'in genome resources in system config.')
        else:
            dbs_dict[db] = cnf['genome'][db]

    if to_exit:
        exit()

    qc_cnf['database_vcfs'] = dbs_dict

    if 'summary_output' in qc_cnf or 'qc_summary_output' in cnf:
        qc_output_fpath = qc_cnf.get('summary_output') or cnf.get('qc_summary_output')
        summary_output_dir = dirname(qc_output_fpath)
        if not isdir(summary_output_dir):
            try:
                makedirs(summary_output_dir)
            except OSError:
                critical('ERROR: cannot create directory for '
                         'qc summary report: ' + summary_output_dir)
        if not verify_dir(summary_output_dir, 'qc_summary_output'):
            exit()


def quality_control(cnf, qc_dir, vcf_fpath):
    if 'quality_control' not in cnf:
        return None, None

    if not isdir(qc_dir):
        mkdir(qc_dir)

    qc_report_fpath = gatk_qc(cnf, qc_dir, vcf_fpath)
    qc_plots_fpaths = bcftools_qc(cnf, qc_dir, vcf_fpath)
    qc_var_distr_plot_fpath = variants_distribution_plot(cnf, qc_dir, vcf_fpath)
    return qc_report_fpath, [qc_var_distr_plot_fpath] + qc_plots_fpaths


def gatk_qc(cnf, qc_dir, vcf_fpath):
    step_greetings(cnf, 'Quality control reports')

    log = cnf['log']
    work_dir = cnf['work_dir']

    qc_cnf = cnf['quality_control']
    databases = qc_cnf.get('database_vcfs')
    novelty = qc_cnf.get('novelty')
    metrics = qc_cnf.get('metrics')

    executable = get_java_tool_cmdline(cnf, 'gatk')
    gatk_opts_line = ' '.join(cnf.get('gatk', {'options': []}).get('options', []))
    ref_fpath = cnf['genome']['seq']
    report_fpath = join(work_dir, cnf['name'] + '_gatk.report')

    cmdline = ('{executable} {gatk_opts_line} -R {ref_fpath} -T VariantEval'
               ' --eval:tmp {vcf_fpath} -o {report_fpath}').format(**locals())

    if 'dbsnp' in databases:
        cmdline += ' -D ' + databases['dbsnp']
    for db_name, db_path in databases.items():
        if not db_name == 'dbsnp':
            cmdline += ' -comp:' + db_name + ' ' + db_path

    call(cnf, cmdline, None, report_fpath, stdout_to_outputfile=False,
         to_remove=[vcf_fpath + '.idx'])

    report = _parse_gatk_report(report_fpath, databases.keys(), novelty, metrics)

    final_report_fpath = join(qc_dir, cnf['name'] + '_qc.report')

    _make_final_report(report, final_report_fpath, cnf['name'],
                       databases.keys(), novelty, metrics)
    return final_report_fpath


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


def variants_distribution_plot(cnf, qc_dir, vcf_fpath):
    step_greetings(cnf, 'Quality control variant distribution plots')

    # step 1: get chr lengths
    chr_len_fpath = cnf.get('chr_lengths')
    if not chr_len_fpath:
        chr_len_fpath = expanduser(cnf['genome'].get('chr_lengths'))
    if not verify_file(chr_len_fpath):
        exit(1)
        # no chromosome lengths file for the genome! TODO: process reference fasta and get lengths from it
    chr_lengths = _get_chr_lengths(chr_len_fpath)

    # step 2: get variants distribution (per chromosome)
    qc_cnf = cnf['quality_control']
    variants_per_kbp = qc_cnf.get('variants_distribution_scale')
    plot_scale = 1000 * variants_per_kbp
    variants_distribution, not_counted = _get_variants_distribution(vcf_fpath, chr_lengths, plot_scale)
    if not_counted:
        info(cnf['log'], 'Warning: some variants were not counted (chromosome names not found): ' + str(not_counted))
    empty_chr = []
    for chr_name in chr_lengths.keys():
        if sum(variants_distribution[chr_name]) == 0:
            empty_chr.append(chr_name)
            del variants_distribution[chr_name]
    if empty_chr:
        info(cnf['log'], 'Chromosomes without variants: ' + ', '.join(human_sorted(empty_chr)))

    # step 3: plotting
    nplots = len(variants_distribution.keys())
    ncols = min(4, nplots)
    nrows = 1 + (nplots - 1) / 4
    fontsize = 15
    mbp = 1000000
    fig = matplotlib.pyplot.figure(figsize=(6 * ncols, 3 * nrows))
    for id, chr_name in enumerate(human_sorted(variants_distribution.keys())):
        #fig = matplotlib.pyplot.figure(figsize=(25,6))

        ax = fig.add_subplot(nrows, ncols, id + 1)
        ax.plot(range(len(variants_distribution[chr_name])), variants_distribution[chr_name], color='#46a246', linewidth=3.0)
        ax.axhline(linewidth=3.0)
        ax.axvline(linewidth=3.0)

        #Axe style###############################
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ax.set_xticklabels('')
        ax.set_ylim(bottom=0)
        if chr_lengths[chr_name] < mbp / 10:
            chr_size = '<0.1 Mbp'
        else:
            chr_size = '%.1f Mbp' % (float(chr_lengths[chr_name]) / mbp)
        ax.set_xlabel(chr_name + ', ' + chr_size)

        for item in ([ax.xaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(fontsize)
            item.set_weight('bold')

        #TODO: individual plots for each chromosome separately
        #fig.savefig(join(qc_dir, cnf['name'] + '_qc_' + chr_name + '.png'))
        #matplotlib.pyplot.close(fig)

    fig.suptitle('Variants per %d kbp' % variants_per_kbp, y=0.93, fontsize=int(1.5 * fontsize), fontweight='bold')

    if empty_chr:
        if len(empty_chr) > 1:
            chr_without_variants = ', '.join(human_sorted(empty_chr))
            renderer = fig.canvas.get_renderer()
            aspect_ratio = 0.5 # This varies with the font!!
            pixels_per_char = aspect_ratio * renderer.points_to_pixels(fontsize)
            new_width = int(fig.dpi * fig.get_figwidth() * 0.7)
            wrap_width = max(1, new_width // pixels_per_char)
            chr_without_variants_wrapped = textwrap.fill(chr_without_variants, wrap_width)
        else:
            chr_without_variants_wrapped = empty_chr[0]
        fig.suptitle('Number of variants from chromosomes not found in reference: ' + str(not_counted),
                     x=0.12, y=0.07, fontsize=fontsize, fontweight='bold', ha='left')
        fig.suptitle('Chromosomes without variants:', x=0.12, y=0.05, fontsize=fontsize, fontweight='bold', ha='left')
        fig.suptitle(chr_without_variants_wrapped, x=0.12, y=0.03, fontsize=fontsize, ha='left')

    variants_distribution_plot_fpath = join(qc_dir, cnf['name'] + '_qc_variant_distribution.png')
    fig.savefig(variants_distribution_plot_fpath, bbox_inches='tight')
    matplotlib.pyplot.close(fig)
    return variants_distribution_plot_fpath


def _parse_gatk_report(report_filename, databases, novelty, metrics):
    database_col_name = 'CompRod'
    database_col_id = None
    novelty_col_name = 'Novelty'
    novelty_col_id = None

    report = dict()
    comments_section = False
    cur_header = []
    cur_metrics_ids = []
    for line in open(report_filename):
        if not line.strip():
            continue
        if line.startswith('#'): # comment line
            comments_section = True
            continue
        elif comments_section:
            comments_section = False
            cur_header = line.split()
            cur_metrics_ids = []
            database_col_id = cur_header.index(database_col_name)
            novelty_col_id = cur_header.index(novelty_col_name)
            for metric in metrics:
                if metric in cur_header:
                    cur_metrics_ids.append(cur_header.index(metric))
                    if metric not in report:
                        report[metric] = dict()
        elif cur_metrics_ids:  # process lines only if there are metrics in current section
            values = line.split()
            cur_database = values[database_col_id]
            cur_novelty = values[novelty_col_id]
            if (cur_database not in databases) or (cur_novelty not in novelty):
                continue
            for metric_id in cur_metrics_ids:
                if cur_database not in report[cur_header[metric_id]]:
                    report[cur_header[metric_id]][cur_database] = dict()
                report[cur_header[metric_id]][cur_database][cur_novelty] = values[metric_id]
    return report


def _make_final_report(report_dict, report_filename, sample_name,
                       databases, novelty, metrics):
    header = ['Metric', 'Novelty'] + databases + ['Average']
    full_report = [header]
    for cur_metric in metrics:
        for cur_novelty in novelty:
            cur_row = [cur_metric, cur_novelty]
            sum = 0.0
            for cur_database in databases:
                if cur_metric == 'variantRatePerBp': # confusing name and value format
                    cur_row[0] = 'basesPerVariant'
                    cur_row.append("%.2f" % float(report_dict[cur_metric][cur_database][cur_novelty]))
                else:
                    cur_row.append(report_dict[cur_metric][cur_database][cur_novelty])
                sum += float(cur_row[-1])
            average = sum / len(databases)
            cur_row.append("%.2f" % average)
            full_report.append(cur_row)

    col_widths = [0] * len(header)
    for row in full_report:
        for id, value in enumerate(row):
            col_widths[id] = max(len(value), col_widths[id])

    out = open(report_filename, 'w')
    out.write('Sample name: ' + sample_name + '\n\n')
    for row in full_report:
        out.write('  '.join('%-*s' % (col_width, value) for col_width, value
                            in zip(col_widths, row)) + "\r\n")
    out.close()


def _get_plots_from_bcftools(cnf, bcftools_report_dir, output_dir):
    original_plots_names = ['indels.0.png', 'substitutions.0.png']
    final_plots_names = [cnf['name'] + '_indels.png', cnf['name'] + '_substitution.png']
    for i, original_plot in enumerate(original_plots_names):
        plot_src_filename = join(bcftools_report_dir, original_plot)
        plot_dst_filename = join(output_dir, final_plots_names[i])
        if file_exists(plot_src_filename):
            shutil.copyfile(plot_src_filename, plot_dst_filename)
    return [join(output_dir, fname) for fname in final_plots_names]


def _get_chr_lengths(chr_len_fpath):
    chr_lengths = dict()
    with open(chr_len_fpath, 'r') as f:
        for line in f:
            if len(line.split()) == 2:
                chr_name = line.split()[0]
                chr_length = int(line.split()[1])
                chr_lengths[chr_name] = chr_length
    return chr_lengths


def _get_variants_distribution(vcf_fpath, chr_lengths, plot_scale):
    variants_distribution = dict()
    for chr_name, chr_length in chr_lengths.items():
        variants_distribution[chr_name] = [0] * max(1, chr_length / plot_scale)
    not_counted = 0

    with open(vcf_fpath, 'r') as f:
        for line in f:
            if line.startswith('#') or len(line.split()) < 8:
                continue
            chr_name = line.split()[0]
            chr_pos = int(line.split()[1])
            if chr_name not in variants_distribution:
                not_counted += 1
            else:
                region_id = min((chr_pos - 1) / plot_scale, len(variants_distribution[chr_name]) - 1)
                variants_distribution[chr_name][region_id] += 1

    # the last region in each chromosome is not exactly equal to plot_scale
    for chr_name, chr_length in chr_lengths.items():
        last_region_length = chr_length % plot_scale + (0 if chr_length < plot_scale else plot_scale)
        variants_distribution[chr_name][-1] = int(variants_distribution[chr_name][-1] * plot_scale /
                                                  float(last_region_length))
    return variants_distribution, not_counted

