# VariantQC script that takes 2 inputs (the vcf file name and the output directory) and 
# computes quality metrics for the vcf using GATK and bcftools pipeline
import os
import subprocess
import sys
import shutil
import tempfile

SUPPORTED_PYTHON_VERSIONS = ['2.5', '2.6', '2.7']
DEBUG_MODE = False  # whether to save temp files or not

databases = ['dbsnp', 'Cosmic', '1KG']
novelty = ['all', 'known', 'novel']
metrics = ['nEvalVariants', 'nSNPs', 'nInsertions', 'nDeletions',
           'nVariantsAtComp', 'compRate', 'nConcordant', 'concordantRate',
           'variantRate', 'variantRatePerBp', 'hetHomRatio', 'tiTvRatio']


def error(msg, prefix="Error"):
    sys.stderr.write(prefix + " " + msg + "\n")
    sys.stderr.flush()
    sys.exit(1)


def check_python_version():
    if sys.version[0:3] not in SUPPORTED_PYTHON_VERSIONS:
        error('Python version ' + sys.version[0:3] + ' is not supported!\n' + \
              'Supported versions are ' + ', '.join(SUPPORTED_PYTHON_VERSIONS))


def _call(cmdline, log_filename=None):
    print ''
    print '*' * 70
    print cmdline
    if log_filename:
        subprocess.call(cmdline.split(), stdout=open(log_filename, 'w'))
    else:
        subprocess.call(cmdline.split())


def run_gatk(input_filename, output_dir, gatk_dir, 
             ref_path, dbsnp_db, cosmic_db, oneKG_db, intervals_bed=None):
    report_filename = os.path.join(output_dir, 'gatk.report')
    log_filename = os.path.join(output_dir, 'gatk.log')

    gatk_jar = os.path.join(gatk_dir, 'GenomeAnalysisTK.jar')
    cmdline = 'java -Xmx2g -jar %s -nt 20 -R %s -T VariantEval ' \
              '-o %s --eval:sample %s -D %s -comp:Cosmic %s ' \
              '-comp:1KG %s ' % \
              (gatk_jar, ref_path, report_filename, input_filename,
               dbsnp_db, cosmic_db, oneKG_db)
    _call(cmdline, log_filename)
    return report_filename


def parse_gatk_report(report_filename):
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
        elif cur_metrics_ids: # process lines only if there are metrics in current section
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


def print_final_report(report_dict, report_filename, sample_name):
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
        out.write('  '.join('%-*s' % (col_width, value) for col_width, value in zip(col_widths, row)) + "\r\n")
    out.close()


def run_bcftools(input_filename, output_dir, bcftools_dir):
    gzipped_sample = os.path.join(output_dir, os.path.basename(input_filename) + '.gz')
    bgzip_cmdline = 'bgzip %s -c ' % \
                    (input_filename)
    _call(bgzip_cmdline, gzipped_sample)
    tabix_cmdline = 'tabix -p vcf %s ' % gzipped_sample
    _call(tabix_cmdline)
    
    bcftools_binary = os.path.join(bcftools_dir, 'bcftools')
    plotter_binary = os.path.join(bcftools_dir, 'plot-vcfstats')
    text_report_filename = os.path.join(output_dir, 'bcftools.report')
    viz_report_dir = os.path.join(output_dir, 'viz_report/')
    
    bcftools_cmdline = '%s stats %s ' % (bcftools_binary, gzipped_sample)
    _call(bcftools_cmdline, text_report_filename)

    plotter_cmdline = '%s -s %s -p %s --no-PDF ' % \
                      (plotter_binary, text_report_filename, viz_report_dir)
    _call(plotter_cmdline) # works only on Python v.2.5 and higher
    return viz_report_dir


def get_plots_from_bcftools(bcftools_report_dir, output_dir):
    original_plots_names = ['indels.0.png', 'substitutions.0.png']
    final_plots_names = ['indels.png', 'substitution.png']
    for id, original_plot in enumerate(original_plots_names):
        plot_src_filename = os.path.join(bcftools_report_dir, original_plot)
        plot_dst_filename = os.path.join(output_dir, final_plots_names[id])
        if os.path.exists(plot_src_filename):
            shutil.copyfile(plot_src_filename, plot_dst_filename)


if __name__ == '__main__':
    args = sys.argv[1:]

    if len(args) < 2:
        error('Usage: python ' + str(sys.argv[0]) + ' input.vcf output_dir [sample_name]', prefix='')
    check_python_version()

    input_filename = args[0]
    output_dir = os.path.abspath(args[1])
    if len(args) > 2:
        sample_name = args[2]
    else:
        sample_name = os.path.splitext(os.path.basename(input_filename))[0]

    ref_name = 'hg19'
    ref_path = '/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa'
    gatk_dir = '/group/ngs/src/CancerAnalysisPackage'
    bcftools_dir = '/group/ngs/src/bcftools'
    dbsnp_db = '/ngs/reference_data/genomes/Hsapiens/hg19/variation/dbsnp_138.vcf'
    cosmic_db = '/ngs/reference_data/genomes/Hsapiens/hg19/variation/cosmic-v67_20131024-hg19.vcf'
    oneKG_db = '/ngs/reference_data/genomes/Hsapiens/hg19/variation/1000G_omni2.5.vcf'
    #intervals_bed = '/ngs/reference_data/genomes/Hsapiens/staging/Ov_SS_bed_files/SureSelect/SS_regions_CP.bed'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    tmp_dir = tempfile.mkdtemp(dir=output_dir, prefix='variant_qc_tmp')

    gatk_report_filename = run_gatk(input_filename, tmp_dir, gatk_dir, ref_path,
                                    dbsnp_db, cosmic_db, oneKG_db)
    
    bcftools_viz_report = run_bcftools(input_filename, tmp_dir, bcftools_dir)

    gatk_report_filename = os.path.join(tmp_dir, 'gatk.report')
    gatk_report_dict = parse_gatk_report(gatk_report_filename)

    final_report_filename = os.path.join(output_dir, 'report.txt')
    print_final_report(gatk_report_dict, final_report_filename, sample_name)

    get_plots_from_bcftools(bcftools_viz_report, output_dir)

    print "\n\nText report is %s" % final_report_filename
    print "Plots are in %s" % output_dir

    if not DEBUG_MODE:
        shutil.rmtree(tmp_dir)

