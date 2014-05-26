import shutil
import textwrap
from os import mkdir, makedirs
from os.path import basename, join, isdir, dirname, expanduser

from source.bcbio_utils import file_exists
from source.utils import info, err, verify_file, verify_dir, verify_module, critical
from source.varqc.stats_gatk import gatk_qc

if verify_module('matplotlib'):
    import matplotlib
    matplotlib.use('Agg')  # non-GUI backend
    from source.varqc.distribution_plots import variants_distribution_plot
    from source.varqc.stats_bcftools import bcftools_qc
else:
    info('Warning: matplotlib is not installed, cannot draw plots.')


def run_qc(cnf, qc_dir, vcf_fpath):
    if 'quality_control' not in cnf:
        return None, None

    if not isdir(qc_dir):
        mkdir(qc_dir)

    qc_report_fpath = gatk_qc(cnf, qc_dir, vcf_fpath)
    if verify_module('matplotlib'):
        qc_plots_fpaths = bcftools_qc(cnf, qc_dir, vcf_fpath)
        qc_var_distr_plot_fpath = variants_distribution_plot(cnf, qc_dir, vcf_fpath)
        return qc_report_fpath, [qc_var_distr_plot_fpath] + qc_plots_fpaths
    else:
        return qc_report_fpath, None


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