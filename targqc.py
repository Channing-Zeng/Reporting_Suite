#!/usr/bin/env python

import sub_scripts.__common  # checking for python version and adding site dirs inside

import sys
import os
from optparse import OptionParser
from source.config import Config, defaults
from source.targetcov.summarize_targqc import summary_reports
from source.prepare_args_and_cnf import add_post_bcbio_args, detect_sys_cnf
from source.logger import info, err, warn, critical
from source.file_utils import verify_dir, safe_mkdir, adjust_path, verify_file, adjust_system_path
from source.main import check_genome_resources

import hashlib
import base64
from os.path import splitext, abspath, basename, join, isfile, dirname
from os import listdir
from source.bcbio_runner import Step, fix_bed_for_qualimap
from source.ngscat.bed_file import verify_bam, verify_bed
from source.bcbio_structure import BCBioStructure, Sample
from source.utils import get_system_path
from source.calling_process import call


def main():
    info(' '.join(sys.argv))
    info()

    description = 'This script generates target QC reports for each BAM provided as an input.'
    parser = OptionParser(description=description)
    add_post_bcbio_args(parser)
    parser.add_option('--work-dir', dest='work_dir', metavar='DIR')
    parser.add_option('--genome', dest='genome')
    parser.add_option('--only-summary', dest='only_summary', action='store_true')

    (opts, args) = parser.parse_args()

    output_dir = args[0] if len(args) > 0 else getcwd()
    output_dir = adjust_path(output_dir)
    info('Running on ' + output_dir)
    if not verify_dir(output_dir): sys.exit(1)
    bam_fpaths = [join(output_dir, fname) for fname in listdir(output_dir) if fname.endswith('.bam')]
    if not bam_fpaths: critical('No BAM files inside ' + output_dir)

    opts.sys_cnf = adjust_path(opts.sys_cnf) if opts.sys_cnf else detect_sys_cnf(opts)
    if not verify_file(opts.sys_cnf): sys.exit(1)
    info('Using ' + opts.sys_cnf)

    opts.run_cnf = adjust_path(opts.run_cnf) if opts.run_cnf else defaults['run_cnf']
    project_run_cnf_fpath = adjust_path(join(config_dirpath, basename(opts.run_cnf)))
    info('Using ' + opts.run_cnf + ', copying to ' + project_run_cnf_fpath)


    cnf = Config(opts.__dict__, opts.sys_cnf, opts.run_cnf)
    cnf.output_dir = output_dir

    cnf.qsub_runner = adjust_system_path(cnf.qsub_runner)
    if not cnf.qsub_runner: critical('Error: qsub-runner is not provided is sys-config.')
    if not verify_file(cnf.qsub_runner): sys.exit(1)

    check_genome_resources(cnf)

    if not cnf.work_dir: cnf.work_dir = join(cnf.output_dir, 'work')
    safe_mkdir(cnf.work_dir)
    cnf.log_dir = join(cnf.work_dir, 'log')
    safe_mkdir(cnf.log_dir)

    if not cnf.project_name: critical('Error: specify --project-name')

    info('*' * 70)
    info()

    run(cnf, bam_fpaths)


def _prep_steps(cnf, max_threads, threads_per_sample, bed_fpath):
    hasher = hashlib.sha1(cnf.output_dir)
    path_hash = base64.urlsafe_b64encode(hasher.digest()[0:4])[:-1]
    run_id = path_hash + '_' + cnf.project_name

    basic_params = \
        ' --sys-cnf ' + cnf.sys_cnf + \
        ' --run-cnf ' + cnf.run_cnf

    params_for_one_sample = basic_params + \
        ' -t ' + str(threads_per_sample) + \
       (' --reuse ' if cnf.reuse_intermediate else '') + \
        ' --genome ' + cnf.genome.name + \
        ' --project-name ' + cnf.project_name + ' ' \
        ' --log-dir ' + cnf.log_dir

    targetcov_params = params_for_one_sample + \
        ' -s {sample}' + \
        ' -o ' + join(cnf.output_dir, '{sample}_' + BCBioStructure.targetseq_name) + \
        ' --work-dir ' + join(cnf.work_dir, '{sample}_' + BCBioStructure.targetseq_name) + \
        ' --bam {bam}' + \
        ' --bed ' + cnf.bed + \
       ('--exons ' + cnf.exons if cnf.exons else '')

    targetcov_step = Step(cnf, run_id,
        name=BCBioStructure.targetseq_name, short_name='tc',
        interpreter='python',
        script=join('sub_scripts', 'targetcov.py'),
        paramln=targetcov_params
    )

    ngscat_params = params_for_one_sample + \
        ' -s {sample} ' + \
        ' -o ' + join(cnf.output_dir, '{sample}_' + BCBioStructure.ngscat_name) + \
        ' --work-dir ' + join(cnf.work_dir, '{sample}_' + BCBioStructure.ngscat_name) + \
        ' --bam {bam}' + \
        ' --bed ' + cnf.bed + \
        ' --saturation y '

    ngscat_step = Step(cnf, run_id,
        name=BCBioStructure.ngscat_name, short_name='nc',
        interpreter='python',
        script=join('sub_scripts', 'ngscat.py'),
        paramln=ngscat_params
    )

    qualimap_bed_fpath = join(cnf.work_dir, 'tmp_qualimap.bed')
    fix_bed_for_qualimap(bed_fpath, qualimap_bed_fpath)

    qualimap_params = \
        ' bamqc' + \
        ' -nt ' + str(threads_per_sample) + \
        ' --java-mem-size=24G' + \
        ' -nr 5000 ' + \
        ' -bam {bam}' + \
        ' -outdir ' + join(cnf.output_dir, '{sample}_' + BCBioStructure.qualimap_name) + \
        ' -gff ' + qualimap_bed_fpath + \
        ' -c' + \
        ' -gd HUMAN'

    qualimap_step = Step(cnf, run_id,
        name=BCBioStructure.qualimap_name, short_name='qm',
        script='qualimap',
        paramln=qualimap_params
    )

    summary_cmdline_params = ' '.join(sys.argv) + '--only-summary'

    targqc_summary_step = Step(
        cnf, run_id,
        name=BCBioStructure.targqc_name, short_name='targqc',
        interpreter='python',
        script=abspath(__file__),
        paramln=summary_cmdline_params
    )

    return targetcov_step, ngscat_step, qualimap_step, targqc_summary_step


def run(cnf, bam_fpaths):
    for bam_fpath in bam_fpaths:
        if not verify_bam(bam_fpath):
            sys.exit(1)

    bed_fpath = cnf.bed
    if not verify_bed(bed_fpath):
        sys.exit(1)

    samples = [Sample(basename(splitext(bam_fpath)[0]), bam=bam_fpath, bed=bed_fpath, genome=cnf.genome.name)
               for bam_fpath in bam_fpaths]

    max_threads = cnf.threads or 40
    threads_per_sample = max(max_threads / len(samples), 1)

    if not cnf.only_summary:
        targetcov_step, ngscat_step, qualimap_step, targqc_summary_step = \
            _prep_steps(cnf, max_threads, threads_per_sample, bed_fpath)

        summary_wait_for_steps = []

        for sample in samples:
            info('Processing "' + basename(sample.bam) + '"')

            info('TargetSeq for "' + basename(sample.bam) + '"')
            _submit_job(cnf, targetcov_step, sample.name, threads=threads_per_sample,
                bam=sample.bam, sample=sample.name)

            info('NgsCat for "' + basename(sample.bam) + '"')
            _submit_job(cnf, ngscat_step, sample.name, threads=threads_per_sample,
                bam=sample.bam, sample=sample.name)

            info('Qualimap for "' + basename(sample.bam) + '"')
            _submit_job(cnf, qualimap_step, sample.name, threads=threads_per_sample,
                bam=sample.bam, sample=sample.name)

            summary_wait_for_steps.append(targetcov_step.job_name(sample.name))
            summary_wait_for_steps.append(ngscat_step.job_name(sample.name))
            summary_wait_for_steps.append(qualimap_step.job_name(sample.name))

        _submit_job(cnf, targqc_summary_step, wait_for_steps=summary_wait_for_steps)

    else:
        info('Making targqc summary')
        _summary(cnf, samples, bed_fpath)


def _submit_job(cnf, step, sample_name='', wait_for_steps=None, threads=1, **kwargs):
    log_fpath = join(cnf.log_dir, (step.name + ('_' + sample_name if sample_name else '') + '.log'))

    if isfile(log_fpath):
        try:
            os.remove(log_fpath)
        except OSError:
            err('Warning: cannot remove log file ' + log_fpath + ', probably permission denied.')

    safe_mkdir(dirname(log_fpath))

    tool_cmdline = get_system_path(cnf, step.interpreter, step.script)
    if not tool_cmdline: sys.exit(1)
    cmdline = tool_cmdline + ' ' + step.param_line.format(**kwargs)

    hold_jid_line = '-hold_jid ' + ','.join(wait_for_steps or ['_'])
    job_name = step.job_name(sample_name)
    qsub = get_system_path(cnf, 'qsub')
    threads = str(threads)
    queue = cnf.queue
    runner_script = cnf.qsub_runner
    qsub_cmdline = (
        '{qsub} -pe smp {threads} -S /bin/bash -q {queue} '
        '-j n -o {log_fpath} -e {log_fpath} {hold_jid_line} '
        '-N {job_name} {runner_script} "{cmdline}"'.format(**locals()))

    info(step.name)
    info(qsub_cmdline)

    call(cnf, qsub_cmdline, silent=True)

    info()


from collections import OrderedDict
from os.path import relpath, join, exists
from os import listdir
import shutil
from source.reporting import SampleReport, FullReport, Metric, MetricStorage, ReportSection, write_tsv_rows, load_records
from source.logger import step_greetings, info, send_email, critical, warn
from source.targetcov import cov
from source.qualimap import report_parser as qualimap_report_parser
from source.ngscat import report_parser as ngscat_report_parser
from source.tools_from_cnf import get_system_path, get_qualimap_type
from source.calling_process import call
from source.file_utils import safe_mkdir, verify_file, verify_dir


def summary(cnf, samples, bed_fpath):
    step_greetings('Coverage statistics for all samples based on TargetSeq, ngsCAT, and Qualimap reports')

    targetcov_metric_storage = cov.header_metric_storage
    for depth in cnf.coverage_reports.depth_thresholds:
        name = 'Part of target covered at least by ' + str(depth) + 'x'
        targetcov_metric_storage.add_metric(
            Metric(name, short_name=str(depth) + 'x', description=name, unit='%'),
            'depth_metrics')

    targetcov_jsons_by_sample = dict(s.name, join(cnf.output_dir, s.name + BCBioStructure.targetseq_name, s.name + '.' + BCBioStructure.targetseq_name + '.json'))
    targetcov_htmls_by_sample = dict(s.name, join(cnf.output_dir, s.name + BCBioStructure.targetseq_name, s.name + '.' + BCBioStructure.targetseq_name + '.html'))
    ngscat_htmls_by_sample = dict(s.name, join(cnf.output_dir, s.name + BCBioStructure.ngscat_name, 'captureQC.html'))
    qualimap_htmls_by_sample = dict(s.name, join(cnf.output_dir, s.name + BCBioStructure.qualimap_name, 'qualimapReport.html'))

    all_htmls_by_sample = OrderedDict()
    for sample in bcbio_structure.samples:
        all_htmls_by_sample[sample.name] = OrderedDict()
        if sample.name in targetcov_htmls_by_sample:
            all_htmls_by_sample[sample.name]['targetcov'] = relpath(targetcov_htmls_by_sample[sample.name], cnf.output_dir)
        if sample.name in ngscat_htmls_by_sample:
            all_htmls_by_sample[sample.name]['ngscat'] = relpath(ngscat_htmls_by_sample[sample.name], cnf.output_dir)
        if sample.name in qualimap_htmls_by_sample:
            all_htmls_by_sample[sample.name]['qualimap'] = relpath(qualimap_htmls_by_sample[sample.name], cnf.output_dir)

    targqc_metric_storage = _get_targqc_metric_storage(OrderedDict(
        targetcov=targetcov_metric_storage,
        ngscat=ngscat_report_parser.metric_storage,
        qualimap=qualimap_report_parser.metric_storage))

    targqc_full_report = FullReport(cnf.name, [
        SampleReport(sample,
                     records=_get_targqc_records(OrderedDict(
                         targetcov=load_records(targetcov_jsons_by_sample[sample.name])
                         if sample.name in targetcov_jsons_by_sample else [],
                         ngscat=ngscat_report_parser.parse_ngscat_sample_report(ngscat_htmls_by_sample[sample.name])
                         if sample.name in ngscat_htmls_by_sample else [],
                         qualimap=qualimap_report_parser.parse_qualimap_sample_report(qualimap_htmls_by_sample[sample.name])
                         if sample.name in qualimap_htmls_by_sample else []
                     )),
                     html_fpath=all_htmls_by_sample[sample.name]) for sample in samples
            if sample.name in set(targetcov_jsons_by_sample.keys() +
                                  ngscat_htmls_by_sample.keys() +
                                  qualimap_htmls_by_sample.keys())],
        metric_storage=targqc_metric_storage)

    # Qualimap2 run for multi-sample plots
    if len(qualimap_htmls_by_sample):
        qualimap = get_system_path(cnf, interpreter=None, name='qualimap')
        if qualimap is not None and get_qualimap_type(qualimap) == "full":
            qualimap_output_dir = join(cnf.work_dir, 'qualimap_multi_bamqc')
            plots_dirpath = join(cnf.output_dir, 'plots')
            _correct_qualimap_genome_results(bcbio_structure)

            safe_mkdir(qualimap_output_dir)
            rows = []
            for sample_name, html_fpath in qualimap_htmls_by_sample.items():
                rows += [[sample_name, html_fpath]]
            data_file = write_tsv_rows(rows, qualimap_output_dir, 'qualimap_results_by_sample')
            cmdline = '{qualimap} multi-bamqc --data {data_file} -outdir {qualimap_output_dir}'.format(**locals())
            ret_code = call(cnf, cmdline, exit_on_error=False, return_err_code=True)
            targqc_full_report.plots = []
            qualimap_plots_dirpath = join(qualimap_output_dir, 'images_multisampleBamQcReport')
            if (ret_code is None or ret_code == 0) and verify_dir(qualimap_plots_dirpath):
                if exists(plots_dirpath):
                    shutil.rmtree(plots_dirpath)
                shutil.move(qualimap_plots_dirpath, plots_dirpath)
                for plot_fpath in listdir(plots_dirpath):
                    plot_fpath = join(plots_dirpath, plot_fpath)
                    if verify_file(plot_fpath) and plot_fpath.endswith('.png'):
                        targqc_full_report.plots.append(relpath(plot_fpath, cnf.output_dir))
            else:
                warn('Warning: Qualimap for multi-sample analysis failed to finish. TargQC will not contain plots.')
        else:
            warn('Warning: Qualimap for multi-sample analysis was not found. TargQC will not contain plots.')

    final_summary_report_fpaths = targqc_full_report.save_into_files(
        cnf.output_dir, BCBioStructure.targqc_name,
        'Coverage statistics for all samples based on TargetSeq, ngsCAT, and Qualimap reports')

    info()
    info('*' * 70)
    info('TargQC summary saved in: ')
    for fpath in final_summary_report_fpaths:
        if fpath: info('  ' + fpath)


qualimap_to_targetcov_dict = {'Number of reads': cov.header_metric_storage.get_metric('Reads'),
                              'Mapped reads': cov.header_metric_storage.get_metric('Mapped reads'),
                              'Unmapped reads': cov.header_metric_storage.get_metric('Unmapped reads'),
                              'Mapped reads (on target)': cov.header_metric_storage.get_metric('Reads mapped on target'),
                              'Coverage Mean': cov.header_metric_storage.get_metric('Average target coverage depth'),
                              'Coverage Standard Deviation': cov.header_metric_storage.get_metric('Std. dev. of target coverage depth')}

ngscat_to_targetcov_dict = {'Number reads': cov.header_metric_storage.get_metric('Mapped reads'),
                            '% target bases with coverage >= 1x': cov.header_metric_storage.get_metric('Percentage of target covered by at least 1 read'),
                            '% reads on target': cov.header_metric_storage.get_metric('Reads mapped on target'),
                            'mean coverage': cov.header_metric_storage.get_metric('Average target coverage depth')}


def _get_targqc_metric(metric, report_type='targetcov'):  # report type is in ['targetcov', 'qualimap', 'ngscat']
    if report_type == 'targetcov':
        return metric
    elif report_type == 'qualimap':
        if metric.name in qualimap_to_targetcov_dict.keys():
            return qualimap_to_targetcov_dict[metric.name]
        return metric
    elif report_type == 'ngscat':
        if metric.name in ngscat_to_targetcov_dict.keys():
            return ngscat_to_targetcov_dict[metric.name]
        return metric
    critical('Incorrect usage of get_targqc_metric(), report_type is %s but should be one of the following: %s' %
             (report_type, ", ".join(['targetcov', 'qualimap', 'ngscat'])))
    return None


def _get_targqc_metric_storage(metric_storages_by_report_type):
    class SectionId:
        def __init__(self, name, title):
            self.name = name
            self.title = title

        def __hash__(self):
            #return hash((self.name, self.title))
            return hash(self.name)  # use title from the first metric_storage

        def __eq__(self, other):
            #return (self.name, self.title) == (other.name, other.title)
            return self.name == other.name  # use title from the first metric_storage

    metrics_by_sections = OrderedDict()
    general_section_id = None
    general_section_metric_list = []
    for report_type, metric_storage in metric_storages_by_report_type.items():
        for section in metric_storage.sections:
            section_id = SectionId(section.name, section.title)
            if section_id not in metrics_by_sections.keys():
                metrics_by_sections[section_id] = []
            metrics_by_sections[section_id] += [metric for metric in
                                                metric_storage.get_metrics(sections=[section],
                                                                           skip_general_section=True)
                                                if metric == _get_targqc_metric(metric, report_type)]

        # specific behaviour for general section
        general_section_metric_list += [metric for metric in metric_storage.general_section.metrics
                                        if metric == _get_targqc_metric(metric, report_type)]
        if not general_section_id:
            general_section_id = SectionId(metric_storage.general_section.name, metric_storage.general_section.title)

    sections = []
    for section_id, metric_list in metrics_by_sections.items():
        sections.append(ReportSection(section_id.name, section_id.title, metric_list))

    return MetricStorage(general_section=ReportSection(general_section_id.name, general_section_id.title,
                                                       general_section_metric_list),
                         sections=sections)


def _get_targqc_records(records_by_report_type):
    targqc_records = []
    filled_metric_names = []
    for report_type, records in records_by_report_type.items():
        for record in records:
            new_metric = _get_targqc_metric(record.metric, report_type)
            if new_metric.name not in filled_metric_names:
                filled_metric_names.append(new_metric.name)
                record.metric = new_metric
                targqc_records.append(record)
    return targqc_records


# fixing java.lang.Double.parseDouble error on entries like "6,082.49"
def _correct_qualimap_genome_results(bcbio_structure):
    qualimap_results_txt_by_sample = bcbio_structure.get_qualimap_results_txt_fpath_by_sample()
    for sample_name, results_txt_fpath in qualimap_results_txt_by_sample.items():
        with open(results_txt_fpath, 'r') as f:
            content = f.readlines()
        with open(results_txt_fpath, 'w') as f:
            metrics_started = False
            for line in content:
                if ">> Reference" in line:
                    metrics_started = True
                if metrics_started:
                    line = line.replace(',', '')
                f.write(line)


if __name__ == '__main__':
    main()