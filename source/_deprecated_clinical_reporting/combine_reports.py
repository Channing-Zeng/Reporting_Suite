from collections import OrderedDict, defaultdict
import json
from os.path import join, abspath, dirname, relpath, basename

import variant_filtering

import source
from source import info, verify_file
from source._deprecated_clinical_reporting.clinical_parser import clinical_sample_info_from_bcbio_structure, Patient, \
    get_sample_info, get_group_num, parse_mutations, get_mutations_fpath_from_bs
from source._deprecated_clinical_reporting.clinical_reporting import Chromosome, BaseClinicalReporting
from source._deprecated_clinical_reporting.combine_clinical_reporting_utils import format_experiment_names
from source.file_utils import file_transaction
from source.logger import err
from source.qsub_utils import wait_for_jobs
from source.reporting.reporting import MetricStorage, Metric, PerRegionSampleReport, write_static_html_report, \
    build_report_html, calc_cell_contents, make_cell_th, make_cell_td
from source.reporting.reporting import ReportSection
from source.targetcov.bam_and_bed_utils import sambamba_depth
from source.targetcov.cov import parse_sambamba_depth_output
from source.utils import OrderedDefaultDict, is_us
from tools.add_jbrowse_tracks import get_jbrowser_link


def run_clinical_target2wgs(cnf, wgs_bs, trg_bs, shared_sample_names, output_dirpath):
    info('Running clinical reporting comparison')

    for sname in shared_sample_names:
        info('Preparing ' + sname + '...')
        trg_sample = next(s for s in trg_bs.samples if s.name == sname)
        wgs_sample = next(s for s in wgs_bs.samples if s.name == sname)

        info('-' * 70)
        clin_trg_info = clinical_sample_info_from_bcbio_structure(cnf, trg_bs, trg_sample, is_target2wgs_comparison=True)
        info('')
        info('-' * 70)
        clin_wgs_info = clinical_sample_info_from_bcbio_structure(cnf, wgs_bs, wgs_sample, is_target2wgs_comparison=True)

        info('')
        info('*' * 70)
        infos_by_key = {'Target': clin_trg_info, 'WGS': clin_wgs_info}
        run_sample_combine_clinreport(cnf, infos_by_key, output_dirpath, is_target2wgs=True)
        info('*' * 70)
        info('Successfully finished.')


def run_combine_clinical_reports(cnf, bcbio_structures, parameters_info, samples_data):
    info('Running clinical reporting comparison')

    infos_by_key = OrderedDict()
    sample_names = [s.name for bs in bcbio_structures for s in bs.samples]
    for i, bs in enumerate(bcbio_structures):
        info()
        info('Preparing ' + bs.project_name + '...')
        info('-' * 70)
        for sample in bs.samples:
            if not cnf.sample_names or (cnf.sample_names and sample.name in cnf.sample_names):
                info('Preparing ' + sample.name + '...')
                info('-' * 70)
                # sample.targetcov_detailed_tsv = None
                clin_info = clinical_sample_info_from_bcbio_structure(cnf, bs, sample)
                group = samples_data[bs.bcbio_project_dirpath][sample.name].group
                uniq_key = get_uniq_sample_key(bs.project_name, sample, sample_names)
                infos_by_key[(group, uniq_key)] = clin_info
                info('')

        rejected_mutations = get_rejected_mutations(cnf, bs, clin_info.key_gene_by_name_chrom, clin_info.genes_collection_type)
        for sample in bs.samples:
            if not cnf.sample_names or (cnf.sample_names and sample.name in cnf.sample_names):
                group = samples_data[bs.bcbio_project_dirpath][sample.name].group
                uniq_key = get_uniq_sample_key(bs.project_name, sample, sample_names)
                infos_by_key[(group, uniq_key)].rejected_mutations = rejected_mutations[sample.name]

    sample_infos = OrderedDict({k: get_sample_info(e.sample.name, e.sample.dirpath, samples_data) for k, e in infos_by_key.iteritems()})
    sorted_sample_infos = sorted(sample_infos.items(), key=lambda x: ([x[1][j] for j in range(len(x[1]))], x[0][1]))
    sorted_experiments = OrderedDict()
    for k, v in sorted_sample_infos:
        sorted_experiments[k] = infos_by_key[k]
    save_all_mutations_depth(cnf, infos_by_key)

    info('*' * 70)
    run_sample_combine_clinreport(cnf, infos_by_key, cnf.output_dir, parameters_info, samples_data)
    info('*' * 70)
    info('Successfully finished.')


def get_uniq_sample_key(project_name, sample, sample_names):
    if sample_names and sample_names.count(sample.name) > 1:
        return sample.name + '_' + project_name
    return sample.name


def get_rejected_mutations_fpaths(pass_mutations_fpath):
    mut_basename = pass_mutations_fpath.split('.' + variant_filtering.mut_pass_suffix)[0]
    mut_reject_ending = variant_filtering.mut_reject_suffix + '.' + variant_filtering.mut_file_ext
    possible_reject_mutations_fpaths = [mut_basename + '.' + mut_reject_ending,
                                        mut_basename + '.' + variant_filtering.mut_paired_suffix + '.' + mut_reject_ending,
                                        mut_basename + '.' + variant_filtering.mut_single_suffix + '.' + mut_reject_ending]
    return possible_reject_mutations_fpaths


def get_rejected_mutations(cnf, bs, key_gene_by_name_chrom, genes_collection_type):
    rejected_mutations = defaultdict(dict)
    rejected_mutations_by_sample = defaultdict(list)

    pass_mutations_fpath, _ = get_mutations_fpath_from_bs(bs)

    for reject_mutations_fpath in get_rejected_mutations_fpaths(pass_mutations_fpath):
        if verify_file(reject_mutations_fpath, silent=True):
            info('Parsing rejected mutations from ' + str(reject_mutations_fpath))
            parse_mutations(cnf, None, key_gene_by_name_chrom, reject_mutations_fpath, genes_collection_type,
                            mutations_dict=rejected_mutations_by_sample)
            for sample, mutations in rejected_mutations_by_sample.iteritems():
                for mut in mutations:
                    rejected_mutations[sample][(mut.gene.name, mut.pos)] = mut
    return rejected_mutations


def save_all_mutations_depth(cnf, infos_by_key):
    mut_bed_fpath = join(cnf.work_dir, 'mutations.bed')

    if not cnf.reuse_intermediate or not verify_file(mut_bed_fpath):
        all_mutations_pos = defaultdict(set)
        for e in infos_by_key.values():
            for mut in e.mutations:
                all_mutations_pos[mut.chrom].add(mut.pos)
        with file_transaction(cnf.work_dir, mut_bed_fpath) as tx:
            with open(tx, 'w') as out_f:
                for chrom, positions in all_mutations_pos.iteritems():
                    for pos in positions:
                        out_f.write('\t'.join([chrom, str(pos - 1), str(pos)]) + '\n')

    sambamba_output_by_experiment = run_sambamba_use_grid(cnf, infos_by_key, mut_bed_fpath)

    for e, sambamba_output_fpath in sambamba_output_by_experiment.iteritems():
        regions = parse_sambamba_depth_output(e.sample.name, sambamba_output_fpath)
        depth_dict = defaultdict()
        for region in regions:
            depth_dict[region.end] = region.avg_depth
        e.mutations_depth = depth_dict


def run_sambamba_use_grid(cnf, infos_by_key, mut_bed_fpath):
    sambamba_output_by_experiment = dict()
    not_submitted_experiments = infos_by_key.values()
    while not_submitted_experiments:
        jobs_to_wait = []
        submitted_experiments = []
        reused_experiments = []

        for (group, uniq_key), e in infos_by_key.iteritems():
            if e not in not_submitted_experiments:
                continue
            sambamba_output_fpath = join(cnf.work_dir, uniq_key + '__mutations.bed')
            sambamba_output_by_experiment[e] = sambamba_output_fpath

            if cnf.reuse_intermediate and verify_file(sambamba_output_fpath, silent=True):
                info(sambamba_output_fpath + ' exists, reusing')
                reused_experiments.append(e)
                continue
            else:
                if not e.sample.bam:
                    err('Sample ' + e.sample.name + ' in ' + str(group) + ', ' + str(uniq_key) + ' has no BAM')
                    continue
                j = sambamba_depth(cnf, mut_bed_fpath, e.sample.bam, output_fpath=sambamba_output_fpath, only_depth=True, silent=True, use_grid=True)
                submitted_experiments.append(e)

                if not j.is_done:
                    jobs_to_wait.append(j)
                if len(jobs_to_wait) >= cnf.threads:
                    break
        if jobs_to_wait:
            info('Submitted ' + str(len(jobs_to_wait)) + ' jobs, waiting...')
            jobs_to_wait = wait_for_jobs(cnf, jobs_to_wait)
        else:
            info('No jobs to submit.')
        not_submitted_experiments = [e for e in not_submitted_experiments if
                                     e not in submitted_experiments and
                                     e not in reused_experiments]

    return sambamba_output_by_experiment


def run_sample_combine_clinreport(cnf, infos_by_key, output_dirpath, parameters_info=None, samples_data=None, is_target2wgs=False):
    report = ComparisonClinicalReporting(cnf, infos_by_key, parameters_info, samples_data)
    report.mutations_parameters = make_parameters_json(parameters_info, samples_data, report.experiment_by_key)
    report_sample_names = get_sample_names_for_report(report.experiment_by_key, samples_data)
    report_sample_names = OrderedDict({experiment: '<a href="' + basename(experiment.sample.clinical_html) + '" target="_blank">' + name + '</a>'
                        for experiment, name in report_sample_names.iteritems()})
    report.write_report(
        join(output_dirpath, 'clinical_report.html'), sample_names=report_sample_names, is_target2wgs=is_target2wgs)
    report.key_genes_report = None
    sample_nums = set([get_group_num(key) for key in infos_by_key.keys()])
    for group_num in sample_nums:
        sample_report = ComparisonClinicalReporting(cnf, dict())
        sample_experiments = OrderedDict()
        for k, e in report.experiment_by_key.items():
            if get_group_num(e.key) == group_num and 'PBMC' not in e.sample.name:
                sample_experiments[k] = e

        sample_report.experiment_by_key = sample_experiments
        sample_report.mutations_report, sample_report.venn_plot_data = report.mutations_reports[group_num]
        sample_report.mutations_parameters = make_parameters_json(parameters_info, samples_data, sample_experiments, group_num)
        sample_report.seq2c_report = report.seq2c_reports[group_num]
        sample_names = get_sample_names_for_report(sample_experiments, samples_data, cur_group_num=group_num)
        sample_report.write_report(join(output_dirpath, 'report_' + str(group_num) + '.html'), sample_names=sample_names)


def get_sample_names_for_report(report_experiments, samples_data, cur_group_num=None):
    mutations_by_experiment = OrderedDict()
    for e in report_experiments.values():
        mutations_by_experiment[e] = e.mutations
    _, sample_names = format_experiment_names(mutations_by_experiment, samples_data, cur_group_num=cur_group_num)
    return sample_names


def make_parameters_json(parameters_info, samples_data, sample_experiments, group_num=None):
    parameters_list = []
    parameters_values = defaultdict(set)
    for k, e in sample_experiments.iteritems():
        if group_num and get_group_num(e.key) != group_num:
            continue
        project_dirpath = dirname(dirname(e.sample.dirpath))
        sample_info = samples_data[project_dirpath][e.sample.name].data
        for parameter_name, value in sample_info.iteritems():
            parameters_values[parameter_name].add(value)
    for parameter_name, parameter in parameters_info.iteritems():
        if len(parameters_values[parameter_name]) < 2:
            continue
        parameter_info = dict()
        parameter_info['parameter'] = parameter_name
        sorted_values = sorted(list(parameter.values))
        parameter_info['values'] = ['all'] + sorted_values + [', '.join(sorted_values)]
        parameter_info['valuesNames'] = ['all'] + sorted_values + ['common']
        parameters_list.append(parameter_info)
    return json.dumps(parameters_list)

# class ComparisonMutation(Mutation):
#     def __init__(self, target_match=None, wgs_match=None, *args, **kwargs):
#         Mutation.__init__(self, *args, **kwargs)
#         self.target_match = target_match
#         self.wgs_match = wgs_match


# class MultiSampleMutation:
#     def __init__(self, muts_by_sample):
#         Mutation.__init__(self, muts_by_sample)


class ComparisonClinicalReporting(BaseClinicalReporting):
    def __init__(self, cnf, experiment_by_key, parameters_info=None, samples_data=None, *args):
        BaseClinicalReporting.__init__(self, cnf, *args)

        self.experiment_by_key = experiment_by_key

        self.mutations_report = None
        self.mutations_plot_data = None
        self.venn_plot_data = None
        self.substitutions_plot_data = None
        self.sv_report = None
        self.actionable_genes_report = None
        self.seq2c_plot_data = None
        self.seq2c_report = None
        self.key_genes_report = None
        self.cov_plot_data = None
        self.mutations_by_experiment = OrderedDict()
        self.mutations_reports = dict()
        self.mutations_parameters = None
        self.seq2c_reports = defaultdict()

        self.sample_names = []
        for k, e in experiment_by_key.items():
            e.key = k
            import re
            e.sample.clinical_html = abspath(join(cnf.output_dir, 'report_' + str(get_group_num(k)) + '.html'))
            e.cnf.work_dir = cnf.work_dir
            if is_us:
                e.sample.clinical_html = re.sub('^/ngs/', '/gpfs/ngs/', e.sample.clinical_html)
                e.project_report_path = re.sub('^/ngs/', '/gpfs/ngs/', e.project_report_path)
            self.sample_names.append(e.sample.name)
        sample_infos = OrderedDict({k: get_sample_info(e.sample.name, e.sample.dirpath, samples_data) for k, e in experiment_by_key.iteritems()})
        sorted_sample_infos = sorted(sample_infos.items(), key=lambda x: (x[0][1], [x[1][j] for j in range(len(x[1]))]))
        sorted_experiments = OrderedDict()
        for k, v in sorted_sample_infos:
            sorted_experiments[k] = experiment_by_key[k]

        self.experiment_by_key = sorted_experiments
        # self.patient = self.merge_patients(self.infos)
        # bed_fpaths = set(experiment.target.bed_fpath for experiment in experiment_by_key.values() if experiment.target.bed_fpath)
        # bed_fnames = [basename(bed_fpath).split('.')[0] + '.bed' for bed_fpath in bed_fpaths]
        jbrowser_link = get_jbrowser_link(self.cnf.genome.name, self.sample_names)

        info('Preparing data...')
        # self.mut_by_key_by_exper = self.arrange_mutations({k: i.mutations for k, i in experiment_by_key.items()})
        for e in sorted_experiments.values():
            if e.mutations:
                self.mutations_by_experiment[e] = e.mutations
        group_nums = set(get_group_num(key) for key in self.experiment_by_key.keys())
        if self.mutations_by_experiment:
            self.mutations_report, self.venn_plot_data = self.make_mutations_report(self.mutations_by_experiment, jbrowser_link,
                                                                                    samples_data=samples_data, parameters_info=parameters_info,
                                                                                    create_venn_diagrams=True)
            info('Preparing data for each sample...')
            for num in group_nums:
                sample_mut_report, venn_plot_data = self.make_mutations_report(self.mutations_by_experiment, jbrowser_link,
                                                                               samples_data=samples_data, parameters_info=parameters_info,
                                                                               create_venn_diagrams=True, cur_group_num=num)
                self.mutations_reports[num] = (sample_mut_report, venn_plot_data)
            # self.mutations_plot_data = self.make_mutations_json(mutations_by_experiment)
            # self.substitutions_plot_data = self.make_substitutions_json(mutations_by_experiment)
        #self.actionable_genes_report = self.make_actionable_genes_report(experiment_by_key.values()[0].actionable_genes_dict)
        seq2c_events_by_experiment = {e: e.seq2c_events_by_gene for e in experiment_by_key.values() if e.seq2c_events_by_gene}
        if seq2c_events_by_experiment:
            for num in group_nums:
                seq2c_report = self.make_seq2c_report(seq2c_events_by_experiment, samples_data=samples_data, cur_group_num=num)
                self.seq2c_reports[num] = seq2c_report
            # self.seq2c_plot_data = self.make_seq2c_plot_json(self.experiment_by_key)
            # self.seq2c_report = self.make_seq2c_report(seq2c_events_by_experiment)
        # self.key_genes_report = self.make_key_genes_cov_report(self.experiment_by_key)
        # self.cov_plot_data = self.make_key_genes_cov_json(self.experiment_by_key)


    @staticmethod
    def merge_patients(patients):
        gender = None
        genders = set(p.gender for p in patients if p.gender)
        if genders:
            if len(genders) > 1:
                err('Different genders detected for the same sample: ' + str(genders))
            gender = next(genders)
        return Patient(gender)

    # @staticmethod
    # def arrange_mutations(muts_by_experiment):
    #     muts_by_key_by_exper = OrderedDefaultDict(OrderedDict)
    #
    #     for experiment_key, muts in muts_by_experiment.items():
    #         for mut in muts:
    #             muts_by_key_by_exper[mut.get_key()][experiment_key] = mut
    #
    #     return muts_by_key_by_exper

    def write_report(self, output_fpath, is_target2wgs=False, sample_names=None):
        info('')

        data = {
            'key_or_target': self.experiment_by_key.values()[0].genes_collection_type,
            'genes_description': self.experiment_by_key.values()[0].genes_description,
            'sample': {
                'experiments': [self.sample_section(e, sample_name=e.sample.name + (', ' + sample_names[e] if
                                sample_names and e in sample_names else '')) for k, e in self.experiment_by_key.items()],
            },
            # 'patient': self.__patient_section(self.patient),
            # 'sample_name': self.sample_name,
            'variants': self.__mutations_section(self.mutations_report, self.experiment_by_key),
            'coverage': self.__coverage_section(self, self.key_genes_report, self.cov_plot_data),
            # 'actionable_genes': self.__actionable_genes_section()
        }

        min_af = self.cnf.min_af or 0
        data['min_af'] = str(float(min_af) * 100)
        if self.seq2c_report:
            data['seq2c'] = {'amp_del': self.seq2c_section()}
        if self.seq2c_plot_data:
            data['seq2c']['plot'] = {'plot_data': self.seq2c_plot_data}
        if data['variants']:
            data['variants']['venn_diagram'] = {'diagram_data': self.venn_plot_data}
            data['variants']['mut_parameters'] = self.mutations_parameters
        write_static_html_report(self.cnf, data, output_fpath,
           tmpl_fpath=join(dirname(abspath(__file__)), 'template_target2wgs.html' if is_target2wgs else 'template_combine.html'),
           extra_js_fpaths=[join(dirname(abspath(__file__)), 'static', 'clinical_report.js'),
                            join(dirname(abspath(__file__)), 'static', 'combined_clinical_report.js'),
                            join(dirname(abspath(__file__)), 'static', 'draw_genes_coverage_plot.js'),
                            join(dirname(abspath(__file__)), 'static', 'draw_mutations_plot.js'),
                            join(dirname(abspath(__file__)), 'static', 'd3.min.js'),
                            join(dirname(abspath(__file__)), 'static', 'venn.js'),
                            join(dirname(abspath(__file__)), 'static', 'draw_venn_diagram.js'),
                            join(dirname(abspath(__file__)), 'static', 'draw_substitutions_plot.js'),
                            join(dirname(abspath(__file__)), 'static', 'draw_seq2c_plot.js')],
           extra_css_fpaths=[join(dirname(abspath(__file__)), 'static', 'clinical_report.css'),
                             join(dirname(abspath(__file__)), 'static', 'header_picture.css')])

        info('Saved clinical report to ' + output_fpath)
        info('-' * 70)
        info()
        return output_fpath

    # def __patient_section(self, patient):
    #     patient_dict = dict()
    #     if self.patient.gender:
    #         patient_dict['gender'] = self.patient.gender
    #     return patient_dict

    def __mutations_section(self, mutations_report, experiment_by_key):
        mutations_dict = dict()
        if self.mutations_report and self.mutations_report.rows:
            # if cnf.debug:
            #     mutations_report.regions = mutations_report.regions[::20]
            mutations_dict['table'] = build_report_html(mutations_report, sortable=True)

            mutations_dict['total_variants'] = ', '.join(Metric.format_value(e.total_variants, is_html=True) + ' ' + k[1]
                for k, e in experiment_by_key.items())

            mutations_dict['total_key_genes'] = '/'.join(
                set(Metric.format_value(len(e.key_gene_by_name_chrom.values()), is_html=True) + ' ' + k[1]
                    for k, e in experiment_by_key.items()))

            mutations_dict['experiments'] = [dict(header=k[1], key=k[1].lower()) for k in self.experiment_by_key.keys()]

            mutations_dict['plot_data'] = self.mutations_plot_data
            mutations_dict['substitutions_plot_data'] = self.substitutions_plot_data

        return mutations_dict

    @staticmethod
    def __coverage_section(self, key_genes_report, cov_plot_data):
        if not key_genes_report:
            return None
        coverage_dict = dict(columns=[])
        GENE_COL_NUM = 2
        genes_in_col = [len(key_genes_report.rows) / GENE_COL_NUM] * GENE_COL_NUM
        for i in range(len(key_genes_report.rows) % GENE_COL_NUM):
            genes_in_col[i] += 1
        calc_cell_contents(key_genes_report, key_genes_report.get_rows_of_records())
        printed_genes = 0
        for i in range(GENE_COL_NUM):
            column_dict = dict()
            # column_dict['table'] = build_report_html(coverage_report)
            column_dict['metric_names'] = [make_cell_th(m) for m in key_genes_report.metric_storage.get_metrics()]
            column_dict['rows'] = [
                dict(records=[make_cell_td(r) for r in region.records])
                    for region in key_genes_report.rows[printed_genes:printed_genes + genes_in_col[i]]]
            coverage_dict['columns'].append(column_dict)
            printed_genes += genes_in_col[i]
        coverage_dict['experiments'] = [dict(header=k, key=k.lower()) for k in self.experiment_by_key.keys()]
        coverage_dict['plot_data'] = cov_plot_data
        return coverage_dict

    @staticmethod
    def make_key_genes_cov_report(experiment_by_key):
        info('Making key genes coverage report...')

        ms = [
            Metric('Gene'),
            Metric('Chr', with_heatmap=False, max_width=20, align='right')]

        for i, (k, e) in enumerate(experiment_by_key.items()):
            ms.extend([
                Metric(k + ' Ave depth', short_name=k + '\nave depth', med=e.ave_depth, class_='shifted_column' if i == 0 else ''),
                Metric(k + ' % cov at {}x'.format(e.depth_cutoff),  short_name='% at {}x'.format(e.depth_cutoff), unit='%', med=1, low_inner_fence=0.5, low_outer_fence=0.1),
                Metric(k + ' CNV', short_name='&nbsp;&nbsp;CNV')]  # short name is hack for IE9 who doesn't have "text-align: left" and tries to stick "CNV" to the previous col header
            )
        clinical_cov_metric_storage = MetricStorage(sections=[ReportSection(metrics=ms)])
        key_genes_report = PerRegionSampleReport(sample=experiment_by_key.values()[0].sample, metric_storage=clinical_cov_metric_storage)

        # Writing records
        hits_by_gene_by_experiment = OrderedDefaultDict(OrderedDict)
        for k, e in experiment_by_key.items():
            for gene in e.key_gene_by_name.values():
                hits_by_gene_by_experiment[gene.name][e] = gene

        for gname, hit_by_experiment in sorted(hits_by_gene_by_experiment.items(), key=lambda (gname, h): gname):
            gene = next((m for m in hit_by_experiment.values() if m is not None), None)

            row = key_genes_report.add_row()
            row.add_record('Gene', gene.name)
            row.add_record('Chr', gene.chrom.replace('chr', ''))

            for e, hit in hit_by_experiment.items():
                row.add_record(e.key + ' Ave depth', hit.ave_depth)
                m = clinical_cov_metric_storage.find_metric(e.key + ' % cov at {}x'.format(e.depth_cutoff))
                row.add_record(m.name, next((cov for cutoff, cov in hit.cov_by_threshs.items() if cutoff == e.depth_cutoff), None))
                if hit.seq2c_event and (hit.seq2c_event.is_amp() or hit.seq2c_event.is_del()):
                    row.add_record(e.key + ' CNV', hit.seq2c_event.amp_del + ', ' + hit.seq2c_event.fragment)

        return key_genes_report


    # def make_key_genes_cov_json(self, experiment_by_key):
    #     chr_cum_lens = Chromosome.get_cum_lengths(self.chromosomes_by_name)
    #     chr_cum_len_by_chrom = dict(zip([c.name for c in self.chromosomes_by_name.values()], chr_cum_lens))
    #
    #     gene_names = []
    #     coord_x = []
    #
    #     ticks_x = [[(chr_cum_lens[i] + chr_cum_lens[i + 1])/2, self.chromosomes_by_name.values()[i].short_name]
    #                for i in range(len(self.chromosomes_by_name.keys()))]
    #
    #     hits = list()
    #     for key, e in experiment_by_key.items():
    #         gene_ave_depths = []
    #         covs_in_thresh = []
    #         cds_cov_by_gene = defaultdict(list)
    #         mut_info_by_gene = dict()
    #
    #         for gene in e.key_gene_by_name.values():
    #             mut_info_by_gene[gene.name] = [('p.' + m.aa_change if m.aa_change else '.') for m in gene.mutations]
    #             if gene.seq2c_event and (gene.seq2c_event.is_amp() or gene.seq2c_event.is_del()):
    #                 mut_info_by_gene[gene.name].append(gene.seq2c_event.amp_del + ', ' + gene.seq2c_event.fragment)
    #
    #         for gene in e.key_gene_by_name.values():
    #             gene_names.append(gene.name)
    #             gene_ave_depths.append(gene.ave_depth)
    #             covs_in_thresh.append(gene.cov_by_threshs.get(e.depth_cutoff))
    #             coord_x.append(chr_cum_len_by_chrom[gene.chrom] + gene.start + (gene.end - gene.start) / 2)
    #             cds_cov_by_gene[gene.name] = [dict(
    #                 start=cds.start,
    #                 end=cds.end,
    #                 aveDepth=cds.ave_depth,
    #                 percentInThreshold=cds.cov_by_threshs.get(e.depth_cutoff),
    #             ) for cds in gene.cdss]
    #
    #         hits.append(dict(
    #             gene_ave_depths=gene_ave_depths,
    #             covs_in_thresh=covs_in_thresh,
    #             cds_cov_by_gene=cds_cov_by_gene,
    #             mut_info_by_gene=mut_info_by_gene
    #         ))
    #
    #     data = dict(
    #         coords_x=coord_x,
    #         gene_names=gene_names,
    #         ticks_x=ticks_x,
    #         lines_x=chr_cum_lens,
    #     )
    #
    #     if len(hits) == 1:
    #         for k, v in hits[0].items():
    #             data[k] = v
    #         data['color'] = None
    #     else:
    #         data['hits'] = hits
    #
    #     return json.dumps(data)
