from collections import OrderedDict, defaultdict
from itertools import izip
import json
from os.path import join, abspath, dirname, relpath, basename
from source import info
from source.clinical_reporting.clinical_parser import clinical_sample_info_from_bcbio_structure, Mutation, Patient
from source.clinical_reporting.clinical_reporting import Chromosome, BaseClinicalReporting
from source.file_utils import safe_mkdir
from source.logger import err
from source.reporting.reporting import MetricStorage, Metric, PerRegionSampleReport, write_static_html_report, \
    build_report_html, calc_cell_contents, make_cell_th, make_cell_td
from source.reporting.reporting import ReportSection
from source.utils import OrderedDefaultDict
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


def run_combine_clinical_reports(cnf, bcbio_structures):
    info('Running clinical reporting comparison')

    infos_by_key = dict()
    for i, bs in enumerate(bcbio_structures):
        for sample in bs.samples:
            if not cnf.sample_names or (cnf.sample_names and sample.name in cnf.sample_names):
                info('Preparing ' + sample.name + '...')
                info('-' * 70)
                clin_info = clinical_sample_info_from_bcbio_structure(cnf, bs, sample)
                infos_by_key[sample.name] = clin_info
                info('')

    info('*' * 70)
    run_sample_combine_clinreport(cnf, infos_by_key, cnf.output_dir)
    info('*' * 70)
    info('Successfully finished.')


def run_sample_combine_clinreport(cnf, infos_by_key, output_dirpath, is_target2wgs=False):
    ComparisonClinicalReporting(cnf, infos_by_key).write_report(
        join(output_dirpath, 'clinical_report.html'), is_target2wgs=is_target2wgs)


# class ComparisonMutation(Mutation):
#     def __init__(self, target_match=None, wgs_match=None, *args, **kwargs):
#         Mutation.__init__(self, *args, **kwargs)
#         self.target_match = target_match
#         self.wgs_match = wgs_match


# class MultiSampleMutation:
#     def __init__(self, muts_by_sample):
#         Mutation.__init__(self, muts_by_sample)


class ComparisonClinicalReporting(BaseClinicalReporting):
    def __init__(self, cnf, experiment_by_key, *args):
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

        sample_names = []
        for k, e in experiment_by_key.items():
            e.key = k
            sample_names.append(e.sample.name)
        # self.patient = self.merge_patients(self.infos)
        # bed_fpaths = set(experiment.target.bed_fpath for experiment in experiment_by_key.values() if experiment.target.bed_fpath)
        # bed_fnames = [basename(bed_fpath).split('.')[0] + '.bed' for bed_fpath in bed_fpaths]
        jbrowser_link = get_jbrowser_link(self.cnf.genome.name, sample_names)

        info('Preparing data...')
        # self.mut_by_key_by_exper = self.arrange_mutations({k: i.mutations for k, i in experiment_by_key.items()})
        mutations_by_experiment = {e: e.mutations for e in experiment_by_key.values() if e.mutations}
        if mutations_by_experiment:
            self.mutations_report, self.venn_plot_data = self.make_mutations_report(mutations_by_experiment, jbrowser_link, create_venn_diagrams=True)
            # self.mutations_plot_data = self.make_mutations_json(mutations_by_experiment)
            # self.substitutions_plot_data = self.make_substitutions_json(mutations_by_experiment)
        #self.actionable_genes_report = self.make_actionable_genes_report(experiment_by_key.values()[0].actionable_genes_dict)
        seq2c_events_by_experiment = {e: e.seq2c_events_by_gene for e in experiment_by_key.values() if e.seq2c_events_by_gene}
        if seq2c_events_by_experiment:
            # self.seq2c_plot_data = self.make_seq2c_plot_json(self.experiment_by_key)
            self.seq2c_report = self.make_seq2c_report(seq2c_events_by_experiment)
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

    def write_report(self, output_fpath, is_target2wgs=False):
        info('')

        data = {
            'key_or_target': self.experiment_by_key.values()[0].genes_collection_type,
            'genes_description': self.experiment_by_key.values()[0].genes_description,
            'sample': {
                'experiments': [self.sample_section(e, use_abs_report_fpath=True)
                                for k, e in self.experiment_by_key.items()],
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
        write_static_html_report(self.cnf, data, output_fpath,
           tmpl_fpath=join(dirname(abspath(__file__)), 'template_target2wgs.html' if is_target2wgs else 'template_combine.html'),
           extra_js_fpaths=[join(dirname(abspath(__file__)), 'static', 'clinical_report.js'),
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

            mutations_dict['total_variants'] = ', '.join(Metric.format_value(e.total_variants, is_html=True) + ' ' + k
                for k, e in experiment_by_key.items())

            mutations_dict['total_key_genes'] = '/'.join(
                set(Metric.format_value(len(e.key_gene_by_name_chrom.values()), is_html=True) + ' ' + k
                    for k, e in experiment_by_key.items()))

            mutations_dict['experiments'] = [dict(header=k, key=k.lower()) for k in self.experiment_by_key.keys()]

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
