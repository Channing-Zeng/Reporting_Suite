from collections import OrderedDict, defaultdict
from itertools import izip, chain
from json import load
import json
from os.path import join, dirname, abspath, relpath

import source
from source import verify_file, info
from source.clinical_reporting.solvebio_mutations import query_mutations
from source.logger import warn, err
from source.reporting.reporting import MetricStorage, Metric, PerRegionSampleReport, ReportSection, SampleReport, \
    calc_cell_contents, make_cell_td, write_static_html_report, make_cell_th, build_report_html
from source.targetcov.Region import SortableByChrom
from source.targetcov.flag_regions import get_depth_cutoff
from source.targetcov.summarize_targetcov import get_float_val, get_val
from source.utils import get_chr_lengths


class ClinicalReporting:
    def __init__(self, clinical_sample_info):
        self.info = clinical_sample_info

    def write_report(self):
        info('')
        info('Building report')
        total_variants = get_total_variants_number(self.sample, self.cnf.varqc_json_fpath)
        mutations_report = self.make_mutations_report(self.mutations)
        actionable_genes_dict = parse_broad_actionable()
        actionable_genes_report = self.make_actionable_genes_report(actionable_genes_dict)

        seq2c_events_by_gene_name, seq2c_plot_data = None, None
        if not self.cnf.seq2c_tsv_fpath:
            warn('No Seq2C results provided by option --seq2c, skipping plotting Seq2C')
        else:
            self.parse_seq2c_report(self.cnf.seq2c_tsv_fpath)
            seq2c_plot_data = self.make_seq2c_plot_json()

        key_genes_report = self.make_key_genes_cov_report(self.key_gene_by_name, self.ave_depth)
        cov_plot_data = self.make_key_genes_cov_json(self.key_gene_by_name)

        sample_dict = dict()
        sample_dict['sample'] = self.sample.name.replace('_', ' ')
        if self.patient.gender:
            sample_dict['patient'] = {'sex': self.patient.gender}
        sample_dict['project_name'] = self.cnf.project_name.replace('_', ' ')
        if self.cnf.project_report_path:
            sample_dict['project_report_rel_path'] = relpath(self.cnf.project_report_path, dirname(self.sample.clinical_html))
        sample_dict['panel'] = self.target.type
        sample_dict['bed_path'] = self.target.bed_fpath or ''
        if self.cnf.debug:
            sample_dict['panel'] = self.cnf.target_type + ', AZ300 IDT panel'
            sample_dict['bed_path'] = 'http://blue.usbod.astrazeneca.net/~klpf990/reference_data/genomes/Hsapiens/hg19/bed/Panel-IDT_PanCancer_AZSpike_V1.bed'

        sample_dict['sample_type'] = self.sample.normal_match if self.sample.normal_match else 'unpaired'  # plasma, unpaired'
        sample_dict['genome_build'] = self.cnf.genome.name
        sample_dict['target_type'] = self.target.type
        sample_dict['target_fraction'] = Metric.format_value(self.target.coverage_percent, is_html=True, unit='%')
        # approach_dict['min_depth'] = Metric.format_value(min_depth, is_html=True)
        sample_dict['ave_depth'] = Metric.format_value(self.ave_depth, is_html=True)

        mutations_dict = dict()
        if mutations_report.rows:
            # if cnf.debug:
            #     mutations_report.regions = mutations_report.regions[::20]
            mutations_dict['table'] = build_report_html(mutations_report, sortable=True)
            mutations_dict['total_variants'] = Metric.format_value(total_variants, is_html=True)
            mutations_dict['total_key_genes'] = Metric.format_value(len(self.key_gene_by_name), is_html=True)

        coverage_dict = dict(depth_cutoff=self.depth_cutoff, columns=[])
        GENE_COL_NUM = 3
        genes_in_col = len(key_genes_report.rows) / GENE_COL_NUM
        calc_cell_contents(key_genes_report, key_genes_report.get_rows_of_records())
        for i in range(GENE_COL_NUM):
            column_dict = dict()
            # column_dict['table'] = build_report_html(coverage_report)
            column_dict['metric_names'] = [make_cell_th(m) for m in key_genes_report.metric_storage.get_metrics()]
            column_dict['rows'] = [
                dict(records=[make_cell_td(r) for r in region.records])
                    for region in key_genes_report.rows[i * genes_in_col:(i+1) * genes_in_col]]
            coverage_dict['columns'].append(column_dict)
        coverage_dict['plot_data'] = cov_plot_data

        seq2c_dict = dict()
        if seq2c_plot_data:
            seq2c_dict['plot_data'] = seq2c_plot_data

        image_by_key = dict()
        # if seq2c_plot_fpath:
        #     image_by_key['seq2c_plot'] = seq2c_plot_fpath

        actionable_genes_dict = dict()
        if actionable_genes_report.rows:
            actionable_genes_dict['table'] = build_report_html(actionable_genes_report, sortable=False)

        write_static_html_report(self.cnf, {
            'sample': sample_dict,
            'variants': mutations_dict,
            'seq2c': seq2c_dict,
            'coverage': coverage_dict,
            'actionable_genes': actionable_genes_dict
        }, self.sample.clinical_html,
           tmpl_fpath=join(dirname(abspath(__file__)), 'template.html'),
           extra_js_fpaths=[join(dirname(abspath(__file__)), 'static', 'clinical_report.js'),
                            join(dirname(abspath(__file__)), 'static', 'draw_genes_coverage_plot.js'),
                            join(dirname(abspath(__file__)), 'static', 'draw_seq2c_plot.js')],
           extra_css_fpaths=[join(dirname(abspath(__file__)), 'static', 'clinical_report.css'),
                             join(dirname(abspath(__file__)), 'static', 'header_picture.css')],
           image_by_key=image_by_key)

        # clin_rep_symlink = adjust_path(join(sample.dirpath, '..', sample.name + '.clinical_report.html'))
        # if islink(clin_rep_symlink):
        #     os.unlink(clin_rep_symlink)
        # os.symlink(sample.clinical_html, clin_rep_symlink)

        info('Saved clinical report to ' + self.sample.clinical_html)
        info('-' * 70)
        info()
        return self.sample.clinical_html

    def make_actionable_genes_report(self, actionable_genes_dict):
        info('Preparing mutations stats for key gene tables')

        clinical_action_metric_storage = MetricStorage(
            sections=[ReportSection(metrics=[
                Metric('Gene', min_width=70, max_width=70),  # Gene & Transcript
                Metric('Variant', min_width=80, max_width=80, style='white-space: pre !important;', class_='long_line'),            # p.Glu82Lys
                Metric('Type', min_width=120, max_width=120, style='white-space: pre; !important', class_='long_line'),               # Frameshift
                Metric('Types of recurrent alterations', short_name='Types of recurrent\nalterations',
                       min_width=130, max_width=130, style='white-space: pre;'),  # Mutation
                Metric('Rationale', style='max-width: 300px !important; white-space: normal;'),          # Translocations predict sensitivity
                Metric('Therapeutic Agents'),  # Sorafenib
            ])])

        report = PerRegionSampleReport(sample=self.sample, metric_storage=clinical_action_metric_storage)
        actionable_gene_names = actionable_genes_dict.keys()

        sv_mutation_types = {'Rearrangement', 'Fusion'}
        cnv_mutation_types = {'Amplification', 'Deletion'}

        for gene in self.key_gene_by_name.values():
            if gene.name not in actionable_gene_names:
                continue
            possible_mutation_types = set(actionable_genes_dict[gene.name][1].split('; '))
            # possible_mutation_types = possible_mutation_types - sv_mutation_types
            # if not possible_mutation_types: continue

            variants = []
            types = []

            vardict_mut_types = possible_mutation_types - sv_mutation_types - cnv_mutation_types
            if vardict_mut_types:
                for mut in self.mutations:
                    if mut.gene.name == gene.name:
                        variants.append(mut.aa_change if mut.aa_change else '.')
                        types.append(mut.var_type)

            if cnv_mutation_types and gene.seq2c_event:
                if 'Amplification' in possible_mutation_types and gene.seq2c_event.amp_del == 'Amp' or \
                        'Deletion' in possible_mutation_types and gene.seq2c_event.amp_del == 'Del':
                    variants.append(gene.seq2c_event.amp_del + ', ' + gene.seq2c_event.fragment)
                    types.append(gene.seq2c_event.amp_del)

            if not variants:
                continue

            reg = report.add_row()
            reg.add_record('Gene', gene.name)
            reg.add_record('Variant', '\n'.join(variants))
            reg.add_record('Type', '\n'.join(types))
            reg.add_record('Types of recurrent alterations', actionable_genes_dict[gene.name][1].replace('; ', '\n'))
            reg.add_record('Rationale', actionable_genes_dict[gene.name][0])
            reg.add_record('Therapeutic Agents', actionable_genes_dict[gene.name][2])

        report.save_tsv(self.sample.clinical_target_tsv, human_readable=True)
        info('Saved report for actionable genes to ' + report.tsv_fpath)
        info('-' * 70)
        info()
        return report


    def make_seq2c_plot_json(self):
        chr_cum_lens = ClinicalReporting.Chromosome.get_cum_lengths(self.chromosomes_by_name)
        chr_cum_len_by_chrom = dict(zip([c.name for c in self.chromosomes_by_name.values()], chr_cum_lens))

        # {'nrm': {'xs': [], 'ys': []},
        #         'amp': {'xs': [], 'ys': [], 'gs': []},
        #         'del': {'xs': [], 'ys': [], 'gs': []},

        data = dict(
            events=[],
            ticksX=[[(chr_cum_lens[i] + chr_cum_lens[i + 1])/2, self.chromosomes_by_name.values()[i].short_name]
                    for i in range(len(self.chromosomes_by_name.keys()))],
            linesX=chr_cum_lens
        )

        for gene in self.key_gene_by_name.values():
            if gene.seq2c_event:
                data['events'].append(dict(
                    x=chr_cum_len_by_chrom[gene.chrom] + gene.start + (gene.end - gene.start) / 2,
                    geneName=gene.name,
                    logRatio=gene.seq2c_event.ab_log2r if gene.seq2c_event.ab_log2r is not None else gene.seq2c_event.log2r,
                    ampDel=gene.seq2c_event.amp_del,
                    fragment=gene.seq2c_event.fragment))

                # if not gene.seq2c_event.ab_log2r or gene.seq2c_event.fragment == 'BP':  # breakpoint, meaning part of exon is not amplified

                # gene.seq2c_event.log2r
                # if gene.seq2c_event.ab_log2r:
                #     if gene.seq2c_event.is_amp() or gene.seq2c_event.is_del():
                #         d = data[gene.seq2c_event.amp_del.lower()]
                #         d['xs'].append(x)
                #         d['ys'].append(gene.seq2c_event.ab_log2r)
                #         d['gs'].append(gene.name)
                #     else:
                #         warn('Event is not Amp or Del, it\'s ' + gene.seq2c_event.amp_del)

        data['maxY'] = max([e['logRatio'] for e in data['events']] + [2])  # max(chain(data['nrm']['ys'], data['amp']['ys'], data['del']['ys'], [2]))
        data['minY'] = min([e['logRatio'] for e in data['events']] + [-2])  # min(chain(data['nrm']['ys'], data['amp']['ys'], data['del']['ys'], [-2]))

        return json.dumps(data)


    def make_key_genes_cov_json(self, key_gene_by_name):
        chr_cum_lens = ClinicalReporting.Chromosome.get_cum_lengths(self.chromosomes_by_name)
        chr_cum_len_by_chrom = dict(zip([c.name for c in self.chromosomes_by_name.values()], chr_cum_lens))

        mut_info_by_gene = dict()
        for gene in key_gene_by_name.values():
            mut_info_by_gene[gene.name] = [('p.' + m.aa_change if m.aa_change else '.') for m in gene.mutations]
            if gene.seq2c_event and (gene.seq2c_event.is_amp() or gene.seq2c_event.is_del()):
                mut_info_by_gene[gene.name].append(gene.seq2c_event.amp_del + ', ' + gene.seq2c_event.fragment)

        ticks_x = [[(chr_cum_lens[i] + chr_cum_lens[i + 1])/2, self.chromosomes_by_name.values()[i].short_name]
                   for i in range(len(self.chromosomes_by_name.keys()))]

        gene_names = []
        gene_ave_depths = []
        cov_in_thresh = []
        coord_x = []
        cds_cov_by_gene = defaultdict(list)

        for gene in key_gene_by_name.values():
            gene_names.append(gene.name)
            gene_ave_depths.append(gene.ave_depth)
            cov_in_thresh.append(gene.cov_by_threshs.get(self.depth_cutoff))
            coord_x.append(chr_cum_len_by_chrom[gene.chrom] + gene.start + (gene.end - gene.start) / 2)
            cds_cov_by_gene[gene.name] = [dict(
                start=cds.start,
                end=cds.end,
                aveDepth=cds.ave_depth,
                percentInThreshold=cds.cov_by_threshs.get(self.depth_cutoff)
            ) for cds in gene.cdss]

        json_txt = '{"coords_x":%s, "ave_depths":%s, "cov_in_thresh":%s, "gene_names":%s, "mutations":%s, "ticks_x":%s, ' \
                    '"lines_x":%s, "cds_cov_by_gene":%s}'\
                   % (json.dumps(coord_x), json.dumps(gene_ave_depths), json.dumps(cov_in_thresh), json.dumps(gene_names),
                      json.dumps(mut_info_by_gene), json.dumps(ticks_x), json.dumps(chr_cum_lens), json.dumps(cds_cov_by_gene))
        return json_txt