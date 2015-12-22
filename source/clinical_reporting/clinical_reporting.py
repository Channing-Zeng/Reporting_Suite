from collections import OrderedDict, defaultdict
import json
from os.path import join, dirname, abspath, relpath

from source import info
from source.logger import warn
from source.reporting.reporting import MetricStorage, Metric, PerRegionSampleReport, ReportSection, calc_cell_contents, make_cell_td, write_static_html_report, make_cell_th, build_report_html
from source.utils import get_chr_lengths, OrderedDefaultDict


def make_clinical_report(cnf, clinical_sample_info, output_fpath):
    return ClinicalReporting(cnf, clinical_sample_info).write_report(output_fpath)


class Chromosome:
    def __init__(self, name, length=None):
        self.name = name
        self.short_name = name[3:] if name.startswith('chr') else name
        self.length = length

    @staticmethod
    def build_chr_by_name(cnf):
        chr_by_name = OrderedDict(
            (chr_name, Chromosome(chr_name, length=l)) for chr_name, l in get_chr_lengths(cnf).items()
            if '_' not in chr_name)  # not drawing extra chromosomes chr1_blablabla
        return chr_by_name

    @staticmethod
    def get_cum_lengths(chromosomes_by_name):
        return [sum(c.length for c in chromosomes_by_name.values()[:i])
            for i in range(len(chromosomes_by_name.keys()) + 1)]


class BaseClinicalReporting:
    def __init__(self, cnf, *args):
        self.cnf = cnf
        self.chromosomes_by_name = Chromosome.build_chr_by_name(self.cnf)

    def make_mutations_report(self, mutations_by_experiment):
        ms = [
            Metric('Gene'),  # Gene & Transcript
            Metric('Transcript'),  # Gene & Transcript
            # Metric('Codon chg', max_width=max_width, class_='long_line'),            # c.244G>A
            Metric('AA len', max_width=50, class_='stick_to_left', with_heatmap=False),          # 128
            Metric('AA chg', short_name='AA change', max_width=70, class_='long_line'),            # p.Glu82Lys
            # Metric('Allele'),             # Het.
            # Metric('Chr', max_width=33, with_heatmap=False),       # chr11
            Metric('Position', with_heatmap=False, align='left', sort_direction='ascending'),       # g.47364249
            Metric('Change', max_width=95, class_='long_line'),       # G>A
            # Metric('COSMIC', max_width=70, style='', class_='long_line'),                 # rs352343, COSM2123
            Metric('Effect', max_width=100, class_='long_line'),               # Frameshift
            Metric('VarDict status', short_name='Pathogenicity,\nReported by VarDict', class_='long_line'),     # Likely
            # Metric('VarDict reason', short_name='VarDict\nreason'),     # Likely
            Metric('Databases'),                 # rs352343, COSM2123
        ]
        # if is_local():
        #     ms.append(Metric('ClinVar', short_name='SolveBio ClinVar'))  # Pathogenic?, URL

        if len(mutations_by_experiment) == 1:
            ms.extend([
                Metric('Freq', short_name='Freq', max_width=55, unit='%', class_='shifted_column', with_heatmap=False),          # .19
                Metric('Depth', short_name='Depth', max_width=48, med=mutations_by_experiment.keys()[0].ave_depth, with_heatmap=False),              # 658
            ])
        else:
            for e in mutations_by_experiment.keys():
                ms.extend([
                    Metric(e.key + ' Freq', short_name=e.key + '\nfreq', max_width=55, unit='%',
                           class_='shifted_column', with_heatmap=False,
                           td_style='background-color: white'),          # .19
                    Metric(e.key + ' Depth', short_name='depth', max_width=48,
                           med=mutations_by_experiment.keys()[0].ave_depth, with_heatmap=False,
                           td_style='background-color: white'),              # 658
                ])

        clinical_mut_metric_storage = MetricStorage(sections=[ReportSection(metrics=ms)])
        report = PerRegionSampleReport(sample=mutations_by_experiment.keys()[0].sample,
            metric_storage=clinical_mut_metric_storage, expandable=True)

        # Writing records
        muts_by_key_by_experiment = OrderedDefaultDict(OrderedDict)
        for e, muts in mutations_by_experiment.items():
            for mut in muts:
                muts_by_key_by_experiment[mut.get_key()][e] = mut

        for mut_key, mut_by_experiment in muts_by_key_by_experiment.items():
            mut = next((m for m in mut_by_experiment.values() if m is not None), None)

            row = report.add_row()
            row.add_record('Gene', mut.gene.name)
            row.add_record('Transcript', mut.transcript)
            row.add_record('AA chg', **self._aa_chg_recargs(mut))
            row.add_record('Position', **self._pos_recargs(mut))
            row.add_record('Change', **self._chg_recargs(mut))
            row.add_record('AA len', mut.aa_len)
            row.add_record('Effect', mut.eff_type)
            row.add_record('VarDict status', **self._status_field(mut))
            # row.add_record('VarDict reason', mut.reason)
            row.add_record('Databases', **self._db_recargs(mut))
            # if is_local():
            #     row.add_record('ClinVar', **self._clinvar_recargs(mut))

            if len(mutations_by_experiment.values()) == 1:
                row.add_record('Freq', mut.freq if mut else None)
                row.add_record('Depth', mut.depth if mut else None)
            else:
                for e, m in mut_by_experiment.items():
                    row.add_record(e.key + ' Freq', m.freq if m else None)
                    row.add_record(e.key + ' Depth', m.depth if m else None)

            self._highlighting_and_hiding_mut_row(row, mut)
            if len(mut_by_experiment.keys()) == len(mutations_by_experiment.keys()):
                row.highlighted_green = True

        return report

    def make_seq2c_plot_json(self, experiment_by_key):
        data = dict()

        for k, e in experiment_by_key.items():
            chr_cum_lens = Chromosome.get_cum_lengths(self.chromosomes_by_name)
            chr_cum_len_by_chrom = dict(zip([c.name for c in self.chromosomes_by_name.values()], chr_cum_lens))

            d = dict(
                events=[],
                ticksX=[[(chr_cum_lens[i] + chr_cum_lens[i + 1])/2, self.chromosomes_by_name.values()[i].short_name]
                        for i in range(len(self.chromosomes_by_name.keys()))],
                linesX=chr_cum_lens
            )

            for gene in e.key_gene_by_name.values():
                if gene.seq2c_event:
                    d['events'].append(dict(
                        x=chr_cum_len_by_chrom[gene.chrom] + gene.start + (gene.end - gene.start) / 2,
                        geneName=gene.name,
                        logRatio=gene.seq2c_event.ab_log2r if gene.seq2c_event.ab_log2r is not None else gene.seq2c_event.log2r,
                        ampDel=gene.seq2c_event.amp_del,
                        fragment=gene.seq2c_event.fragment))

                    # if not gene.seq2c_event.ab_log2r or gene.seq2c_event.fragment == 'BP':  # breakpoint, meaning part of exon is not amplified

            d['maxY'] = max([e['logRatio'] for e in d['events']] + [2])  # max(chain(data['nrm']['ys'], data['amp']['ys'], data['del']['ys'], [2]))
            d['minY'] = min([e['logRatio'] for e in d['events']] + [-2])  # min(chain(data['nrm']['ys'], data['amp']['ys'], data['del']['ys'], [-2]))

            data[k.lower()] = d

        if len(experiment_by_key.keys()) == 1:
            return json.dumps(data.values()[0])
        else:
            return json.dumps(data)

    def make_key_genes_cov_json(self, experiment_by_key):
        chr_cum_lens = Chromosome.get_cum_lengths(self.chromosomes_by_name)
        chr_cum_len_by_chrom = dict(zip([c.name for c in self.chromosomes_by_name.values()], chr_cum_lens))

        gene_names = []
        coord_x = []

        ticks_x = [[(chr_cum_lens[i] + chr_cum_lens[i + 1])/2, self.chromosomes_by_name.values()[i].short_name]
                   for i in range(len(self.chromosomes_by_name.keys()))]

        hits = list()
        for key, e in experiment_by_key.items():
            gene_ave_depths = []
            covs_in_thresh = []
            cds_cov_by_gene = defaultdict(list)
            mut_info_by_gene = dict()

            for gene in e.key_gene_by_name.values():
                mut_info_by_gene[gene.name] = [('p.' + m.aa_change if m.aa_change else '.') for m in gene.mutations]
                if gene.seq2c_event and (gene.seq2c_event.is_amp() or gene.seq2c_event.is_del()):
                    mut_info_by_gene[gene.name].append(gene.seq2c_event.amp_del + ', ' + gene.seq2c_event.fragment)

            for gene in e.key_gene_by_name.values():
                gene_names.append(gene.name)
                gene_ave_depths.append(gene.ave_depth)
                covs_in_thresh.append(gene.cov_by_threshs.get(e.depth_cutoff))
                coord_x.append(chr_cum_len_by_chrom[gene.chrom] + gene.start + (gene.end - gene.start) / 2)
                cds_cov_by_gene[gene.name] = [dict(
                    start=cds.start,
                    end=cds.end,
                    aveDepth=cds.ave_depth,
                    percentInThreshold=cds.cov_by_threshs.get(e.depth_cutoff),
                ) for cds in gene.cdss]

            hits.append(dict(
                gene_ave_depths=gene_ave_depths,
                covs_in_thresh=covs_in_thresh,
                cds_cov_by_gene=cds_cov_by_gene,
                mut_info_by_gene=mut_info_by_gene
            ))

        data = dict(
            coords_x=coord_x,
            gene_names=gene_names,
            ticks_x=ticks_x,
            lines_x=chr_cum_lens,
        )

        if len(hits) == 1:
            for k, v in hits[0].items():
                data[k] = v
            data['color'] = None
        else:
            data['hits'] = hits

        return json.dumps(data)

    def sample_section(self, experiment):
        d = dict()
        d['patient'] = {'sex': 'unknown'}
        d['project_report_rel_path'] = 'not generagted'
        d['panel'] = 'unknown'
        d['bed_path'] = 'unknown'
        d['target_type'] = 'unknown'
        d['target_fraction'] = 'unknown'
        d['ave_depth'] = 'unknown'

        d['key'] = experiment.key
        d['sample'] = experiment.sample.name.replace('_', ' ')
        if experiment.patient and experiment.patient.gender:
            d['patient'] = {'sex': experiment.patient.gender}
        d['project_name'] = experiment.project_name.replace('_', ' ')
        if experiment.project_report_path:
            d['project_report_rel_path'] = relpath(experiment.project_report_path, dirname(experiment.sample.clinical_html))
        if experiment.target:
            d['panel'] = experiment.target.type
            d['bed_path'] = experiment.target.bed_fpath or ''
            d['target_type'] = experiment.target.type
            d['target_fraction'] = Metric.format_value(experiment.target.coverage_percent, is_html=True, unit='%')
            if self.cnf.debug:
                d['panel'] = experiment.target.type + ', AZ300 IDT panel'
                d['bed_path'] = 'http://blue.usbod.astrazeneca.net/~klpf990/reference_data/genomes/Hsapiens/hg19/bed/Panel-IDT_PanCancer_AZSpike_V1.bed'

        d['sample_type'] = experiment.sample.normal_match if experiment.sample.normal_match else 'unpaired'  # plasma, unpaired'
        d['genome_build'] = self.cnf.genome.name  # TODO: get genome build from the relevant project, not from the default config for this new run
        # approach_dict['min_depth'] = Metric.format_value(min_depth, is_html=True)
        if experiment.ave_depth is not None:
            d['ave_depth'] = Metric.format_value(experiment.ave_depth, is_html=True)
        return d

    @staticmethod
    def _db_recargs(mut):
        db = ''
        if mut.dbsnp_ids:
            db += ', '.join(('<a href="http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?searchType=adhoc_search'
                             '&type=rs&rs=' + rs_id + '">dbSNP</a>') for rs_id in mut.dbsnp_ids)
        if db and mut.cosmic_ids:
            db += ', '
        if mut.cosmic_ids:
            db += ', '.join('<a href="http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=' +
                            cid + '">COSM</a>' for cid in mut.cosmic_ids)
        return dict(value=db)

    @staticmethod
    def _aa_chg_recargs(mut):
        aa_chg = ''.join([gray(c) if c.isdigit() else c for c in (mut.aa_change or '')])
        return dict(value=aa_chg)

    @staticmethod
    def _pos_recargs(mut):
        c = (mut.chrom.replace('chr', '')) if mut.chrom else ''
        p = Metric.format_value(mut.pos, human_readable=True, is_html=True) if mut.pos else ''
        return dict(value=gray(c + ':') + p, num=mut.get_chrom_key() * 100000000000 + mut.pos)

    @staticmethod
    def _chg_recargs(mut):
        chg = mut.ref + '>' + mut.alt
        if mut.var_type:
            t = mut.var_type
            if t in ['Insertion', 'Deletion']:
                t = t[:3]
            chg = gray(t) + ' ' + chg
        return dict(value=chg)

    @staticmethod
    def _status_field(mut):
        status = mut.status
        if mut.reason:
            status += gray(' (' + mut.reason + ')')
        return dict(value=status)

    @staticmethod
    def _clinvar_recargs(mut):
        if mut.solvebio:
            return dict(value=mut.solvebio.clinsig, url=mut.solvebio.url)
        else:
            return dict(value='')

    @staticmethod
    def _highlighting_and_hiding_mut_row(row, mut):
        if not mut.status or mut.status.lower() == 'unknown':
            if mut.solvebio and 'Pathogenic' in mut.solvebio.clinsig:
                warn('Mutation ' + str(mut) + ' is unknown, but found in SolveBio')
            row.hidden = True
        else:
            if mut.status and mut.status.lower() == 'known':
                row.highlighted = True
        return row


class ClinicalReporting(BaseClinicalReporting):
    def __init__(self, cnf, clinical_experiment_info):
        BaseClinicalReporting.__init__(self, cnf)

        self.experiment = clinical_experiment_info
        self.sample = clinical_experiment_info.sample

        self.mutations_report = None
        self.actionable_genes_report = None
        self.seq2c_plot_data = None
        self.key_genes_report = None
        self.cov_plot_data = None

        info('Preparing data...')
        if self.experiment.mutations:
            self.mutations_report = self.make_mutations_report({self.experiment: self.experiment.mutations})
        if self.experiment.actionable_genes_dict:
            self.actionable_genes_report = self.make_actionable_genes_report(self.experiment.actionable_genes_dict)
        if self.experiment.seq2c_events_by_gene_name:
            self.seq2c_plot_data = self.make_seq2c_plot_json({self.experiment.key: self.experiment})
        if self.experiment.ave_depth:
            self.key_genes_report = self.make_key_genes_cov_report(self.experiment.key_gene_by_name, self.experiment.ave_depth)
            self.cov_plot_data = self.make_key_genes_cov_json({self.experiment.key: self.experiment})

    def write_report(self, output_fpath):
        info('')

        write_static_html_report(self.cnf, {
            'key_or_target': self.experiment.key_or_target_genes,
            'genes_description': self.experiment.genes_description,
            'sample': self.sample_section(self.experiment),
            'variants': self.__mutations_section(),
            'seq2c': {'plot_data': self.seq2c_plot_data},
            'coverage': self.__coverage_section(),
            'actionable_genes': self.__actionable_genes_section()
        }, output_fpath,
           tmpl_fpath=join(dirname(abspath(__file__)), 'template.html'),
           extra_js_fpaths=[join(dirname(abspath(__file__)), 'static', 'clinical_report.js'),
                            join(dirname(abspath(__file__)), 'static', 'draw_genes_coverage_plot.js'),
                            join(dirname(abspath(__file__)), 'static', 'draw_seq2c_plot.js')],
           extra_css_fpaths=[join(dirname(abspath(__file__)), 'static', 'clinical_report.css'),
                             join(dirname(abspath(__file__)), 'static', 'header_picture.css')])

        info('Saved clinical report to ' + output_fpath)
        info('-' * 70)
        info()
        return output_fpath

    def __mutations_section(self):
        mutations_dict = dict()
        if self.mutations_report and self.mutations_report.rows:
            # if cnf.debug:
            #     mutations_report.regions = mutations_report.regions[::20]
            mutations_dict['table'] = build_report_html(self.mutations_report, sortable=True)
            mutations_dict['total_variants'] = Metric.format_value(self.experiment.total_variants, is_html=True)
            mutations_dict['total_key_genes'] = Metric.format_value(len(self.experiment.key_gene_by_name), is_html=True)
        return mutations_dict

    def __coverage_section(self):
        if self.experiment.depth_cutoff is not None:
            coverage_dict = dict(depth_cutoff=self.experiment.depth_cutoff, columns=[])
            GENE_COL_NUM = 3
            genes_in_col = len(self.key_genes_report.rows) / GENE_COL_NUM
            calc_cell_contents(self.key_genes_report, self.key_genes_report.get_rows_of_records())
            for i in range(GENE_COL_NUM):
                column_dict = dict()
                # column_dict['table'] = build_report_html(coverage_report)
                column_dict['metric_names'] = [make_cell_th(m) for m in self.key_genes_report.metric_storage.get_metrics()]
                column_dict['rows'] = [
                    dict(records=[make_cell_td(r) for r in region.records])
                        for region in self.key_genes_report.rows[i * genes_in_col:(i+1) * genes_in_col]]
                coverage_dict['columns'].append(column_dict)
            coverage_dict['plot_data'] = self.cov_plot_data
            return coverage_dict
        else:
            return dict()

    def __actionable_genes_section(self):
        actionable_genes_dict = dict()
        if self.actionable_genes_report and self.actionable_genes_report.rows:
            actionable_genes_dict['table'] = build_report_html(self.actionable_genes_report, sortable=False)
        return actionable_genes_dict

    def make_key_genes_cov_report(self, key_gene_by_name, ave_depth):
        if self.experiment.depth_cutoff is None:
            return None

        info('Making ' + self.experiment.key_or_target_genes + ' genes coverage report...')
        clinical_cov_metrics = [
            Metric('Gene'),
            Metric('Chr', with_heatmap=False, max_width=20, align='right'),
            Metric('Ave depth', med=ave_depth),
            Metric('% cov at {}x'.format(self.experiment.depth_cutoff), unit='%', med=1, low_inner_fence=0.5, low_outer_fence=0.1),
            Metric('CNV', short_name='&nbsp;&nbsp;CNV')]  # short name is hack for IE9 who doesn't have "text-align: left" and tries to stick "CNV" to the previous col header

        clinical_cov_metric_storage = MetricStorage(sections=[ReportSection(metrics=clinical_cov_metrics)])

        key_genes_report = PerRegionSampleReport(sample=self.sample, metric_storage=clinical_cov_metric_storage)

        for gene in sorted(key_gene_by_name.values(), key=lambda g: g.name):
            reg = key_genes_report.add_row()
            reg.add_record('Gene', gene.name)
            reg.add_record('Chr', gene.chrom.replace('chr', ''))
            reg.add_record('Ave depth', gene.ave_depth)
            m = clinical_cov_metric_storage.find_metric('% cov at {}x'.format(self.experiment.depth_cutoff))
            reg.add_record(m.name, next((cov for cutoff, cov in gene.cov_by_threshs.items() if cutoff == self.experiment.depth_cutoff), None))
            if gene.seq2c_event and (gene.seq2c_event.is_amp() or gene.seq2c_event.is_del()):
                reg.add_record('CNV', gene.seq2c_event.amp_del + ', ' + gene.seq2c_event.fragment)

        return key_genes_report

    def make_actionable_genes_report(self, actionable_genes_dict):
        info('Preparing mutations stats for ' + self.experiment.key_or_target_genes + ' gene tables')

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

        for gene in self.experiment.key_gene_by_name.values():
            if gene.name not in actionable_gene_names:
                continue
            possible_mutation_types = set(actionable_genes_dict[gene.name][1].split('; '))
            # possible_mutation_types = possible_mutation_types - sv_mutation_types
            # if not possible_mutation_types: continue

            variants = []
            types = []

            vardict_mut_types = possible_mutation_types - sv_mutation_types - cnv_mutation_types
            if vardict_mut_types:
                for mut in self.experiment.mutations:
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

        return report

    # def make_seq2c_plot_json(self):
    #     chr_cum_lens = Chromosome.get_cum_lengths(self.chromosomes_by_name)
    #     chr_cum_len_by_chrom = dict(zip([c.name for c in self.chromosomes_by_name.values()], chr_cum_lens))
    #
    #     data = dict(
    #         events=[],
    #         ticksX=[[(chr_cum_lens[i] + chr_cum_lens[i + 1])/2, self.chromosomes_by_name.values()[i].short_name]
    #                 for i in range(len(self.chromosomes_by_name.keys()))],
    #         linesX=chr_cum_lens
    #     )
    #
    #     for gene in self.experiment.key_gene_by_name.values():
    #         if gene.seq2c_event:
    #             data['events'].append(dict(
    #                 x=chr_cum_len_by_chrom[gene.chrom] + gene.start + (gene.end - gene.start) / 2,
    #                 geneName=gene.name,
    #                 logRatio=gene.seq2c_event.ab_log2r if gene.seq2c_event.ab_log2r is not None else gene.seq2c_event.log2r,
    #                 ampDel=gene.seq2c_event.amp_del,
    #                 fragment=gene.seq2c_event.fragment))
    #
    #             # if not gene.seq2c_event.ab_log2r or gene.seq2c_event.fragment == 'BP':  # breakpoint, meaning part of exon is not amplified
    #
    #     data['maxY'] = max([e['logRatio'] for e in data['events']] + [2])  # max(chain(data['nrm']['ys'], data['amp']['ys'], data['del']['ys'], [2]))
    #     data['minY'] = min([e['logRatio'] for e in data['events']] + [-2])  # min(chain(data['nrm']['ys'], data['amp']['ys'], data['del']['ys'], [-2]))
    #
    #     return json.dumps(data)


def tooltip_long(string, max_len=30):
    if len(string) < max_len:
        return string
    else:
        return '<a class="tooltip-link" rel="tooltip" title="' + string + '">' + string[:max_len - 2] + '...</a>'


def gray(text):
    return '<span style="color: gray">' + text + '</span>'
