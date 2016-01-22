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
            Metric('Change', max_width=95, class_='long_line', description='cDNA change'),       # G>A
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

    def make_mutations_json(self, mutations_by_experiment):
        data = dict()

        for e, muts in mutations_by_experiment.items():
            d = dict(
                mutations=[]
            )

            for i, mut in enumerate(sorted(muts, key=lambda m: m.freq, reverse=True)):
                d['mutations'].append(dict(
                    x=i+1,
                    geneName=mut.gene.name,
                    chrom=mut.chrom, position=mut.pos, freq=mut.freq * 100,
                    mutType=mut.eff_type, aaChg=mut.aa_change, cdnaChange=mut.codon_change))
            d['minY'] = 0

            data[e.key.lower()] = d

        if len(mutations_by_experiment.keys()) == 1:
            return json.dumps(data.values()[0])
        else:
            return json.dumps(data)

    def make_substitutions_json(self, mutations_by_experiment):
        data = dict()
        nucleotides = ['A', 'C', 'G', 'T']

        def _add_nuc(nuc, substitutions):
            for nuc2 in nucleotides:
                if nuc != nuc2:
                    substitutions[nuc + '>' + nuc2] = 0

        for e, muts in mutations_by_experiment.items():
            d = dict(
                substitutions = OrderedDict()
            )
            for nuc in nucleotides:
                _add_nuc(nuc, d['substitutions'])
            for mut in muts:
                if mut.var_type == 'SNV' and mut.alt in nucleotides:
                    d['substitutions'][mut.ref + '>' + mut.alt] += 1
            d['maxY'] = max([d['substitutions'][s] for s in d['substitutions']])
            d['maxRate'] = d['maxY'] * 100 / sum([d['substitutions'][s] for s in d['substitutions']])
            d['minY'] = 0

            data[e.key.lower()] = d

        if len(mutations_by_experiment.keys()) == 1:
            return json.dumps(data.values()[0])
        else:
            return json.dumps(data)

    @staticmethod
    def make_sv_report(svs_by_experiment):
        ms = [
            Metric('Chr', with_heatmap=False, max_width=50, align='left', sort_direction='ascending'),
            Metric('Type'),
            Metric('Location', with_heatmap=False, align='left', sort_direction='ascending'),
            Metric('Genes', max_width=200, class_='long_line'),
            Metric('Status', with_heatmap=False, align='left', sort_direction='ascending'),
        ]

        if len(svs_by_experiment) == 1:
            ms.extend([
                Metric('Split read support', short_name='Split /', with_heatmap=False),
                Metric('Paired read support', short_name='paired read support', with_heatmap=False),
            ])
        else:
            for e in svs_by_experiment.keys():
                ms.extend([
                    Metric(e.key + ' Split read support', short_name='Split /', with_heatmap=False),
                    Metric(e.key + ' Paired read support', short_name='paired read support', with_heatmap=False),
                ])

        metric_storage = MetricStorage(sections=[ReportSection(name='main_sv_section', metrics=ms)])
        report = PerRegionSampleReport(sample=svs_by_experiment.keys()[0].sample,
                                       metric_storage=metric_storage)

        # Writing records
        svs_by_key_by_experiment = OrderedDefaultDict(OrderedDict)
        for e, svs in svs_by_experiment.items():
            for sv in svs:
                svs_by_key_by_experiment[sv.get_key()][e] = sv

        for sv_key, sv_by_experiment in svs_by_key_by_experiment.items():
            sv = next((s for s in sv_by_experiment.values() if s is not None), None)

            row = report.add_row()

            type_dict = {
                'BND': 'Fusion',
                'INV': 'Invertion',
                'INS': 'Insertion',
                'DEL': 'Deletion',
                'DUP': 'Duplication',
            }
            row.add_record('Chr', sv.chrom)
            row.add_record('Type', type_dict.get(sv.type, sv.type) if sv else None)
            row.add_record('Location', Metric.format_value(sv.start) + (('..' + Metric.format_value(sv.end)) if sv.end else '') if sv else None, num=sv.start if sv else None)
            row.add_record('Genes', ', '.join('/'.join(a.genes) for a in sv.annotations))
            row.add_record('Status', 'known' if any(a.known for a in sv.annotations) else None)

            if len(sv_by_experiment.values()) == 1:
                row.add_record('Split read support', sv.split_read_support if sv else None)
                row.add_record('Paired read support', sv.paired_end_support if sv else None)
            else:
                for e, m in sv_by_experiment.items():
                    row.add_record(e.key + ' Split read support', sv.split_read_support if sv else None)
                    row.add_record(e.key + ' Paired read support', sv.paired_end_support if sv else None)

            if any(a.known for a in sv.annotations):
                row.highlighted = True

            if len(sv_by_experiment.keys()) == len(svs_by_experiment.keys()):
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
                for se in gene.seq2c_events:
                    d['events'].append(dict(
                        x=chr_cum_len_by_chrom[gene.chrom] + gene.start + (gene.end - gene.start) / 2,
                        geneName=gene.name,
                        logRatio=se.ab_log2r if se.ab_log2r is not None else se.log2r,
                        ampDel=se.amp_del,
                        fragment=se.fragment))

                    # if not gene.seq2c_event.ab_log2r or gene.seq2c_event.fragment == 'BP':  # breakpoint, meaning part of exon is not amplified

            d['maxY'] = max([e['logRatio'] for e in d['events']] + [2])  # max(chain(data['nrm']['ys'], data['amp']['ys'], data['del']['ys'], [2]))
            d['minY'] = min([e['logRatio'] for e in d['events']] + [-2])  # min(chain(data['nrm']['ys'], data['amp']['ys'], data['del']['ys'], [-2]))

            data[k.lower()] = d

        if len(experiment_by_key.keys()) == 1:
            return json.dumps(data.values()[0])
        else:
            return json.dumps(data)

    def make_seq2c_report(self, seq2c_by_experiment):
        ms = [
            Metric('Gene', with_heatmap=False, max_width=50, align='left', sort_direction='ascending'),  # Gene
            Metric('Chr', with_heatmap=False, max_width=20, align='right'),
            Metric('Log ratio', max_width=80),
            Metric('Amp/Del', max_width=70),
            Metric('BP/Whole', max_width=50),
        ]

        reports = []

        # Writing records
        for e, seq2c in seq2c_by_experiment.items():
            metric_storage = MetricStorage(sections=[ReportSection(name='seq2c_section', metrics=ms)])
            report = PerRegionSampleReport(sample=seq2c_by_experiment.keys()[0].sample,
                                           metric_storage=metric_storage)
            seq2c_events = seq2c.values()
            for event in sorted(seq2c_events, key=lambda e: e.gene.name):
                if event.is_amp() or event.is_del():
                    row = report.add_row()
                    row.add_record('Gene', event.gene.name)
                    row.add_record('Chr', event.gene.chrom.replace('chr', ''))
                    row.add_record('Log ratio', '%.2f' % (event.ab_log2r if event.ab_log2r is not None else event.log2r))
                    row.add_record('Amp/Del', event.amp_del)
                    row.add_record('BP/Whole', event.fragment)
            reports.append(report)

        if len(seq2c_by_experiment.keys()) == 1:
            return reports[0]
        else:
            return reports


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
                for se in gene.seq2c_events:
                    if se and (se.is_amp() or se.is_del()):
                        mut_info_by_gene[gene.name].append(se.amp_del + ', ' + se.fragment)

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
        # d['panel'] = 'unknown'
        # d['bed_path'] = 'unknown'
        # d['target_type'] = 'unknown'
        # d['target_fraction'] = 'unknown'
        # d['ave_depth'] = 'unknown'

        d['key'] = experiment.key
        d['sample'] = experiment.sample.name.replace('_', ' ')
        if experiment.patient and experiment.patient.gender:
            d['patient'] = {'sex': experiment.patient.gender}
        d['project_name'] = experiment.project_name.replace('_', ' ')
        if experiment.project_report_path:
            d['project_report_rel_path'] = relpath(experiment.project_report_path, dirname(experiment.sample.clinical_html))
        if experiment.target:
            d['target_section'] = dict()
            d['target_section']['panel'] = experiment.target.type
            d['target_section']['bed_path'] = experiment.target.bed_fpath
            d['target_section']['targqc_link'] = experiment.target.targqc_link
            d['target_section']['target_type'] = experiment.target.type
            d['target_section']['target_fraction'] = Metric.format_value(experiment.target.coverage_percent, is_html=True, unit='%')
            # if self.cnf.debug:
            #     d['target_section']['panel'] = experiment.target.type + ', AZ300 IDT panel'
            #     d['target_section']['bed_path'] = 'http://blue.usbod.astrazeneca.net/~klpf990/reference_data/genomes/Hsapiens/hg19/bed/Panel-IDT_PanCancer_AZSpike_V1.bed'

        d['sample_type'] = experiment.sample.normal_match
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
        self.sv_report = None
        self.actionable_genes_report = None
        self.seq2c_plot_data = None
        self.seq2c_report = None
        self.key_genes_report = None
        self.cov_plot_data = None

        info('Preparing data...')
        if self.experiment.mutations:
            self.mutations_report = self.make_mutations_report({self.experiment: self.experiment.mutations})
            self.mutations_plot_data = self.make_mutations_json({self.experiment: self.experiment.mutations})
            self.substitutions_plot_data = self.make_substitutions_json({self.experiment: self.experiment.mutations})
        if self.experiment.sv_events:
            self.sv_report = self.make_sv_report({self.experiment: self.experiment.sv_events})
        if self.experiment.seq2c_events_by_gene_name:
            self.seq2c_plot_data = self.make_seq2c_plot_json({self.experiment.key: self.experiment})
            self.seq2c_report = self.make_seq2c_report({self.experiment: self.experiment.seq2c_events_by_gene_name})
        if self.experiment.actionable_genes_dict and \
                (self.experiment.mutations or self.experiment.seq2c_events_by_gene_name or self.experiment.sv_events):
            self.actionable_genes_report = self.make_actionable_genes_report(self.experiment.actionable_genes_dict)
        if self.experiment.ave_depth:
            self.key_genes_report = self.make_key_genes_cov_report(self.experiment.key_gene_by_name, self.experiment.ave_depth)
            self.cov_plot_data = self.make_key_genes_cov_json({self.experiment.key: self.experiment})

    def write_report(self, output_fpath):
        info('')

        data = {
            'key_or_target': self.experiment.genes_collection_type,
            'genes_description': self.experiment.genes_description,
            'sample': self.sample_section(self.experiment),
            'variants': self.__mutations_section(),
            'coverage': self.__coverage_section(),
            'actionable_genes': self.__actionable_genes_section()
        }
        if self.sv_report:
            data['sv'] = {}
            section = self.__sv_section()
            if section:
                data['sv'] = {'report': section}
        if self.seq2c_plot_data:
            data['seq2c'] = {'plot_data': self.seq2c_plot_data}
            if self.seq2c_report:
                data['seq2c']['amp_del'] = self.__seq2c_section()
        write_static_html_report(self.cnf, data, output_fpath,
           tmpl_fpath=join(dirname(abspath(__file__)), 'template.html'),
           extra_js_fpaths=[join(dirname(abspath(__file__)), 'static', 'clinical_report.js'),
                            join(dirname(abspath(__file__)), 'static', 'draw_mutations_plot.js'),
                            join(dirname(abspath(__file__)), 'static', 'draw_substitutions_plot.js'),
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
            mutations_dict['plot_data'] = self.mutations_plot_data
            mutations_dict['substitutions_plot_data'] = self.substitutions_plot_data
        return mutations_dict

    def __sv_section(self):
        sv_dict = dict()
        if self.sv_report and self.sv_report.rows:
            sv_dict['table'] = build_report_html(self.sv_report, sortable=True)
            sv_dict['total_key_genes'] = Metric.format_value(len(self.experiment.key_gene_by_name), is_html=True)
        return sv_dict

    def __coverage_section(self):
        if self.experiment.depth_cutoff is not None and self.key_genes_report is not None:
            coverage_dict = dict(depth_cutoff=self.experiment.depth_cutoff, columns=[])
            GENE_COL_NUM = 3
            genes_in_col = [len(self.key_genes_report.rows) / GENE_COL_NUM] * GENE_COL_NUM
            for i in range(len(self.key_genes_report.rows) % GENE_COL_NUM):
                genes_in_col[i] += 1
            calc_cell_contents(self.key_genes_report, self.key_genes_report.get_rows_of_records())
            printed_genes = 0
            for i in range(GENE_COL_NUM):
                column_dict = dict()
                # column_dict['table'] = build_report_html(coverage_report)
                column_dict['metric_names'] = [make_cell_th(m) for m in self.key_genes_report.metric_storage.get_metrics()]
                column_dict['rows'] = [
                    dict(records=[make_cell_td(r) for r in region.records])
                        for region in self.key_genes_report.rows[printed_genes:printed_genes + genes_in_col[i]]]
                coverage_dict['columns'].append(column_dict)
                printed_genes += genes_in_col[i]
            coverage_dict['plot_data'] = self.cov_plot_data
            return coverage_dict
        else:
            return dict()

    def __seq2c_section(self):
        seq2c_dict = dict()
        if self.seq2c_report and self.seq2c_report.rows:
            seq2c_dict = dict(columns=[])
            GENE_COL_NUM = min(3, len(self.seq2c_report.rows))
            genes_in_col = [len(self.seq2c_report.rows) / GENE_COL_NUM] * GENE_COL_NUM
            for i in range(len(self.seq2c_report.rows) % GENE_COL_NUM):
                genes_in_col[i] += 1
            calc_cell_contents(self.seq2c_report, self.seq2c_report.get_rows_of_records())
            printed_genes = 0
            for i in range(GENE_COL_NUM):
                column_dict = dict()
                column_dict['metric_names'] = [make_cell_th(m) for m in self.seq2c_report.metric_storage.get_metrics()]
                column_dict['rows'] = [
                    dict(records=[make_cell_td(r) for r in region.records])
                        for region in self.seq2c_report.rows[printed_genes:printed_genes + genes_in_col[i]]]
                seq2c_dict['columns'].append(column_dict)
                printed_genes += genes_in_col[i]
        return seq2c_dict

    def __actionable_genes_section(self):
        actionable_genes_dict = dict()
        if self.actionable_genes_report and self.actionable_genes_report.rows:
            actionable_genes_dict['table'] = build_report_html(self.actionable_genes_report, sortable=False)
        return actionable_genes_dict

    def make_key_genes_cov_report(self, key_gene_by_name, ave_depth):
        if self.experiment.depth_cutoff is None:
            return None

        info('Making ' + self.experiment.genes_collection_type + ' genes coverage report...')
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
            if any(se.is_amp() or se.is_del() for se in gene.seq2c_events):
                reg.add_record('CNV', '; '.join([se.amp_del + ', ' + se.fragment for se in gene.seq2c_events if se.is_amp() or se.is_del()]))

        return key_genes_report

    def make_actionable_genes_report(self, actionable_genes_dict):
        info('Preparing mutations stats for ' + self.experiment.genes_collection_type + ' gene tables')

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

        sv_mutation_types = {'Rearrangement', 'Fusion', 'Amplification', 'Deletion'}
        cnv_mutation_types = {'Amplification', 'Deletion'}

        for gene in self.experiment.key_gene_by_name.values():
            if gene.name not in actionable_gene_names:
                continue
            possible_mutation_types = set(actionable_genes_dict[gene.name][1].split('; '))
            # possible_mutation_types = possible_mutation_types - sv_mutation_types
            # if not possible_mutation_types: continue

            variants = []
            types = []

            if self.experiment.mutations:
                vardict_mut_types = possible_mutation_types - sv_mutation_types - cnv_mutation_types
                if vardict_mut_types:
                    for mut in self.experiment.mutations:
                        if mut.gene.name == gene.name:
                            variants.append(mut.aa_change if mut.aa_change else '.')
                            types.append(mut.var_type)

            if cnv_mutation_types:
                for se in gene.seq2c_events:
                    if 'Amplification' in possible_mutation_types and se.amp_del == 'Amp' or \
                            'Deletion' in possible_mutation_types and se.amp_del == 'Del':
                        variants.append(se.amp_del + ', ' + se.fragment)
                        types.append(se.amp_del)

            if sv_mutation_types:
                for se in gene.sv_events:
                    if ('Fusion' in possible_mutation_types or 'Rearrangement' in possible_mutation_types) and se.type == 'BND' or \
                       'Deletion' in possible_mutation_types and se.type == 'DEL' or \
                       'Amplification' in possible_mutation_types and se.type == 'DUP':
                        variants.append(se.type)
                        types.append(se.type)

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
