from os.path import join
from source import info
from source.clinical_reporting.clinical_parser import clinical_sample_info_from_bcbio_structure
from source.clinical_reporting.clinical_reporting import Chromosome
from source.reporting.reporting import MetricStorage, Metric, PerRegionSampleReport
from source.reporting.reporting import ReportSection


def run_clinical_target2wgs(cnf, wgs_bs, target_bs, shared_sample_names, output_dirpath):
    info('Running clinical reporting comparison')

    for sname in shared_sample_names:
        info('Preparing ' + sname + '...')
        target_sample = next(s for s in target_bs.samples if s.name == sname)
        wgs_sample = next(s for s in wgs_bs.samples if s.name == sname)

        info('-' * 70)
        clinsample_target_info = clinical_sample_info_from_bcbio_structure(cnf, target_bs, target_sample)
        info('')
        info('-' * 70)
        clinsample_wgs_info = clinical_sample_info_from_bcbio_structure(cnf, wgs_bs, wgs_sample)

        info('')
        info('*' * 70)
        run_sample_clinreport_target2wgs(cnf, clinsample_target_info, clinsample_wgs_info, output_dirpath)
        info('*' * 70)
        info('Successfully finished.')


def run_sample_clinreport_target2wgs(cnf, clinsample_target_info, clinsample_wgs_info, output_dirpath):
    ComparisonClinicalReporting(cnf, clinsample_target_info, clinsample_wgs_info).write_report(
        join(output_dirpath, 'clinical_report.html'))


class ComparisonClinicalReporting:
    def __init__(self, cnf, target_info, wgs_info):
        self.target_info = target_info
        self.wgs_info = wgs_info
        self.target_sample = target_info.sample
        self.wgs_sample = wgs_info.sample
        self.cnf = cnf
        self.chromosomes_by_name = Chromosome.build_chr_section(self.cnf)

        info('Preparing data...')
        self.mutations_report = self.make_mutations_report()
        self.actionable_genes_report = self.make_actionable_genes_report(self.info.actionable_genes_dict)
        self.seq2c_plot_data = self.make_seq2c_plot_json() if self.info.seq2c_events_by_gene_name is not None else None
        self.key_genes_report = self.make_key_genes_cov_report(self.info.key_gene_by_name, self.info.ave_depth)
        self.cov_plot_data = self.make_key_genes_cov_json(self.info.key_gene_by_name)

    def write_report(self, output_fpath):
        info('')

        write_static_html_report(self.cnf, {
            'sample': self.__sample_section(),
            'variants': self.__mutations_section(),
            'seq2c': self.__seq2c_section(),
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

    def __sample_section(self):
        sample_dict = dict()
        sample_dict['sample'] = self.sample.name.replace('_', ' ')
        if self.info.patient.gender:
            sample_dict['patient'] = {'sex': self.info.patient.gender}
        sample_dict['project_name'] = self.info.project_name.replace('_', ' ')
        if self.info.project_report_path:
            sample_dict['project_report_rel_path'] = relpath(self.info.project_report_path, dirname(self.sample.clinical_html))
        sample_dict['panel'] = self.info.target.type
        sample_dict['bed_path'] = self.info.target.bed_fpath or ''
        if self.cnf.debug:
            sample_dict['panel'] = self.info.target.type + ', AZ300 IDT panel'
            sample_dict['bed_path'] = 'http://blue.usbod.astrazeneca.net/~klpf990/reference_data/genomes/Hsapiens/hg19/bed/Panel-IDT_PanCancer_AZSpike_V1.bed'

        sample_dict['sample_type'] = self.sample.normal_match if self.sample.normal_match else 'unpaired'  # plasma, unpaired'
        sample_dict['genome_build'] = self.cnf.genome.name
        sample_dict['target_type'] = self.info.target.type
        sample_dict['target_fraction'] = Metric.format_value(self.info.target.coverage_percent, is_html=True, unit='%')
        # approach_dict['min_depth'] = Metric.format_value(min_depth, is_html=True)
        sample_dict['ave_depth'] = Metric.format_value(self.info.ave_depth, is_html=True)
        return sample_dict

    def __mutations_section(self):
        mutations_dict = dict()
        if self.mutations_report.rows:
            # if cnf.debug:
            #     mutations_report.regions = mutations_report.regions[::20]
            mutations_dict['table'] = build_report_html(self.mutations_report, sortable=True)
            mutations_dict['total_variants'] = Metric.format_value(self.info.total_variants, is_html=True)
            mutations_dict['total_key_genes'] = Metric.format_value(len(self.info.key_gene_by_name), is_html=True)
        return mutations_dict

    def __coverage_section(self):
        coverage_dict = dict(depth_cutoff=self.info.depth_cutoff, columns=[])
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

    def __seq2c_section(self):
        seq2c_dict = dict()
        if self.seq2c_plot_data:
            seq2c_dict['plot_data'] = self.seq2c_plot_data
        return seq2c_dict

    def __actionable_genes_section(self):
        actionable_genes_dict = dict()
        if self.actionable_genes_report.rows:
            actionable_genes_dict['table'] = build_report_html(self.actionable_genes_report, sortable=False)
        return actionable_genes_dict

    def make_key_genes_cov_report(self, key_gene_by_name, ave_depth):
        info('Making key genes coverage report...')
        clinical_cov_metrics = [
            Metric('Gene'),
            Metric('Chr', with_heatmap=False, max_width=20, align='right'),
            Metric('Ave depth', med=ave_depth),
            Metric('% cov at {}x'.format(self.info.depth_cutoff), unit='%', med=1, low_inner_fence=0.5, low_outer_fence=0.1),
            Metric('CNV', short_name='&nbsp;&nbsp;CNV')]  # short name is hack for IE9 who doesn't have "text-align: left" and tries to stick "CNV" to the previous col header

        clinical_cov_metric_storage = MetricStorage(sections=[ReportSection(metrics=clinical_cov_metrics)])

        key_genes_report = PerRegionSampleReport(sample=self.sample, metric_storage=clinical_cov_metric_storage)

        for gene in sorted(key_gene_by_name.values(), key=lambda g: g.name):
            reg = key_genes_report.add_row()
            reg.add_record('Gene', gene.name)
            reg.add_record('Chr', gene.chrom.replace('chr', ''))
            reg.add_record('Ave depth', gene.ave_depth)
            m = clinical_cov_metric_storage.find_metric('% cov at {}x'.format(self.info.depth_cutoff))
            reg.add_record(m.name, next((cov for cutoff, cov in gene.cov_by_threshs.items() if cutoff == self.info.depth_cutoff), None))
            if gene.seq2c_event and (gene.seq2c_event.is_amp() or gene.seq2c_event.is_del()):
                reg.add_record('CNV', gene.seq2c_event.amp_del + ', ' + gene.seq2c_event.fragment)

        return key_genes_report

    @staticmethod
    def make_db_html(mut):
        db = ''
        if mut.dbsnp_ids:
            db += ', '.join(('<a href="http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?searchType=adhoc_search'
                             '&type=rs&rs=' + rs_id + '">dbSNP</a>') for rs_id in mut.dbsnp_ids)
        if db and mut.cosmic_ids:
            db += ', '
        if mut.cosmic_ids:
            db += ', '.join('<a href="http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=' +
                            cid + '">COSM</a>' for cid in mut.cosmic_ids)
        return db

    def make_mutations_report(self):
        clinical_mut_metric_storage = MetricStorage(
            sections=[ReportSection(metrics=[
                Metric('Gene'),  # Gene & Transcript
                Metric('Transcript'),  # Gene & Transcript
                # Metric('Codon chg', max_width=max_width, class_='long_line'),            # c.244G>A
                Metric('AA chg', max_width=70, class_='long_line'),            # p.Glu82Lys
                # Metric('Allele'),             # Het.
                # Metric('Chr', max_width=33, with_heatmap=False),       # chr11
                Metric('Position',
                       # sort_by=lambda v: (v.split(':')[0], int(''.join(ch for ch in v.split(':')[1] if ch.isdigit()))),
                       with_heatmap=False, align='left', sort_direction='ascending'),       # g.47364249
                Metric('Change', max_width=95, class_='long_line'),       # G>A
                Metric('DP Target', max_width=48),              # 658
                Metric('DP WGS', max_width=48),              # 658
                Metric('AF Target', max_width=45, unit='%', with_heatmap=False),          # .19
                Metric('AF WGS', max_width=45, unit='%', with_heatmap=False),          # .19
                Metric('AA len', max_width=50, with_heatmap=False),          # 128
                # Metric('COSMIC', max_width=70, style='', class_='long_line'),                 # rs352343, COSM2123
                Metric('Effect', max_width=100, class_='long_line'),               # Frameshift
                Metric('VarDict status', short_name='Pathogenic,\nreported by VarDict'),     # Likely
                # Metric('VarDict reason', short_name='VarDict\nreason'),     # Likely
                Metric('Databases'),                 # rs352343, COSM2123
                Metric('ClinVar', short_name='SolveBio ClinVar'),    # Pathogenic?, URL
            ])])

        report = PerRegionSampleReport(sample=self.target_info.sample, metric_storage=clinical_mut_metric_storage, expandable=True)

        for mut in mutations:
            row = report.add_row()
            row.add_record('Gene', mut.gene.name)
            row.add_record('Transcript', mut.transcript)
            # row.add_record('Codon chg', mut.codon_change)

            aa_chg = ''.join([gray(c) if c.isdigit() else c for c in (mut.aa_change or '')])
            row.add_record('AA chg', aa_chg)

            c = (mut.chrom.replace('chr', '')) if mut.chrom else ''
            p = Metric.format_value(mut.pos, human_readable=True, is_html=True) if mut.pos else ''
            row.add_record('Position', gray(c + ':') + p,
                           num=mut.get_chrom_key() * 100000000000 + mut.pos)

            chg = mut.ref + '>' + mut.alt
            if mut.var_type:
                t = mut.var_type
                if t in ['Insertion', 'Deletion']:
                    t = t[:3]
                chg = gray(t) + ' ' + chg
            row.add_record('Change', chg)

            row.add_record('Depth', mut.depth)
            row.add_record('Freq', mut.freq)
            row.add_record('AA len', mut.aa_len)
            row.add_record('Effect', mut.eff_type)

            status = mut.status
            if mut.reason:
                status += gray(' (' + mut.reason + ')')
            row.add_record('VarDict status', status)
            # row.add_record('VarDict reason', mut.reason)
            row.add_record('Databases', self.make_db_html(mut), parse=False)
            row.add_record('ClinVar', mut.solvebio.clinsig if mut.solvebio else '', url=mut.solvebio.url if mut.solvebio else None)

            # Highlighting or hiding:
            if not mut.status or mut.status.lower() == 'unknown':
                if mut.solvebio and 'Pathogenic' in mut.solvebio.clinsig:
                    warn('Mutation ' + str(mut) + ' is unknown, but found in SolveBio')
                row.hidden = True
            else:
                if mut.status and mut.status.lower() == 'known':
                    row.highlighted = True
                # if mut.solvebio and 'Pathogenic' in mut.solvebio.clinsig:
                #     row.highlighted = True

        return report

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

        for gene in self.info.key_gene_by_name.values():
            if gene.name not in actionable_gene_names:
                continue
            possible_mutation_types = set(actionable_genes_dict[gene.name][1].split('; '))
            # possible_mutation_types = possible_mutation_types - sv_mutation_types
            # if not possible_mutation_types: continue

            variants = []
            types = []

            vardict_mut_types = possible_mutation_types - sv_mutation_types - cnv_mutation_types
            if vardict_mut_types:
                for mut in self.info.mutations:
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

    def make_seq2c_plot_json(self):
        chr_cum_lens = Chromosome.get_cum_lengths(self.chromosomes_by_name)
        chr_cum_len_by_chrom = dict(zip([c.name for c in self.chromosomes_by_name.values()], chr_cum_lens))

        data = dict(
            events=[],
            ticksX=[[(chr_cum_lens[i] + chr_cum_lens[i + 1])/2, self.chromosomes_by_name.values()[i].short_name]
                    for i in range(len(self.chromosomes_by_name.keys()))],
            linesX=chr_cum_lens
        )

        for gene in self.info.key_gene_by_name.values():
            if gene.seq2c_event:
                data['events'].append(dict(
                    x=chr_cum_len_by_chrom[gene.chrom] + gene.start + (gene.end - gene.start) / 2,
                    geneName=gene.name,
                    logRatio=gene.seq2c_event.ab_log2r if gene.seq2c_event.ab_log2r is not None else gene.seq2c_event.log2r,
                    ampDel=gene.seq2c_event.amp_del,
                    fragment=gene.seq2c_event.fragment))

                # if not gene.seq2c_event.ab_log2r or gene.seq2c_event.fragment == 'BP':  # breakpoint, meaning part of exon is not amplified

        data['maxY'] = max([e['logRatio'] for e in data['events']] + [2])  # max(chain(data['nrm']['ys'], data['amp']['ys'], data['del']['ys'], [2]))
        data['minY'] = min([e['logRatio'] for e in data['events']] + [-2])  # min(chain(data['nrm']['ys'], data['amp']['ys'], data['del']['ys'], [-2]))

        return json.dumps(data)


    def make_key_genes_cov_json(self, key_gene_by_name):
        chr_cum_lens = Chromosome.get_cum_lengths(self.chromosomes_by_name)
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
            cov_in_thresh.append(gene.cov_by_threshs.get(self.info.depth_cutoff))
            coord_x.append(chr_cum_len_by_chrom[gene.chrom] + gene.start + (gene.end - gene.start) / 2)
            cds_cov_by_gene[gene.name] = [dict(
                start=cds.start,
                end=cds.end,
                aveDepth=cds.ave_depth,
                percentInThreshold=cds.cov_by_threshs.get(self.info.depth_cutoff)
            ) for cds in gene.cdss]

        json_txt = '{"coords_x":%s, "ave_depths":%s, "cov_in_thresh":%s, "gene_names":%s, "mutations":%s, "ticks_x":%s, ' \
                    '"lines_x":%s, "cds_cov_by_gene":%s}'\
                   % (json.dumps(coord_x), json.dumps(gene_ave_depths), json.dumps(cov_in_thresh), json.dumps(gene_names),
                      json.dumps(mut_info_by_gene), json.dumps(ticks_x), json.dumps(chr_cum_lens), json.dumps(cds_cov_by_gene))
        return json_txt
