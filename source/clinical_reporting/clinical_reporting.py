from collections import OrderedDict, defaultdict
import json
from genericpath import exists
from os.path import join, dirname, abspath, relpath, basename

import re

import source
from ext_modules.genologics import lims
from source import info, verify_file
from source.calling_process import call
from source.logger import warn
from source.reporting.reporting import MetricStorage, Metric, PerRegionSampleReport, ReportSection, calc_cell_contents, make_cell_td, write_static_html_report, make_cell_th, build_report_html
from source.tools_from_cnf import get_script_cmdline
from source.utils import get_chr_lengths, OrderedDefaultDict
from tools.add_jbrowse_tracks import get_jbrowser_link


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
            (chr_name, Chromosome(chr_name, length=l)) for chr_name, l in get_chr_lengths(cnf)
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

    def make_mutations_report(self, mutations_by_experiment, jbrowser_link):
        ms = [
            Metric('Gene'),  # Gene & Transcript
            Metric('Transcript'),  # Gene & Transcript
            # Metric('Codon chg', max_width=max_width, class_='long_line'),            # c.244G>A
            Metric('AA len', max_width=50, class_='stick_to_left', with_heatmap=False),          # 128
            Metric('AA chg', short_name='AA change', max_width=70, class_='long_line'),            # p.Glu82Lys
            # Metric('Allele'),             # Het.
            # Metric('Chr', max_width=33, with_heatmap=False),       # chr11
            Metric('Position', with_heatmap=False, align='left', sort_direction='ascending'),       # g.47364249
            Metric('Change', max_width=100, class_='long_line', description='Genomic change'),       # G>A
            Metric('cDNA change', class_='long_line', description='cDNA change'),       # G>A
            # Metric('COSMIC', max_width=70, style='', class_='long_line'),                 # rs352343, COSM2123
            Metric('Status', short_name='Status'),     # Somatic
            Metric('Effect', max_width=150, class_='long_line'),               # Frameshift
            Metric('VarDict status', short_name='Status', max_width=230, class_='long_line'),     # Likely
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

        clinical_mut_metric_storage = MetricStorage(sections=[ReportSection(metrics=ms, name='mutations')])
        report = PerRegionSampleReport(sample=mutations_by_experiment.keys()[0].sample,
            metric_storage=clinical_mut_metric_storage, expandable=True)

        # Writing records
        print_cdna = False
        muts_by_key_by_experiment = OrderedDefaultDict(OrderedDict)
        for e, muts in mutations_by_experiment.items():
            for mut in muts:
                muts_by_key_by_experiment[mut.get_key()][e] = mut
                if mut.cdna_change.strip():
                    print_cdna = True

        mut_canonical = [[m.is_canonical if m is not None else False for m in muts] for e, muts in mutations_by_experiment.items()]
        mut_positions = [m.pos for i, (e, muts) in enumerate(mutations_by_experiment.items()) for j, m in enumerate(muts) if m is not None and mut_canonical[i][j]]

        for mut_key, mut_by_experiment in muts_by_key_by_experiment.items():
            mut = next((m for m in mut_by_experiment.values() if m is not None), None)
            if mut.pos not in mut_positions:
                mut_positions.append(mut.pos)
                row = report.add_row()
                row.add_record('Gene', mut.gene.name)
                row_class = ' expandable_gene_row collapsed'
                row.add_record('Position',
                    **self._pos_recargs(mut.chrom, mut.get_chrom_key(), mut.pos, mut.pos, jbrowser_link))
                row.add_record('Change', **self._g_chg_recargs(mut))

                if len(mutations_by_experiment.values()) == 1:
                    row.add_record('Freq', mut.freq if mut else None)
                    row.add_record('Depth', mut.depth if mut else None)
                else:
                    for e, m in mut_by_experiment.items():
                        row.add_record(e.key + ' Freq', m.freq if m else None)
                        row.add_record(e.key + ' Depth', m.depth if m else None)
                row.class_ = row_class
                self._highlighting_and_hiding_mut_row(row, mut)
                if len(mut_by_experiment.keys()) == len(mutations_by_experiment.keys()):
                    row.highlighted_green = True

            row = report.add_row()
            row.add_record('Gene', mut.gene.name, show_content=mut.is_canonical)
            if mut.is_canonical:
                row.add_record('Transcript', mut.transcript)
                row_class = ' expandable_gene_row collapsed'
            else:
                row.add_record('Transcript', mut.transcript)
                row_class = ' row_to_hide row_hidden'
            row.add_record('AA chg', **self._aa_chg_recargs(mut))
            row.add_record('Position', show_content=mut.is_canonical,
               **self._pos_recargs(mut.chrom, mut.get_chrom_key(), mut.pos, mut.pos, jbrowser_link))
            row.add_record('Change', show_content=mut.is_canonical, **self._g_chg_recargs(mut))
            if print_cdna:
                row.add_record('cDNA change', **self._cdna_chg_recargs(mut))
            row.add_record('AA len', mut.aa_len)
            if mut.status:
                row.add_record('Status', mut.status)
            row.add_record('Effect', mut.eff_type)
            row.add_record('VarDict status', **self._signif_field(mut))
            # row.add_record('VarDict reason', mut.reason)
            row.add_record('Databases', **self._db_recargs(mut))
            # if is_local():
            #     row.add_record('ClinVar', **self._clinvar_recargs(mut))

            if len(mutations_by_experiment.values()) == 1:
                row.add_record('Freq', mut.freq if mut else None, show_content=mut.is_canonical)
                row.add_record('Depth', mut.depth if mut else None, show_content=mut.is_canonical)
            else:
                for e, m in mut_by_experiment.items():
                    row.add_record(e.key + ' Freq', m.freq if m else None, show_content=mut.is_canonical)
                    row.add_record(e.key + ' Depth', m.depth if m else None, show_content=mut.is_canonical)

            row.class_ = row_class
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
                if mut.is_canonical:
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
            substitutions_sum = sum([d['substitutions'][s] for s in d['substitutions']])
            d['maxRate'] = d['maxY'] * 100 / substitutions_sum if substitutions_sum > 0 else 0
            d['minY'] = 0

            data[e.key.lower()] = d

        if len(mutations_by_experiment.keys()) == 1:
            return json.dumps(data.values()[0])
        else:
            return json.dumps(data)

    sv_type_dict = {
        'BND': 'Fusion',
        'INV': 'Invertion',
        'INS': 'Insertion',
        'DEL': 'Deletion',
        'DUP': 'Duplication',
    }

    def make_sv_report(self, svs_by_experiment, jbrowser_link):
        ms = [
            Metric('Genes'),
            # Metric('Chr', with_heatmap=False, max_width=50, align='left', sort_direction='ascending'),
            Metric('Type'),
            Metric('Location', with_heatmap=False, align='left', sort_direction='ascending', style="width: 500px"),
            Metric('Transcript', align='left', style="width: 300px"),
            # Metric('Status', with_heatmap=False, align='left', sort_direction='ascending'),
            # Metric('Effects', with_heatmap=False, align='left', class_='long_line'),
        ]

        if len(svs_by_experiment) == 1:
            ms.extend([
                Metric('Split read support', short_name='Split /'),
                Metric('Paired read support', short_name='paired read support'),
            ])
        else:
            for e in svs_by_experiment.keys():
                ms.extend([
                    Metric(e.key + ' Split read support', short_name='Split /'),
                    Metric(e.key + ' Paired read support', short_name='paired read support'),
                ])

        metric_storage = MetricStorage(sections=[ReportSection(name='main_sv_section', metrics=ms)])
        report = PerRegionSampleReport(sample=svs_by_experiment.keys()[0].sample,
                                       metric_storage=metric_storage)

        # Writing records
        svanns_by_key_by_experiment = OrderedDefaultDict(OrderedDict)
        for e, svs in svs_by_experiment.items():
            exon_dels_cnt = 0
            known_cnt = 0
            affecting_2_genes_cnt = 0
            other_cnt = 0
            for sv_event in svs:
                for an in sv_event.key_annotations:
                    # reporting all whole exon deletions
                    if sv_event.is_deletion() and ('exon_del' in an.effect.lower() or 'exon_loss' in an.effect.lower()):
                        svanns_by_key_by_experiment[an.get_key()][e] = an
                        exon_dels_cnt += 1

                    # reporting all known (fusion)
                    elif an.known:
                        svanns_by_key_by_experiment[an.get_key()][e] = an
                        known_cnt += 1

                    # # reporting all non-fusion events affecting 2 or more genes (start and end should not be the same gene. handling overlapping gene cases.)
                    # elif sv_event.end_genes and all(ann_g not in sv_event.end_genes for ann_g in an.genes):
                    #     svanns_by_key_by_experiment[an.get_key()][e] = an
                    #     affecting_2_genes_cnt += 1

                    else:
                        other_cnt += 1

            info('exon_dels_cnt: ' + str(exon_dels_cnt))
            info('known_cnt: ' + str(known_cnt))
            info('affecting_2_genes_cnt: ' + str(affecting_2_genes_cnt))
            info('other_cnt: ' + str(other_cnt))

        for sv_key, svann_by_experiment in svanns_by_key_by_experiment.items():
            sv_ann = next((s for s in svann_by_experiment.values() if s is not None), None)

            row = report.add_row()

            row.add_record('Genes', '/'.join(set(sv_ann.genes)) if sv_ann else None)

            type_str = None
            if sv_ann:
                type_str = BaseClinicalReporting.sv_type_dict.get(sv_ann.event.type, sv_ann.event.type)
                if sv_ann.effect == 'EXON_DEL':
                    type_str += ' ' + sv_ann.exon_info
                if type_str == 'Fusion':
                    type_str = 'Known fusion'

            row.add_record('Type', type_str)

            if sv_ann and sv_ann.event.start:
                row.add_record('Location',
                    **self._pos_recargs(sv_ann.event.chrom, sv_ann.event.get_chrom_key(),
                                        sv_ann.event.start, sv_ann.event.end, jbrowser_link))
            else:
                row.add_record('Location', None)

            row.add_record('Transcript', sv_ann.transcript if sv_ann else None)
            # row.add_record('Status', 'known' if any(a.known for a in sv.annotations) else None)
            # row.add_record('Effects', ', '.join(set(a.effect.lower().replace('_', ' ') for a  in sv.annotations if a in ['EXON_DEL', 'FUSION'])))

            if len(svann_by_experiment.values()) == 1:
                row.add_record('Split read support', sv_ann.event.split_read_support if sv_ann else None)
                row.add_record('Paired read support', sv_ann.event.paired_end_support if sv_ann else None)
            else:
                for _sv_ann, m in svann_by_experiment.items():
                    row.add_record(e.key + ' Split read support', _sv_ann.event.split_read_support if _sv_ann else None)
                    row.add_record(e.key + ' Paired read support', _sv_ann.event.paired_end_support if _sv_ann else None)

            if any(an.known for an in svann_by_experiment.values()):
                row.highlighted = True

            if len(svann_by_experiment.keys()) == len(svs_by_experiment.keys()):
                row.highlighted_green = True

        return report

    def make_seq2c_plot_json(self, experiment_by_key):
        data = dict()

        for k, e in experiment_by_key.items():
            # if len(e.seq2c_events_by_gene.values()) == 0:
            #     data[k.lower()] = None
            #     continue

            chr_cum_lens = Chromosome.get_cum_lengths(self.chromosomes_by_name)
            chr_cum_len_by_chrom = dict(zip([c.name for c in self.chromosomes_by_name.values()], chr_cum_lens))

            d = dict(
                events=[],
                ticksX=[[(chr_cum_lens[i] + chr_cum_lens[i + 1])/2, self.chromosomes_by_name.values()[i].short_name]
                        for i in range(len(self.chromosomes_by_name.keys()))],
                linesX=chr_cum_lens
            )

            for gene, se in e.seq2c_events_by_gene.items():
                if gene.chrom not in chr_cum_len_by_chrom:
                   warn('Gene ' + gene.name + ' chromosome ' + gene.chrom + ' not found')
                else:
                    d['events'].append(dict(
                        x=chr_cum_len_by_chrom[gene.chrom] + gene.start + (gene.end - gene.start) / 2,
                        geneName=gene.name,
                        logRatio=se.ab_log2r if se.ab_log2r is not None else se.log2r,
                        ampDel=se.amp_del,
                        fragment=se.fragment,
                        isKeyGene=gene.key in e.key_gene_by_name_chrom))

                    # if not gene.seq2c_event.ab_log2r or gene.seq2c_event.fragment == 'BP':  # breakpoint, meaning part of exon is not amplified

            d['maxY'] = max([e['logRatio'] for e in d['events']] + [2])  # max(chain(data['nrm']['ys'], data['amp']['ys'], data['del']['ys'], [2]))
            d['minY'] = min([e['logRatio'] for e in d['events']] + [-2])  # min(chain(data['nrm']['ys'], data['amp']['ys'], data['del']['ys'], [-2]))

            if all(e['ampDel'] is None for e in d['events']):
                d = None
            data[k.lower()] = d

        if len(experiment_by_key.keys()) == 1:
            return json.dumps(data.values()[0]) if data.values()[0] else None
        else:
            return json.dumps(data) if all(d is not None for d in data.values()) else None

    def make_seq2c_report(self, seq2c_by_experiment):
        ms = [
            Metric('Gene', align='left', sort_direction='ascending'),  # Gene
            Metric('Chr', with_heatmap=False, max_width=20, align='right'),
            Metric('Log ratio'),
            Metric('Amp/Del'),
            Metric('BP/Whole'),
        ]

        reports = []

        # Writing records
        for e, seq2c in seq2c_by_experiment.items():
            report = None
            if len(e.seq2c_events_by_gene.values()) > 0:
                metric_storage = MetricStorage(sections=[ReportSection(name='seq2c_section', metrics=ms)])
                report = PerRegionSampleReport(sample=seq2c_by_experiment.keys()[0].sample,
                                               metric_storage=metric_storage, expandable=True)
                seq2c_events = seq2c.values()
                is_whole_genomic_profile = len(e.seq2c_events_by_gene.values()) > len(e.key_gene_by_name_chrom.values())
                for event in sorted(seq2c_events, key=lambda e: e.gene.name):
                    if event.is_amp() or event.is_del():
                        row = report.add_row()
                        row.add_record('Gene', event.gene.name)
                        row.add_record('Chr', event.gene.chrom.replace('chr', ''))
                        row.add_record('Log ratio', '%.2f' % (event.ab_log2r if event.ab_log2r is not None else event.log2r))
                        row.add_record('Amp/Del', event.amp_del)
                        row.add_record('BP/Whole', event.fragment)
                        if is_whole_genomic_profile and event.gene.key not in e.key_gene_by_name_chrom:
                            row.hidden = True
            reports.append(report)

        if len(seq2c_by_experiment.keys()) == 1:
            return reports[0]
        else:
            return reports

    def make_key_genes_cov_json(self, experiment_by_key):
        chr_cum_lens = Chromosome.get_cum_lengths(self.chromosomes_by_name)
        chr_cum_len_by_chrom = dict(zip([c.name for c in self.chromosomes_by_name.values()], chr_cum_lens))

        gene_names = []
        transcript_names = []
        strands = []
        coord_x = []

        ticks_x = [[(chr_cum_lens[i] + chr_cum_lens[i + 1])/2, self.chromosomes_by_name.values()[i].short_name]
                   for i in range(len(self.chromosomes_by_name.keys()))]

        hits = list()
        for key, e in experiment_by_key.items():
            gene_ave_depths = []
            covs_in_thresh = []
            cds_cov_by_gene = defaultdict(list)
            mut_info_by_gene = dict()

            for gene in e.key_gene_by_name_chrom.values():
                mut_info_by_gene[gene.name] = [('p.' + m.aa_change if m.aa_change else '.') for m in gene.mutations if m.is_canonical]
                for se in gene.seq2c_events:
                    if se and (se.is_amp() or se.is_del()):
                        mut_info_by_gene[gene.name].append(se.amp_del + ', ' + se.fragment)

            for gene in e.key_gene_by_name_chrom.values():
                gene_names.append(gene.name)
                transcript_names.append(gene.transcript_id)
                strands.append(gene.strand)
                gene_ave_depths.append(gene.ave_depth)
                covs_in_thresh.append(gene.cov_by_threshs.get(e.depth_cutoff))
                # if not gene.start or not gene.end:
                #     print gene.name, gene.chrom
                #     continue
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
            transcript_names=transcript_names,
            strands=strands,
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
        d['project_report_rel_path'] = 'not generated'
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

        d['analysis_type'] = experiment.sample.normal_match

        try:
            sample_type = self.get_data_from_lims(experiment.cnf, experiment.project_name, experiment.sample.name)
        except:
            pass
        else:
            if sample_type:
                d['sample_type'] = ', '.join(sample_type)

        d['genome_build'] = self.cnf.genome.name  # TODO: get genome build from the relevant project, not from the default config for this new run
        # approach_dict['min_depth'] = Metric.format_value(min_depth, is_html=True)
        if experiment.ave_depth is not None:
            d['ave_depth'] = Metric.format_value(experiment.ave_depth, is_html=True)
        return d

    @staticmethod
    def get_data_from_lims(cnf, project_name, sample_name):
        # Create the LIMS interface instance
        lims_instance = lims.Lims()
        lims_sample = None
        sample_type = None
        jira_url = cnf.jira_url
        if not jira_url:
            return None
        from source.jira_utils import retrieve_jira_info
        jira_case = retrieve_jira_info(jira_url)
        if not jira_case:
            return None
        container_id = re.findall(r'_[A-Z0-9]+XX', jira_case.data_hub)  # '/ngs/oncology/datasets/hiseq4000/160115_K00172_0044_H3GCYBBXX'
        if container_id:
            container_id = container_id[0][1:]
            lims_artifacts = lims_instance.get_artifacts(containername=container_id)
            for artifact in lims_artifacts:
                for sample in artifact.samples:
                    if sample.name.replace(' ', '_') == sample_name.replace(' ', '_'):
                        lims_sample = sample
        if lims_sample:
            sample_tissue, collection_type, tumor_type = None, None, None
            for key, value in lims_sample.udf.items():
                if key == 'Sample Tissue':
                    sample_tissue = value
                if key == 'Collection Type' and value == 'FFPE':
                    collection_type = value
                if key == 'Normal Primary Metastasis':
                    tumor_type = value
            sample_type = [sample_tissue, collection_type, tumor_type]
            sample_type = [s for s in sample_type if s]
        else:
            info('Project was not found in LIMS')
        return sample_type

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
    def _pos_recargs(chrom=None, chrom_key=None, start=None, end=None, jbrowser_link=None):
        c = (chrom.replace('chr', '')) if chrom else ''
        p = Metric.format_value(str(start) + (('...' + str(end)) if end else ''), human_readable=True, is_html=True) if start else ''
        if jbrowser_link:
            p = ('<a href="' + jbrowser_link + '&loc=chr' + c + ':' + str(start) + '...' + str(end or start) +
                 '" target="_blank">' + p + '</a>')
        return dict(value=gray(c + ':') + p, num=chrom_key * 100000000000 + start)

    cdna_chg_regexp = re.compile(r'([c,n]\.)([-\d_+*]+)(.*)')
    @staticmethod
    def _cdna_chg_recargs(mut):
        chg = mut.cdna_change
        if chg:
            p1, num, p3 = BaseClinicalReporting.cdna_chg_regexp.match(chg).groups()
            chg = (gray(str(p1)) if p1 is not None else '') + \
                  (gray(str(num)) if num is not None else '') + \
                  (str(p3) if p3 is not None else '')

        return dict(value=chg)

    @staticmethod
    def _g_chg_recargs(mut):
        chg = mut.ref + '>' + mut.alt
        if mut.var_type:
            t = mut.var_type
            if t in ['Insertion', 'Deletion']:
                t = t[:3]
            chg = gray(t) + ' ' + chg
        return dict(value=chg)

    @staticmethod
    def _hotspot_recargs(mut):
        p = Metric.format_value(mut.pos, human_readable=True, is_html=True) if mut.pos else ''
        chg = mut.ref + '>' + mut.alt
        return dict(value=gray(p + ':') + chg)

    @staticmethod
    def _signif_field(mut):
        signif = mut.signif
        if mut.reason:
            signif += gray(' (' + mut.reason + ')')
        return dict(value=signif)

    @staticmethod
    def _clinvar_recargs(mut):
        if mut.solvebio:
            return dict(value=mut.solvebio.clinsig, url=mut.solvebio.url)
        else:
            return dict(value='')

    @staticmethod
    def _highlighting_and_hiding_mut_row(row, mut):
        if not mut.signif or mut.signif.lower() in ['unknown']:
            if mut.solvebio and 'Pathogenic' in mut.solvebio.clinsig:
                warn('Mutation ' + str(mut) + ' is unknown, but found in SolveBio')
            row.hidden = True
        else:
            if mut.signif and mut.signif.lower() == 'known':
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

        bed_fname = basename(clinical_experiment_info.target.bed_fpath).split('.')[0] + '.bed' if clinical_experiment_info.target.bed_fpath else None
        jbrowser_link = get_jbrowser_link(self.cnf.genome.name, self.cnf.sample, bed_fname)

        info('Preparing data...')
        if self.experiment.mutations:
            self.mutations_report = self.make_mutations_report({self.experiment: self.experiment.mutations}, jbrowser_link)
            self.mutations_plot_data = self.make_mutations_json({self.experiment: self.experiment.mutations})
            self.substitutions_plot_data = self.make_substitutions_json({self.experiment: self.experiment.mutations})
        if self.experiment.sv_events:
            self.sv_report = self.make_sv_report({self.experiment: self.experiment.sv_events}, jbrowser_link)
        if self.experiment.seq2c_events_by_gene:
            self.seq2c_plot_data = self.make_seq2c_plot_json({self.experiment.key: self.experiment})
            if self.seq2c_plot_data:
                self.seq2c_report = self.make_seq2c_report({self.experiment: self.experiment.seq2c_events_by_gene})
        if self.experiment.actionable_genes_dict and \
                (self.experiment.mutations or self.experiment.seq2c_events_by_gene or self.experiment.sv_events):
            self.actionable_genes_report = self.make_actionable_genes_report(self.experiment.actionable_genes_dict)
        if self.experiment.ave_depth and self.experiment.depth_cutoff and self.experiment.sample.targetcov_detailed_tsv:
            self.key_genes_report = self.make_key_genes_cov_report(self.experiment.key_gene_by_name_chrom, self.experiment.ave_depth)
            self.cov_plot_data = self.make_key_genes_cov_json({self.experiment.key: self.experiment})

    def write_report(self, output_fpath):
        info('')

        data = {
            'key_or_target': self.experiment.genes_collection_type,
            'genes_description': self.experiment.genes_description,
            'sample': self.sample_section(self.experiment),
            'variants': self.__mutations_section(),
            'coverage': self.__coverage_section(),
            'actionable_genes': self.__actionable_genes_section(),
            'total_key_genes': Metric.format_value(len(self.experiment.key_gene_by_name_chrom), is_html=True)
        }
        if self.sv_report:
            data['sv'] = {}
            section = self.__sv_section()
            if section:
                data['sv'] = {'report': section, 'sv_link': relpath(self.experiment.sv_fpath, start=dirname(output_fpath))}
        if self.seq2c_plot_data:
            data['seq2c'] = {'plot_data': self.seq2c_plot_data}
            if self.seq2c_report:
                data['seq2c']['amp_del'] = self.__seq2c_section()
                if len(self.experiment.seq2c_events_by_gene.values()) > len(self.experiment.key_gene_by_name_chrom.values()):
                    data['seq2c']['description_for_whole_genomic_profile'] = {'key_or_target': self.experiment.genes_collection_type}
                    data['seq2c']['amp_del']['seq2c_switch'] = {'key_or_target': self.experiment.genes_collection_type}

        circos_plot_fpath = make_circos_plot(self.cnf, output_fpath)
        image_by_key = None
        if circos_plot_fpath:
            data['circos'] = {'circos_img': basename(circos_plot_fpath)}
            image_by_key = {'circos': circos_plot_fpath}

        write_static_html_report(self.cnf, data, output_fpath,
           tmpl_fpath=join(dirname(abspath(__file__)), 'template.html'),
           extra_js_fpaths=[join(dirname(abspath(__file__)), 'static', 'clinical_report.js'),
                            join(dirname(abspath(__file__)), 'static', 'easyzoom.js'),
                            join(dirname(abspath(__file__)), 'static', 'draw_mutations_plot.js'),
                            join(dirname(abspath(__file__)), 'static', 'draw_substitutions_plot.js'),
                            join(dirname(abspath(__file__)), 'static', 'draw_genes_coverage_plot.js'),
                            join(dirname(abspath(__file__)), 'static', 'draw_seq2c_plot.js')],
           extra_css_fpaths=[join(dirname(abspath(__file__)), 'static', 'clinical_report.css'),
                             join(dirname(abspath(__file__)), 'static', 'header_picture.css')],
           image_by_key=image_by_key)

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
            mutations_dict['plot_data'] = self.mutations_plot_data
            mutations_dict['substitutions_plot_data'] = self.substitutions_plot_data
        return mutations_dict

    def __sv_section(self):
        sv_dict = dict()
        if self.sv_report and self.sv_report.rows:
            sv_dict['table'] = build_report_html(self.sv_report, sortable=True)
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
            not_hidden_rows = [r for r in self.seq2c_report.rows if not r.hidden]
            if not_hidden_rows:
                seq2c_dict['short_table'] = self.__seq2c_create_tables(not_hidden_rows)
            seq2c_dict['full_table'] = self.__seq2c_create_tables(self.seq2c_report.rows)
        return seq2c_dict

    def __seq2c_create_tables(self, rows):
        table_dict = dict(columns=[])
        GENE_COL_NUM = min(3, len(rows))
        genes_in_col = [len(rows) / GENE_COL_NUM] * GENE_COL_NUM
        for i in range(len(rows) % GENE_COL_NUM):
            genes_in_col[i] += 1
        calc_cell_contents(self.seq2c_report, rows)
        printed_genes = 0
        for i in range(GENE_COL_NUM):
            column_dict = dict()
            column_dict['metric_names'] = [make_cell_th(m) for m in self.seq2c_report.metric_storage.get_metrics()]
            column_dict['rows'] = [
                dict(records=[make_cell_td(r) for r in region.records])
                    for region in rows[printed_genes:printed_genes + genes_in_col[i]]]
            table_dict['columns'].append(column_dict)
            printed_genes += genes_in_col[i]
        return table_dict

    def __actionable_genes_section(self):
        actionable_genes_dict = dict()
        if self.actionable_genes_report and self.actionable_genes_report.rows:
            actionable_genes_dict['table'] = build_report_html(self.actionable_genes_report, sortable=False)
        return actionable_genes_dict

    def make_key_genes_cov_report(self, key_gene_by_name_chrom, ave_depth):
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

        for gene in sorted(key_gene_by_name_chrom.values(), key=lambda g: g.name):
            reg = key_genes_report.add_row()
            reg.add_record('Gene', gene.name)
            reg.add_record('Chr', gene.chrom.replace('chr', '') if gene.chrom else None)
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
                Metric('Variant', min_width=80, max_width=120, style='white-space: pre !important;', class_='long_line'),            # p.Glu82Lys
                Metric('Type', min_width=120, max_width=120, style='white-space: pre; !important', class_='long_line'),               # Frameshift
                Metric('Types of recurrent alterations', short_name='Types of recurrent\nalterations',
                    min_width=130, max_width=130, style='white-space: pre;'),  # Mutation
                Metric('Rationale', style='max-width: 300px !important; white-space: normal;'),          # Translocations predict sensitivity
                Metric('Therapeutic Agents', max_width=120, style='white-space: normal;'),  # Sorafenib
                Metric('Freq', short_name='Freq', max_width=55, class_='shifted_column', style='white-space: pre;', with_heatmap=False)
            ])])

        report = PerRegionSampleReport(sample=self.sample, metric_storage=clinical_action_metric_storage)
        actionable_gene_names = actionable_genes_dict.keys()

        sv_mutation_types = {'Rearrangement', 'Fusion', 'Amplification', 'Deletion'}
        cnv_mutation_types = {'Amplification', 'Deletion'}

        for gene in self.experiment.key_gene_by_name_chrom.values():
            if gene.name not in actionable_gene_names:
                continue
            possible_mutation_types = set(actionable_genes_dict[gene.name][1].split('; '))
            # possible_mutation_types = possible_mutation_types - sv_mutation_types
            # if not possible_mutation_types: continue

            variants = []
            types = []
            frequencies = []

            if self.experiment.mutations:
                vardict_mut_types = possible_mutation_types - sv_mutation_types - cnv_mutation_types
                if vardict_mut_types:
                    for mut in self.experiment.mutations:
                        if mut.gene.name == gene.name:
                            if mut.signif not in ['unknown'] and mut.is_canonical:
                                variants.append(mut.aa_change if mut.aa_change else '.')
                                types.append(mut.var_type)
                                frequencies.append(Metric.format_value(mut.freq, unit='%'))

            if cnv_mutation_types:
                for se in gene.seq2c_events:
                    if 'Amplification' in possible_mutation_types and se.amp_del == 'Amp' or \
                            'Deletion' in possible_mutation_types and se.amp_del == 'Del':
                        variants.append(se.amp_del + ', ' + se.fragment)
                        types.append(se.amp_del)
                        frequencies.append('')

            if sv_mutation_types:
                svs_by_key = OrderedDict()
                for sv in gene.sv_events:
                    if any(a.known or a.effect == 'EXON_DEL' for a in sv.annotations):
                        svs_by_key[sv.type, tuple(tuple(sorted(a.genes)) for a in sv.annotations)] = sv

                for se in svs_by_key.values():
                    if ('Fusion' in possible_mutation_types or 'Rearrangement' in possible_mutation_types) and se.type == 'BND' or \
                       'Deletion' in possible_mutation_types and se.type == 'DEL' or \
                       'Amplification' in possible_mutation_types and se.type == 'DUP':
                        variants.append(', '.join(set('/'.join(set(a.genes)) for a in se.key_annotations if a.genes)))
                        types.append(BaseClinicalReporting.sv_type_dict.get(se.type, se.type))
                        frequencies.append('')

            if not variants:
                continue

            reg = report.add_row()
            reg.add_record('Gene', gene.name)
            reg.add_record('Variant', '\n'.join(variants))
            reg.add_record('Type', '\n'.join(types))
            reg.add_record('Types of recurrent alterations', actionable_genes_dict[gene.name][1].replace('; ', '\n'))
            reg.add_record('Rationale', actionable_genes_dict[gene.name][0])
            reg.add_record('Therapeutic Agents', actionable_genes_dict[gene.name][2])
            reg.add_record('Freq', '\n'.join(frequencies))

        return report


def make_circos_plot(cnf, output_fpath):
    info('Making circos plot')
    circos_py_executable = get_script_cmdline(cnf, 'python', join('scripts', 'post', 'circos.py'))
    cmdline = '{circos_py_executable} '

    required_files = [cnf.mutations_fpath, cnf.seq2c_tsv_fpath, cnf.sv_vcf_fpath]
    for file, desc in zip(required_files, ['Vardict results', 'Seq2C results', 'SV calling results']):
        if not file or not verify_file(file):
            warn('File with ' + desc + ' is not found. Circos plot cannot be created.')
            return None

    vardict2mut_raw_fpath = cnf.mutations_fpath.replace(source.mut_pass_suffix + '.', '')
    output_dir = dirname(output_fpath)

    if cnf.bed_fpath:
        cmdline += ' --bed ' + cnf.bed_fpath
    cmdline += ' --mutations ' + vardict2mut_raw_fpath
    cmdline += ' --seq2c ' + cnf.seq2c_tsv_fpath
    cmdline += ' --sv ' + cnf.sv_vcf_fpath

    cmdline += ' --sample ' + cnf.sample
    cmdline += ' --genome ' + cnf.genome.name
    cmdline += ' -o ' + output_dir

    cmdline = cmdline.format(**locals())
    res = call(cnf, cmdline, stdout_to_outputfile=False, exit_on_error=False)
    circos_plot_fpath = join(output_dir, cnf.sample + '.png')
    if not exists(circos_plot_fpath):
        return None

    return circos_plot_fpath

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
