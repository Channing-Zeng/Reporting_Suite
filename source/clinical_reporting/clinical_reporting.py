from collections import OrderedDict, defaultdict
import json
from genericpath import exists
from os.path import join, dirname, abspath, relpath, basename

import re

import source
from source import info, verify_file
from source.calling_process import call
from source.clinical_reporting.clinical_parser import get_group_num
from source.clinical_reporting.utils import SVEvent
from source.clinical_reporting.combine_clinical_reporting_utils import get_vcf_readers, add_freq_depth_records,\
    group_for_venn_diagram, update_venn_diagram_data, save_venn_diagram_data, format_experiment_names, add_tooltip
from source.logger import warn, err, debug
from source.reporting.reporting import MetricStorage, Metric, PerRegionSampleReport, ReportSection, calc_cell_contents, make_cell_td, write_static_html_report, make_cell_th, build_report_html
from source.tools_from_cnf import get_script_cmdline
from source.utils import get_chr_lengths, OrderedDefaultDict, is_us, is_uk, is_local, gray, get_version
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


class ActionableVariant:
    def __init__(self):
        self.gene = None
        self.var_type = None
        self.status = ''
        self.freq = ''


class BaseClinicalReporting:
    def __init__(self, cnf, *args):
        self.cnf = cnf
        self.chromosomes_by_name = Chromosome.build_chr_by_name(self.cnf)

    def make_mutations_report(self, mutations_by_experiment, jbrowser_link, samples_data=None, parameters_info=None,
                              create_venn_diagrams=False, cur_group_num=None):
        full_names = []
        if len(mutations_by_experiment) == 1:
            ms = [
                Metric('Gene'),  # Gene & Transcript
                Metric('AA len', max_width=50, class_='stick_to_left', with_heatmap=False),          # 128
                Metric('AA chg', short_name='AA change', max_width=70, class_='long_line'),            # p.Glu82Lys
                Metric('Position', with_heatmap=False, align='left', sort_direction='ascending'),       # g.47364249
                Metric('Change', max_width=100, class_='long_line', description='Genomic change'),       # G>A
                Metric('cDNA change', class_='long_line', description='cDNA change'),       # G>A
                Metric('MSI', short_name='HP', description='Microsatellite instability length', quality='Less is better', with_heatmap=False),
                Metric('Status', short_name='Status'),     # Somatic
                Metric('Effect', max_width=100, class_='long_line'),               # Frameshift
                Metric('VarDict status', short_name='Significance', max_width=230, class_='long_line'),     # Likely
                Metric('Databases'),                 # rs352343, COSM2123, SolveBio
                Metric('Samples', with_heatmap=False),          # 128
                Metric('Other occurrences', class_='long_line', with_heatmap=False),          # 128
                # Metric('ClinVar', short_name='SolveBio ClinVar'),
                Metric('Freq', short_name='Freq', max_width=55, unit='%', with_heatmap=False),          # .19
                Metric('Depth', short_name='Depth', max_width=48, med=mutations_by_experiment.keys()[0].ave_depth, with_heatmap=False),              # 658
                Metric('Indicentalome', short_name='Callability issues')
            ]

        else:
            ms = [
                Metric('Gene'),  # Gene & Transcript
                Metric('AA chg', short_name='AA change', max_width=70, class_='long_line'),            # p.Glu82Lys
                Metric('Effect', max_width=100, class_='long_line'),               # Frameshift
            ]
            short_names, full_names = format_experiment_names(mutations_by_experiment, samples_data, cur_group_num)
            used_full_names = set()
            for index in range(len(short_names.values())):
                short_name = short_names.values()[index]
                full_name = full_names.values()[index]
                next_short_name = short_names.values()[index + 1] if index < len(short_names.values()) - 1 else ''
                if full_name in used_full_names:
                    continue
                used_full_names.add(full_name)
                col_width = 20 + 15 * len(next_short_name.split())

                _get_class = lambda _name: ' '.join(['td_' + n for n in _name.lower().split(' ')])
                freq_name = full_name + ' Freq'
                depth_name = full_name + ' Depth'
                ms.extend([
                    Metric(freq_name, short_name='\n'.join(short_name.split()) + '\nfreq', max_width=45,  min_width=45,
                           align='left', unit='%', with_heatmap=False,
                           td_class=_get_class(freq_name)),          # .19
                    Metric(depth_name, short_name='depth', min_width=col_width, align='left',
                           med=mutations_by_experiment.keys()[0].ave_depth, with_heatmap=False,
                           td_class=_get_class(depth_name)),              # 658
                ])
            ms.extend([
                Metric('Samples', with_heatmap=False),          # 128
                Metric('VarDict status', short_name='Significance', max_width=230, class_='long_line'),     # Likely
                Metric('Other occurrences', class_='long_line', with_heatmap=False),          # 128
                # Metric('ClinVar', short_name='SolveBio ClinVar'),
                Metric('Indicentalome', short_name='Callability issues', class_='long_line'),
                Metric('Databases'),                 # rs352343, COSM2123, SolveBio
                Metric('Status', short_name='Status'),     # Somatic
                Metric('Position', with_heatmap=False, align='left', sort_direction='ascending'),       # g.47364249
                Metric('AA len', max_width=50, class_='stick_to_left', with_heatmap=False),          # 128
                Metric('Change', max_width=100, class_='long_line', description='Genomic change'),       # G>A
                Metric('cDNA change', class_='long_line', description='cDNA change'),       # G>A
                Metric('MSI', short_name='HP', description='Microsatellite instability length', quality='Less is better', with_heatmap=False),
            ])
            if parameters_info:
                for parameter_name in parameters_info.keys():
                    ms.append(Metric(parameter_name, is_hidden=True))

        venn_sets = OrderedDefaultDict(int)

        if create_venn_diagrams:
            samples_by_index, set_labels = group_for_venn_diagram(mutations_by_experiment, full_names, parameters_info, samples_data)

        clinical_mut_metric_storage = MetricStorage(sections=[ReportSection(metrics=ms, name='mutations')])
        report = PerRegionSampleReport(sample=mutations_by_experiment.keys()[0].sample,
            metric_storage=clinical_mut_metric_storage, expandable=True)

        # Writing records
        print_cdna = False
        muts_by_key_by_experiment = OrderedDefaultDict(OrderedDict)
        sample_experiments = []
        if samples_data:
            vcf_readers, filt_vcf_readers = get_vcf_readers(mutations_by_experiment, cur_group_num)
        for e, muts in mutations_by_experiment.items():
            if samples_data and (not cur_group_num or get_group_num(e.key) == cur_group_num):
                sample_experiments.append(e)
            for mut in muts:
                muts_by_key_by_experiment[mut.get_key()][e] = mut
                if mut.cdna_change.strip():
                    print_cdna = True

        # canonical_mutations = [
        #     [m.is_canonical if m is not None else False for m in muts]
        #      for e, muts in mutations_by_experiment.items()]
        #
        # mut_positions = [m.pos for i, (e, muts) in enumerate(mutations_by_experiment.items())
        #                  for j, m in enumerate(muts) if m is not None and canonical_mutations[i][j]]

        for mut_key, mut_by_experiment in muts_by_key_by_experiment.items():
            mut = next((m for m in mut_by_experiment.values() if m is not None), None)
            row_class = ''

            # if mut.pos not in mut_positions:
            #     mut_positions.append(mut.pos)
            #     row = report.add_row()
            #     row.add_record('Gene', mut.gene.name)
            #     row_class = ' expandable_gene_row collapsed'
            #     row.add_record('Position',
            #         **self._pos_recargs(mut.chrom, mut.get_chrom_key(), mut.pos, mut.pos, jbrowser_link))
            #     row.add_record('Change', **self._g_chg_recargs(mut))
            #
            #     if len(mutations_by_experiment.values()) == 1:
            #         row.add_record('Freq', mut.freq if mut else None)
            #         row.add_record('Depth', mut.depth if mut else None)
            #     else:
            #         for e, m in mut_by_experiment.items():
            #             row.add_record(e.key + ' Freq', m.freq if m else None)
            #             row.add_record(e.key + ' Depth', m.depth if m else None)
            #     row.class_ = row_class
            #     self._highlighting_and_hiding_mut_row(row, mut, freq_in_samples)
            #     if len(mutations_by_experiment) > 1 and len(mut_by_experiment.keys()) == len(mutations_by_experiment.keys()):
            #         row.highlighted_green = True

            if not mut.is_canonical:
                continue
            if mut.is_silent:
                continue
            if len(mutations_by_experiment.values()) > 1 and \
                    (not mut.signif or mut.signif.lower() in ['unknown']):
                continue
            if sample_experiments and all(e not in sample_experiments for e in mut_by_experiment.keys()):
                continue
            elif sample_experiments:
                cur_experiments = [e for e in mut_by_experiment.keys() if e in sample_experiments]
            row = report.add_row()
            row.add_record('Gene', **self._gene_recargs(mut))
            # if mut.is_canonical:
            #     row.add_record('Transcript', mut.transcript)
            #     row_class = ' expandable_gene_row collapsed'
            # else:
            #     row.add_record('Transcript', mut.transcript)
            #     row_class = ' row_to_hide row_hidden'
            # row.add_record('AA len', )
            row.add_record('AA chg', **self._aa_chg_recargs(mut))
            row.add_record('Position', show_content=mut.is_canonical,
               **self._pos_recargs(mut.chrom, mut.get_chrom_key(), mut.pos, None, jbrowser_link))
            row.add_record('Change', show_content=mut.is_canonical, **self._g_chg_recargs(mut))
            if print_cdna:
                row.add_record('cDNA change', **self._cdna_chg_recargs(mut))
            if mut.msi > 3:
                row.add_record('MSI', mut.msi)
            if mut.status:
                row.add_record('Status', mut.status)
            row.add_record('Effect', mut.eff_type.replace(' variant', '') if mut.eff_type else None)
            row.add_record('VarDict status', **self._significance_field(mut))
            if mut.incidentalome_reason:
                row.add_record('Indicentalome', value=mut.incidentalome_reason)
            # row.add_record('VarDict reason', mut.reason)
            row.add_record('Databases', **self._db_recargs(mut))
            # row.add_record('ClinVar', **self._clinvar_recargs(mut))

            if len(mutations_by_experiment.values()) == 1:
                row.add_record('Freq', mut.freq if mut else None, show_content=mut.is_canonical)
                row.add_record('Depth', mut.depth if mut else None, show_content=mut.is_canonical)
            else:
                add_freq_depth_records(row, mut, mut_by_experiment, full_names, cur_group_num, samples_data,
                                            parameters_info, vcf_readers, filt_vcf_readers)

            self._highlighting_and_hiding_mut_row(row, mut)

            if len(mutations_by_experiment.values()) > 1:
                row.class_ = row.class_.replace('incidentalome', '')
            if len(mut_by_experiment.keys()) > 1:
                # k = float(len(mut_by_experiment.keys())) / len(mutations_by_experiment.keys())
                # row.color = 'hsl(100, 100%, ' + str(70 + int(30 * (1 - k))) + '%)'
                row.class_ += ' multiple_occurences_row'
            if create_venn_diagrams:
                update_venn_diagram_data(venn_sets, mut_by_experiment, samples_by_index, full_names)

            row.class_ += ' ' + row_class

        if create_venn_diagrams:
            venn_data = save_venn_diagram_data(venn_sets, set_labels)
            return report, venn_data

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
                    chrom=mut.chrom, position=mut.pos, freq=mut.freq * 100, status=mut.signif,
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
                if mut.signif in ['known', 'likely']:
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
            Metric('Genes', min_width=100),
            # Metric('Chr', with_heatmap=False, max_width=50, align='left', sort_direction='ascending'),
            Metric('Type', min_width=50),
            Metric('Location', with_heatmap=False, align='left', sort_direction='ascending', style="width: 550px"),
            Metric('Transcript', align='left', style="width: 300px"),
            Metric('Priority', is_hidden=True),
            Metric('Reads', is_hidden=True),
            # Metric('Status', with_heatmap=False, align='left', sort_direction='ascending'),
            # Metric('Effects', with_heatmap=False, align='left', class_='long_line'),
        ]

        if len(svs_by_experiment) == 1:
            ms.extend([
                Metric('Split read support', short_name='Split /', min_width=30),
                Metric('Paired read support', short_name='paired read support', min_width=50),
            ])
        else:
            for e in svs_by_experiment.keys():
                ms.extend([
                    Metric(e.key + ' Split read support', short_name='Split /', min_width=30),
                    Metric(e.key + ' Paired read support', short_name='paired read support', min_width=50),
                ])

        metric_storage = MetricStorage(sections=[ReportSection(name='main_sv_section', metrics=ms)])
        report = PerRegionSampleReport(sample=svs_by_experiment.keys()[0].sample,
                                       metric_storage=metric_storage, expandable=True, large_table=True)

        # Writing records
        svanns_by_key_by_experiment = OrderedDefaultDict(lambda : OrderedDefaultDict(SVEvent.Annotation))
        for e, svs in svs_by_experiment.items():
            known_cnt = 0
            exon_dels_cnt = 0
            fusions_cnt = 0
            other_cnt = 0
            for sv_event in svs:
                for an in sv_event.key_annotations:
                    # reporting all known (fusion) by default
                    if an.known:
                        if e not in svanns_by_key_by_experiment[an.get_key()]:
                            known_cnt += 1
                        svanns_by_key_by_experiment[an.get_key()][e].update_annotation(an)

                    # reporting all whole exon deletions
                    elif sv_event.is_deletion() and ('exon_del' in an.effect.lower() or 'exon_loss' in an.effect.lower()) \
                            and (not an.priority or an.priority == SVEvent.Annotation.ON_PRIORITY_LIST):
                        if e not in svanns_by_key_by_experiment[an.get_key()]:
                            exon_dels_cnt += 1
                        svanns_by_key_by_experiment[an.get_key()][e].update_annotation(an)

                    # reporting fusions in the AZ priority genes
                    elif sv_event.is_fusion() and (not an.priority or an.priority == SVEvent.Annotation.ON_PRIORITY_LIST):
                        if e not in svanns_by_key_by_experiment[an.get_key()]:
                            fusions_cnt += 1
                        svanns_by_key_by_experiment[an.get_key()][e].update_annotation(an)

                    # # reporting all non-fusion events affecting 2 or more genes (start and end should not be the same gene. handling overlapping gene cases.)
                    # elif sv_event.end_genes and all(ann_g not in sv_event.end_genes for ann_g in an.genes):
                    #     svanns_by_key_by_experiment[an.get_key()][e] = an
                    #     affecting_2_genes_cnt += 1

                    else:
                        other_cnt += 1

            info('known_cnt: ' + str(known_cnt))
            info('exon_dels_cnt: ' + str(exon_dels_cnt))
            info('fusions_cnt: ' + str(fusions_cnt))
            # info('affecting_2_genes_cnt: ' + str(affecting_2_genes_cnt))
            info('other_cnt: ' + str(other_cnt))

        for sv_key, svann_by_experiment in svanns_by_key_by_experiment.items():
            sv_ann = next((s for s in svann_by_experiment.values() if s is not None), None)

            if not sv_ann or not sv_ann.event:
                continue

            row = report.add_row()

            row.add_record('Genes', '/'.join(set(sv_ann.genes)))

            type_str = BaseClinicalReporting.sv_type_dict.get(sv_ann.event.type, sv_ann.event.type)
            if sv_ann.effect == 'EXON_DEL':
                type_str += ' ' + sv_ann.exon_info
            if sv_ann.priority == SVEvent.Annotation.KNOWN:
                type_str = 'Known fusion'
            else:
                row.class_ += ' depth_filterable'
                if sv_ann.event.read_support < SVEvent.min_sv_depth:
                    row.class_ += ' less_threshold'
            if type_str == 'Fusion':
                row.hidden = True

            row.add_record('Type', type_str)
            row.add_record('Priority', sv_ann.priority)

            if sv_ann.event.start:
                row.add_record('Location',
                    **self._pos_recargs(sv_ann.event.chrom, sv_ann.event.get_chrom_key(),
                                        sv_ann.event.start, sv_ann.event.end, jbrowser_link, end_chrom=sv_ann.event.chrom2))
            else:
                row.add_record('Location', None)

            row.add_record('Transcript', sv_ann.transcript)
            # row.add_record('Status', 'known' if any(a.known for a in sv.annotations) else None)
            # row.add_record('Effects', ', '.join(set(a.effect.lower().replace('_', ' ') for a  in sv.annotations if a in ['EXON_DEL', 'FUSION'])))

            row.add_record('Reads', sv_ann.event.read_support)
            if len(svann_by_experiment.values()) == 1:
                row.add_record('Split read support', **self._reads_recargs(sv_ann.event.split_read_support))
                row.add_record('Paired read support', **self._reads_recargs(sv_ann.event.paired_end_support))
            else:
                for _sv_ann, m in svann_by_experiment.items():
                    if _sv_ann:
                        row.add_record(e.key + ' Split read support', **self._reads_recargs(_sv_ann.event.split_read_support))
                        row.add_record(e.key + ' Paired read support', **self._reads_recargs(_sv_ann.event.paired_end_support))

            if any(an.known for an in svann_by_experiment.values()):
                row.highlighted = True

            if len(svann_by_experiment.keys()) == len(svs_by_experiment.keys()):
                row.highlighted_green = True

        return report

    def make_seq2c_plot_json(self, experiment_by_key):
        data = dict()

        for k, e in experiment_by_key.items():
            if not e.seq2c_events_by_gene:
                continue
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

            # if all(e['ampDel'] is None for e in d['events']):
            #     d = None
            data[k.lower()] = d

        if len(experiment_by_key.keys()) == 1:
            return json.dumps(data.values()[0]) if data.values()[0] else None
        else:
            return json.dumps(data) if all(d is not None for d in data.values()) else None

    def make_seq2c_report(self, seq2c_by_experiments, samples_data=None, cur_group_num=None):
        ms = [
            Metric('Gene', align='left', sort_direction='ascending'),  # Gene
            Metric('Chr', with_heatmap=False, max_width=20, align='right'),
            Metric('Log ratio', med=0, quality='Less is better'),  # for consistency with plot: red for amplifications, blue for deletions
            Metric('Amp/Del'),
            Metric('BP/Whole'),
        ]
        if len(seq2c_by_experiments.values()) > 1:
            if cur_group_num:
                ms = [
                    Metric('Gene', align='left', sort_direction='ascending'),  # Gene
                    Metric('Chr', with_heatmap=False, max_width=20, align='right'),
                    Metric('Amp/Del'),
                    Metric('BP/Whole'),
                ]
                short_names, full_names = format_experiment_names(seq2c_by_experiments, samples_data, cur_group_num)
                for index in range(len(short_names.values())):
                    short_name = short_names.values()[index]
                    full_name = full_names.values()[index]
                    ms.append(Metric(full_name + ' log ratio', short_name='\n'.join(short_name.split()) + '\nlog ratio',
                                            med=0, quality='Less is better', align='left'))
            else:
                ms.append(Metric('Samples', max_width=25, align='left'))
        metric_storage = MetricStorage(sections=[ReportSection(name='seq2c_section', metrics=ms)])

        report = PerRegionSampleReport(sample=seq2c_by_experiments.keys()[0].sample,
            metric_storage=metric_storage, expandable=True)

        # Writing records
        seq2c_by_key_by_experiment = OrderedDefaultDict(OrderedDict)
        key_gene_by_name_chrom = []
        for e, seq2c in seq2c_by_experiments.items():
            key_gene_by_name_chrom = e.key_gene_by_name_chrom
            is_whole_genomic_profile = len(e.seq2c_events_by_gene.values()) > len(key_gene_by_name_chrom.values())
            events = seq2c.values()
            for event in events:
                if event.is_amp() or event.is_del():
                    seq2c_by_key_by_experiment[(event.gene.name, event.amp_del)][e] = event

        seq2c_by_key_by_experiment = OrderedDict(sorted(seq2c_by_key_by_experiment.iteritems(), key=lambda x: x[0][0]))
        for seq2c_by_experiment_key, seq2c_by_experiment in seq2c_by_key_by_experiment.items():
            event = next((e for e in seq2c_by_experiment.values() if e is not None), None)
            if len(seq2c_by_experiments.values()) > 1:
                is_event_in_this_report = False
                for e, event in seq2c_by_experiment.iteritems():
                    if e in full_names:
                        is_event_in_this_report = True
                        break
                if not is_event_in_this_report:
                    continue
            row = report.add_row()
            row.add_record('Gene', event.gene.name)
            row.add_record('Chr', event.gene.chrom.replace('chr', ''))
            if len(seq2c_by_experiments.values()) == 1:
                row.add_record('Log ratio', '%.2f' % (event.ab_log2r if event.ab_log2r is not None else event.log2r))
            else:
                for e, event in seq2c_by_experiment.iteritems():
                    if e not in full_names:
                        continue
                    full_name = full_names[e]
                    row.add_record(full_name + ' log ratio',
                                   '%.2f' % (event.ab_log2r if event.ab_log2r is not None else event.log2r))

            row.add_record('Amp/Del', event.amp_del)
            row.add_record('BP/Whole', event.fragment)
            if is_whole_genomic_profile and event.gene.key not in key_gene_by_name_chrom:
                row.hidden = True

            if len(seq2c_by_experiments.values()) > 1 and not cur_group_num:
                num_samples = len(seq2c_by_experiment.keys())
                row.add_record('Samples', num_samples)

        return report

    def make_key_genes_cov_json(self, experiment_by_key):
        chr_cum_lens = Chromosome.get_cum_lengths(self.chromosomes_by_name)
        chr_cum_len_by_chrom = dict(zip([c.name for c in self.chromosomes_by_name.values()], chr_cum_lens))

        gene_names = []
        transcript_names = []
        strands = []

        ticks_x = [[(chr_cum_lens[i] + chr_cum_lens[i + 1])/2, self.chromosomes_by_name.values()[i].short_name]
                   for i in range(len(self.chromosomes_by_name.keys()))]

        hits = list()
        for key, e in experiment_by_key.items():
            gene_ave_depths = []
            covs_in_thresh = []
            coords_x = []
            cds_cov_by_gene = defaultdict(list)
            mut_info_by_gene = dict()

            for gene in e.key_gene_by_name_chrom.values():
                mut_info_by_gene[gene.name] = [('p.' + m.aa_change if m.aa_change else '.') for m in gene.mutations if m.is_canonical]
                for se in gene.seq2c_events:
                    if se and (se.is_amp() or se.is_del()):
                        mut_info_by_gene[gene.name].append(se.amp_del + ', ' + se.fragment)

            for gene in e.key_gene_by_name_chrom.values():
                if not gene.cov_by_threshs:
                    # err('Gene ' + gene.name + ' has no cov_by_threshs')
                    continue
                gene_names.append(gene.name)
                transcript_names.append(gene.transcript_id)
                strands.append(gene.strand)
                gene_ave_depths.append(gene.ave_depth)
                c = gene.cov_by_threshs.get(e.region_depth_cutoff)
                assert c is not None
                covs_in_thresh.append(c)
                if not gene.start or not gene.end or not gene.chrom:
                    err('Gene ' + gene.name + ': no chrom or start or end specified')
                    continue
                coords_x.append(chr_cum_len_by_chrom[gene.chrom] + gene.start + (gene.end - gene.start) / 2)
                cds_cov_by_gene[gene.name] = [dict(
                    start=cds.start,
                    end=cds.end,
                    aveDepth=cds.ave_depth,
                    percentInThreshold=cds.cov_by_threshs.get(e.depth_cutoff),
                ) for cds in gene.cdss]

            hits.append(dict(
                key=key.lower(),
                coords_x=coords_x,
                gene_ave_depths=gene_ave_depths,
                covs_in_thresh=covs_in_thresh,
                cds_cov_by_gene=cds_cov_by_gene,
                mut_info_by_gene=mut_info_by_gene
            ))

        data = dict(
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

    def sample_section(self, experiment, use_abs_report_fpath=False, sample_name=None):
        d = dict()
        d['patient'] = {'sex': 'unknown'}
        d['project_report_rel_path'] = 'not generated'
        d['project_dirpath'] = 'not generated'
        # d['panel'] = 'unknown'
        # d['bed_path'] = 'unknown'
        # d['target_type'] = 'unknown'
        # d['target_fraction'] = 'unknown'
        # d['ave_depth'] = 'unknown'

        d['key'] = experiment.key
        d['sample'] = sample_name or experiment.sample.name.replace('_', ' ')
        if experiment.patient and experiment.patient.gender:
            d['patient'] = {'sex': experiment.patient.gender}
        d['project_name'] = experiment.project_name.replace('_', ' ')
        if experiment.project_report_path:
            if use_abs_report_fpath:
                d['project_report_rel_path'] = experiment.project_report_path
            else:
                d['project_report_rel_path'] = relpath(experiment.project_report_path, dirname(experiment.sample.clinical_html))
            d['project_dirpath'] = experiment.cnf.output_dir
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

    def seq2c_section(self):
        seq2c_dict = dict()
        seq2c_report = self.seq2c_report
        rows = seq2c_report.get_rows_of_records()
        if seq2c_report and rows:
            not_hidden_rows = [r for r in rows if not r.hidden]
            if not_hidden_rows:
                seq2c_dict['short_table'] = self.seq2c_create_tables(not_hidden_rows)
            seq2c_dict['full_table'] = self.seq2c_create_tables(seq2c_report.rows)
        return seq2c_dict

    def seq2c_create_tables(self, rows):
        table_dict = dict(columns=[])
        GENE_COL_NUM = min(3, len(rows))
        genes_in_col = [len(rows) / GENE_COL_NUM] * GENE_COL_NUM
        for i in range(len(rows) % GENE_COL_NUM):
            genes_in_col[i] += 1
        calc_cell_contents(self.seq2c_report, rows)
        printed_genes = 0
        for i in range(GENE_COL_NUM):
            column_dict = dict()
            metrics = self.seq2c_report.metric_storage.get_metrics()
            column_dict['metric_names'] = [make_cell_th(m) for m in metrics]
            column_dict['rows'] = [
                dict(records=[make_cell_td(r) for r in region.records])
                    for region in rows[printed_genes:printed_genes + genes_in_col[i]]]
            table_dict['columns'].append(column_dict)
            printed_genes += genes_in_col[i]
        return table_dict

    @staticmethod
    def get_data_from_lims(cnf, project_name, sample_name):
        # Create the LIMS interface instance
        try:
            from genologics import lims
        except:
            return None
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
        def filter_digits(s):
            return ''.join(c for c in s if c.isdigit())

        val = OrderedDict()

        if mut.cosmic_id:
            val['COSM'] = 'http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=' + \
                                 filter_digits(mut.cosmic_id)
        if mut.dbsnp_id:
            val['dbSNP'] = 'http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs=' + \
                                filter_digits(mut.dbsnp_id)

        return dict(value=None, html_fpath=val)

    @staticmethod
    def _aa_chg_recargs(mut):
        aa_chg = ''.join([gray(c) if c.isdigit() else c for c in (mut.aa_change or '')])
        return dict(value=aa_chg)

    @staticmethod
    def _gene_recargs(mut):
        t = ''
        if is_us() or is_uk() or is_local():  # add button to comment mutation
            t += '<div style="position:relative;"><div class="comment_div" onclick="commentMutation($(this))"></div></div> '

        if mut.transcript:
            tooltip = mut.transcript
            if mut.aa_len:
                tooltip += '<br>AA length: ' + str(mut.aa_len)
            if mut.exon:
               tooltip += '<br>Exon: ' + str(mut.exon)
            t += add_tooltip(mut.gene.name, tooltip)
            str(mut.aa_len)
        else:
            t += mut.gene.name
        return dict(value=t, show_content=mut.is_canonical)

    @staticmethod
    def _pos_recargs(chrom=None, chrom_key=None, start=None, end=None, jbrowser_link=None, end_chrom=None):
        c = (chrom.replace('chr', '')) if chrom else ''
        if not end or not end_chrom:
            p_html = Metric.format_value(start, human_readable=True, is_html=True) + \
                ('-' + Metric.format_value(end, human_readable=True, is_html=True) if end else '') if start else ''
            p_html = gray(c + ':') + p_html
            if jbrowser_link:
                p_html = ('<a href="' + jbrowser_link + '&loc=chr' + c + ':' + str(start) + '...' + str(end or start) +
                     '" target="_blank">' + p_html + '</a>')
        else:
            end_c = (end_chrom.replace('chr', ''))
            start_html = gray(c + ':') + Metric.format_value(start, human_readable=True, is_html=True)
            end_html = gray(end_c + ':') + Metric.format_value(end, human_readable=True, is_html=True)
            if jbrowser_link:
                p_html = ('<a href="' + jbrowser_link + '&loc=chr' + c + ':' + str(start) + '" target="_blank">' + start_html + '</a>')
                p_html += '-'
                p_html += '<a href="' + jbrowser_link + '&loc=chr' + end_c + ':' + str(end) + '" target="_blank">' + end_html + '</a>'
        return dict(value=p_html, num=chrom_key * 100000000000 + start)

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

        if mut.solvebio_url:
            chg = '<a href="' + mut.solvebio_url + '" target="_blank">' + chg + '</a>'

        return dict(value=chg)

    @staticmethod
    def _hotspot_recargs(mut):
        p = Metric.format_value(mut.pos, human_readable=True, is_html=True) if mut.pos else ''
        chg = mut.ref + '>' + mut.alt
        return dict(value=gray(p + ':') + chg)

    @staticmethod
    def _significance_field(mut):
        txt = mut.signif
        if mut.reason:
            txt += '<span class="reason"> (' + mut.reason + ') </span>'
        txt = '<span class="span_status_' + mut.signif + '">' + txt + '</span>'
        return dict(value=txt)

    @staticmethod
    def _highlighting_and_hiding_mut_row(row, mut):
        if mut.incidentalome_reason:
            row.class_ += ' incidentalome'
        if not mut.signif or mut.signif.lower() in ['unknown']:
            # if mut.solvebio and 'Pathogenic' in mut.solvebio.clinsig:
            #     warn('Mutation ' + str(mut) + ' is unknown, but found in SolveBio as Pathogenic')
            row.hidden = True
        else:
            if mut.signif:
                row.class_ += ' ' + mut.signif.lower()
        return row

    @staticmethod
    def _reads_recargs(read_support_dict):
        tooltip = ''
        for caller, reads_num in read_support_dict.iteritems():
            reads = str(reads_num) if reads_num is not None else 'no reads'
            tooltip += caller + ': ' + reads + '<br>'
        read_support_num = read_support_dict['manta'] if 'manta' in read_support_dict else read_support_dict['lumpy']
        read_support = str(read_support_num) if read_support_num is not None else '.'
        t = add_tooltip(read_support, tooltip)
        return dict(value=t, num=read_support_num)


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
        jbrowser_link = get_jbrowser_link(self.cnf.genome.name, self.cnf.sample, [bed_fname])

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
        if self.experiment.ave_depth is not None and self.experiment.depth_cutoff is not None and self.experiment.sample.targetcov_detailed_tsv:
            self.key_genes_report = self.make_key_genes_cov_report(self.experiment.key_gene_by_name_chrom, self.experiment.ave_depth)
            self.cov_plot_data = self.make_key_genes_cov_json({self.experiment.key: self.experiment})

    def write_report(self, output_fpath):
        info('')

        data = {
            'software_version': get_version(),
            'key_or_target': self.experiment.genes_collection_type,
            'genes_description': self.experiment.genes_description,
            'sample': self.sample_section(self.experiment),
            'variants': self.__mutations_section(),
            'coverage': self.__coverage_section(),
            'actionable_genes': self.__actionable_genes_section(),
            'total_key_genes': Metric.format_value(self.experiment.key_genes_number, is_html=True)
        }
        if self.sv_report:
            data['sv'] = {}
            section = self.__sv_section()
            if section:
                data['sv'] = {'report': section, 'sv_link': relpath(self.experiment.sv_fpath, start=dirname(output_fpath)),
                              'default_sv_depth': SVEvent.min_sv_depth}
        if self.seq2c_plot_data:
            data['seq2c'] = {'plot_data': self.seq2c_plot_data}
            if self.seq2c_report:
                data['seq2c']['amp_del'] = self.seq2c_section()
                if len(self.experiment.seq2c_events_by_gene.values()) > len(self.experiment.key_gene_by_name_chrom.values()):
                    data['seq2c']['description_for_whole_genomic_profile'] = {'key_or_target': self.experiment.genes_collection_type}
                    data['seq2c']['amp_del']['seq2c_switch'] = {'key_or_target': self.experiment.genes_collection_type}

        min_freq = self.cnf.variant_filtering.min_freq
        act_min_freq = self.cnf.variant_filtering.act_min_freq
        data['min_af'] = str(float(min_freq) * 100)
        data['act_min_af'] = str(float(act_min_freq) * 100)

        circos_plot_fpath = make_circos_plot(self.cnf, output_fpath)
        image_by_key = None
        if circos_plot_fpath:
            data['circos'] = {'circos_img': basename(circos_plot_fpath)}
            image_by_key = {'circos': circos_plot_fpath}

        comment_php_path = 'http://ngs.usbod.astrazeneca.net/save_comment.php'
        if is_uk():
            comment_php_path = '/ngs/reports/save_comment.php'
        data['comment_php_path'] = json.dumps(comment_php_path)

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
                Metric('Freq', short_name='Freq', max_width=55, unit='%', style='white-space: pre;', with_heatmap=False),
                Metric('VarDict status', short_name='Significance', with_heatmap=False)
            ], name='actionable')])

        report = PerRegionSampleReport(sample=self.sample, metric_storage=clinical_action_metric_storage, keep_order=True)
        actionable_gene_names = actionable_genes_dict.keys()

        sv_mutation_types = {'Rearrangement', 'Fusion', 'Amplification', 'Deletion'}
        cnv_mutation_types = {'Amplification', 'Deletion'}

        for gene in sorted(self.experiment.key_gene_by_name_chrom.values(), key=lambda x: x.name):
            if gene.name not in actionable_gene_names:
                continue
            possible_mutation_types = set(actionable_genes_dict[gene.name][1].split('; '))
            # possible_mutation_types = possible_mutation_types - sv_mutation_types
            # if not possible_mutation_types: continue
            actionable_variants = []

            if self.experiment.mutations:
                vardict_mut_types = possible_mutation_types - sv_mutation_types - cnv_mutation_types
                if vardict_mut_types:
                    for mut in self.experiment.mutations:
                        if mut.gene.name == gene.name:
                            if mut.signif in ['known', 'likely']:
                                actionable_var = ActionableVariant()
                                actionable_var.mut = mut.aa_change if mut.aa_change else '.'
                                actionable_var.var_type = mut.var_type
                                actionable_var.freq = mut.freq
                                actionable_var.status = mut.signif
                                actionable_variants.append(actionable_var)

            if cnv_mutation_types:
                for se in gene.seq2c_events:
                    if 'Amplification' in possible_mutation_types and se.amp_del == 'Amp' or \
                            'Deletion' in possible_mutation_types and se.amp_del == 'Del':
                        actionable_var = ActionableVariant()
                        actionable_var.mut = se.amp_del + ', ' + se.fragment
                        actionable_var.var_type = se.amp_del
                        actionable_variants.append(actionable_var)

            if sv_mutation_types:
                svs_by_key = OrderedDict()
                for sv in gene.sv_events:
                    if any(a.known or a.effect == 'EXON_DEL' for a in sv.annotations):
                        svs_by_key[sv.type, tuple(tuple(sorted(a.genes)) for a in sv.annotations)] = sv

                for se in svs_by_key.values():
                    if ('Fusion' in possible_mutation_types or 'Rearrangement' in possible_mutation_types) and se.type == 'BND' or \
                       'Deletion' in possible_mutation_types and se.type == 'DEL' or \
                       'Amplification' in possible_mutation_types and se.type == 'DUP':
                        actionable_var = ActionableVariant()
                        actionable_var.mut = ', '.join(set('/'.join(set(a.genes)) for a in se.key_annotations if a.genes))
                        actionable_var.var_type = BaseClinicalReporting.sv_type_dict.get(se.type, se.type)
                        actionable_variants.append(actionable_var)

            if not actionable_variants:
                continue

            def get_signif_order(status, freq):
                if status == 'known' and freq >= self.cnf.variant_filtering.act_min_freq:
                    return 2
                if status == 'known' or status == 'likely':
                    return 1
                return 0

            sorted_variants = sorted(actionable_variants, reverse=True,
                                     key=lambda x: (get_signif_order(x.status, x.freq), x.freq))
            for i, actionable_var in enumerate(sorted_variants):
                reg = report.add_row()
                gene_name = gene.name if i == 0 else ''
                rationale = actionable_genes_dict[gene.name][0] if i == 0 else ''
                types_alterations = actionable_genes_dict[gene.name][1].replace('; ', '\n') if i == 0 else ''
                therapeutic_agents = actionable_genes_dict[gene.name][2] if i == 0 else ''
                reg.add_record('Gene', gene_name, rowspan=str(len(sorted_variants)))
                reg.add_record('Variant', actionable_var.mut)
                reg.add_record('Type', actionable_var.var_type)
                reg.add_record('Types of recurrent alterations', types_alterations, rowspan=str(len(sorted_variants)))
                reg.add_record('Rationale', rationale, rowspan=str(len(sorted_variants)))
                reg.add_record('Therapeutic Agents', therapeutic_agents, rowspan=str(len(sorted_variants)))
                reg.add_record('Freq', actionable_var.freq)
                reg.add_record('VarDict status', actionable_var.status)

        return report


def make_circos_plot(cnf, output_fpath):
    info('Making circos plot')
    circos_py_executable = get_script_cmdline(cnf, 'python', join('scripts', 'post', 'circos.py'))
    cmdline = '{circos_py_executable} '

    required_files = [cnf.mutations_fpath, cnf.seq2c_tsv_fpath, cnf.sv_vcf_fpath]
    for fpath, desc in zip(required_files, ['Vardict results', 'Seq2C results', 'SV calling results']):
        if not fpath or not verify_file(fpath):
            warn('File with ' + desc + ' is not found (' + str(fpath) + '). Circos plot cannot be created.')
            return None

    vardict2mut_raw_fpath = cnf.mutations_fpath.replace(source.mut_pass_suffix + '.', '')
    output_dirpath = dirname(output_fpath)

    if cnf.bed_fpath:
        cmdline += ' --bed ' + cnf.bed_fpath
    cmdline += ' --mutations ' + vardict2mut_raw_fpath
    cmdline += ' --seq2c ' + cnf.seq2c_tsv_fpath
    cmdline += ' --sv ' + cnf.sv_vcf_fpath

    cmdline += ' --sample ' + cnf.sample
    cmdline += ' --genome ' + cnf.genome.name
    cmdline += ' -o ' + output_dirpath

    cmdline = cmdline.format(**locals())
    circos_plot_fpath = join(output_dirpath, cnf.sample + '.png')
    res = call(cnf, cmdline, stdout_to_outputfile=False, output_fpath=circos_plot_fpath, exit_on_error=False)
    if not verify_file(circos_plot_fpath):
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

