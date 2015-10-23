from collections import OrderedDict
from itertools import izip
from json import load
import json
from os.path import join, dirname, abspath

import source
from source import verify_file, info
from source.clinical_reporting import seq2c_plot
from source.file_utils import add_suffix, verify_module
from source.logger import warn, err
from source.reporting.reporting import MetricStorage, Metric, PerRegionSampleReport, ReportSection, SampleReport, \
    calc_cell_contents, make_cell_td, write_static_html_report, make_cell_th, build_report_html
from source.targetcov.flag_regions import get_depth_cutoff
from source.targetcov.summarize_targetcov import get_float_val, get_val
from source.utils import get_chr_lengths


# def make_clinical_reports_bcbio(cnf, bcbio_structure):
#     for sample in bcbio_structure.samples:
#         info('Building clinical report for AZ 300 key genes ' + str(cnf.key_genes))
#         cnf.work_dir = join(bcbio_structure.work_dir, '{sample}_' + source.clinreport_name)
#
#         info('Building clinical report for AZ 300 key genes ' + str(cnf.key_genes))
#         cnf.key_genes = verify_file(cnf.key_genes, is_critical=True, description='300 AZ key genes')
#         with open(cnf.key_genes) as f:
#             key_gene_names = set([l.strip() for l in f.readlines() if l.strip() != ''])
#
#
#         make_sample_clinical_report(cnf, sample)
#
#         cnf.work_dir = bcbio_structure.work_dir


def run_sample_clinical_reporting(cnf):
    cr = ClinicalReporting(cnf)
    return cr.make()


class ClinicalReporting:
    ACTIONABLE_GENES_FPATH = join(__file__, '..', 'db', 'broad_db.tsv')

    def __init__(self, cnf):
        self.cnf = cnf
        self.sample = source.BaseSample(cnf.sample, cnf.output_dir, targqc_dirpath=cnf.targqc_dirpath, clinical_report_dirpath=cnf.output_dir)
        self.key_gene_names = None
        self.depth_cutoff = None

    def make(self):
        info('Building clinical report for AZ 300 key genes ' + str(self.cnf.key_genes) + ', sample ' + self.sample.name)
        self.key_gene_names = get_key_genes(self.cnf.key_genes)

        ave_depth = get_ave_coverage(self.sample, self.sample.targetcov_json_fpath)
        self.depth_cutoff = get_depth_cutoff(ave_depth, self.cnf.coverage_reports.depth_thresholds)

        target_fraction = get_target_fraction(self.sample, self.sample.targetcov_json_fpath)
        gender = get_gender(self.sample, self.sample.targetcov_json_fpath)
        total_variants = get_total_variants_number(self.sample, self.cnf.varqc_json_fpath)

        mutations = self.parse_mutations(self.cnf.mutations_fpath)

        mutations_report = self.make_mutations_report(mutations)


        mut_info_by_gene = dict()
        if gene not in mut_info_by_gene:
            mut_info_by_gene[gene] = []
        mut_infos = ['p.' + aa_change if aa_change else '', mut_type]
        mut_info_by_gene[gene].append([mutation for mutation in mut_infos if mutation])


        key_gene_by_name = self.parse_targetseq_detailed()

        key_genes_report = self.make_key_genes_cov_report(key_gene_by_name, ave_depth)

        actionable_genes_report = self.proc_actionable_genes(mutations_report, seq2c_data_by_genename)

        seq2c_events_by_gene_name, seq2c_plot_fpath = None, None
        if not self.cnf.seq2c_tsv_fpath:
            warn('No Seq2C results provided by option --seq2c, skipping plotting Seq2C')
        else:
            seq2c_events_by_gene_name = self.parse_seq2c_report(self.cnf.seq2c_tsv_fpath)
            for gene_name, gene in key_gene_by_name.items():
                seq2c_events_by_gene_name[gene_name].seq2c = 

            if not verify_module('matplotlib'):
                warn('No matplotlib, skipping plotting Seq2C')
            else:
                seq2c_plot_fpath = seq2c_plot.draw_seq2c_plot(self.cnf, self.cnf.seq2c_tsv_fpath,
                    self.sample.name, self.cnf.output_dir, self.key_gene_names)


        cov_plot_data = self.save_key_genes_cov_json(key_gene_by_name, mutations, seq2c_data_by_genename)



        if seq2c_data_by_genename:
            clinical_cov_metrics.append(Metric('CNV', short_name='&nbsp;&nbsp;CNV'))  # short name is hack for IE9 who doesn't have "text-align: left" and tries to stick "CNV" to the previous col header

        sample_dict = dict()
        sample_dict['sample'] = sample.name.replace('_', ' ')
        sample_dict['project_name'] = cnf.project_name.replace('_', ' ')
        sample_dict['panel'] = cnf.target_type
        sample_dict['bed_path'] = ''
        if cnf.debug:
            sample_dict['panel'] = cnf.target_type + ', AZ300 IDT panel'
            sample_dict['bed_path'] = 'http://blue.usbod.astrazeneca.net/~klpf990/reference_data/genomes/Hsapiens/hg19/bed/Panel-IDT_PanCancer_AZSpike_V1.bed'

        sample_dict['sex'] = 'Undetermined'
        if gender:
            if gender.startswith('M'):
                sample_dict['sex'] = 'Male'
            elif gender.startswith('F'):
                sample_dict['sex'] = 'Female'

        sample_dict['sample_type'] = 'unpaired'
        if cnf.debug:
            sample_dict['sample_type'] = 'plasma, unpaired'

        sample_dict['genome_build'] =  cnf.genome.name
        sample_dict['target_type'] = target_type
        sample_dict['target_fraction'] = Metric.format_value(target_fraction, is_html=True, unit='%')
        # approach_dict['min_depth'] = Metric.format_value(min_depth, is_html=True)
        sample_dict['ave_depth'] = Metric.format_value(ave_depth, is_html=True)

        mutations_dict = dict()
        if mutations_report.regions:
            # if cnf.debug:
            #     mutations_report.regions = mutations_report.regions[::20]
            mutations_dict['table'] = build_report_html(mutations_report, sortable=True)
            mutations_dict['total_variants'] = Metric.format_value(total_variants, is_html=True)
            mutations_dict['total_key_genes'] = Metric.format_value(len(key_gene_names), is_html=True)

        coverage_dict = dict(depth_cutoff=depth_cutoff, columns=[])
        GENE_COL_NUM = 3
        genes_in_col = len(coverage_report.regions) / GENE_COL_NUM
        calc_cell_contents(coverage_report, coverage_report.get_rows_of_records())
        for i in range(GENE_COL_NUM):
            column_dict = dict()
            # column_dict['table'] = build_report_html(coverage_report)
            column_dict['metric_names'] = [make_cell_th(m) for m in coverage_report.metric_storage.get_metrics()]
            column_dict['rows'] = [
                dict(records=[make_cell_td(r) for r in region.records])
                    for region in coverage_report.regions[i * genes_in_col:(i+1) * genes_in_col]]
            coverage_dict['columns'].append(column_dict)
        coverage_dict['plot_data'] = cov_plot_data

        image_by_key = dict()
        if seq2c_plot_fpath:
            image_by_key['seq2c_plot'] = seq2c_plot_fpath

        actionable_genes_dict = dict()
        if actionable_genes_report.regions:
            actionable_genes_dict['table'] = build_report_html(actionable_genes_report, sortable=False)

        sample.clinical_html = write_static_html_report(cnf, {
            'sample': sample_dict,
            'variants': mutations_dict,
            'coverage': coverage_dict,
            'actionable_genes': actionable_genes_dict,
            'seq2c_plot': True,
        }, sample.clinical_html,
           tmpl_fpath=join(dirname(abspath(__file__)), 'template.html'),
           extra_js_fpaths=[join(dirname(abspath(__file__)), 'static', 'clinical_report.js'),
                            join(dirname(abspath(__file__)), 'static', 'draw_genes_coverage_plot.js')],
           extra_css_fpaths=[join(dirname(abspath(__file__)), 'static', 'clinical_report.css'),
                             join(dirname(abspath(__file__)), 'static', 'header_picture.css')],
           image_by_key=image_by_key)

        # clin_rep_symlink = adjust_path(join(sample.dirpath, '..', sample.name + '.clinical_report.html'))
        # if islink(clin_rep_symlink):
        #     os.unlink(clin_rep_symlink)
        # os.symlink(sample.clinical_html, clin_rep_symlink)

        info('Saved clinical report to ' + sample.clinical_html)
        info('-' * 70)
        info()
        return sample.clinical_html


    class Seq2CEvent:
        def __init__(self):
            self.gene_name = None
            self.type = None  # Del, Amp
            self.fragment = None  # Whole, BP

    def parse_seq2c_report(self, seq2c_tsv_fpath):
        seq2c_events_by_gene_name = dict()
        if not verify_file(seq2c_tsv_fpath, silent=True):
            return None

        with open(seq2c_tsv_fpath) as f_inp:
            for i, l in enumerate(f_inp):
                if i == 0:
                    continue
                fs = l.strip().split('\t')
                sample_name = fs[0]
                gene_name = fs[1]
                fragment = fs[8]
                type_ = fs[9]

                if sample_name == self.sample.name and gene_name in self.key_gene_names and type_ in ['Del', 'Amp']:
                    event = ClinicalReporting.Seq2CEvent()
                    event.gene_name = gene_name
                    event.type = type_
                    event.fragment = fragment
                    seq2c_events_by_gene_name[gene_name] = event

        return seq2c_events_by_gene_name


    class KeyGene:
        def __init__(self, name, ave_depth=None):
            self.name = name
            self.chrom = None
            self.start = None
            self.end = None
            self.cdss = []
            self.ave_depth = ave_depth
            self.cov_by_threshs = []

    class CDS:
        def __init__(self):
            self.start = None
            self.end = None
            self.ave_depth = None
            self.cov_by_threshs = []

    def parse_targetseq_detailed(self):
        info('Preparing coverage stats key gene tables')

        key_gene_by_name = dict()

        with open(self.sample.targetcov_detailed_tsv) as f_inp:
            for l in f_inp:
                if l.startswith('#'):
                    continue

                fs = l.split('\t')  #Chr	Start	End	Size	Gene	Strand	Feature	Biotype	Min depth	Ave depth	Std dev	W/n 20% of ave depth	1x	5x	10x	25x	50x	100x	500x	1000x	5000x	10000x	50000x
                chrom, start, end, size, symbol, strand, feature, biotype, min_depth, ave_depth, std_dev, wn20pcnt = fs[:12]
                pcnt_val_by_thresh = fs[12:]
                symbol = get_val(symbol)
                if symbol not in self.key_gene_names:
                    continue

                chrom = get_val(chrom)
                if not chrom:
                    err('For gene ' + str(symbol) + ', chrom value if empty in line ' + l + ' in ' + self.sample.targetcov_detailed_tsv)
                    continue

                start, end = int(start), int(end)
                ave_depth = get_float_val(ave_depth)
                cov_by_threshs = [get_float_val(f) for t, f in izip(self.cnf.coverage_reports.depth_thresholds, pcnt_val_by_thresh)]

                if feature in ['Whole-Gene', 'Gene-Exon']:
                    gene = key_gene_by_name.get(symbol)
                    if not gene:
                        gene = ClinicalReporting.KeyGene(symbol)
                        key_gene_by_name[symbol] = gene
                    gene.ave_depth = ave_depth
                    gene.chrom = chrom
                    gene.start = start
                    gene.end = end
                    gene.cov_by_threshs = cov_by_threshs
                    key_gene_by_name[symbol] = gene

                if feature in ['CDS', 'Exon']:
                    cds = ClinicalReporting.CDS()
                    cds.start = start
                    cds.end = end
                    cds.ave_depth = ave_depth
                    cds.cov_by_threshs = [get_float_val(f) for t, f in izip(self.cnf.coverage_reports.depth_thresholds, fs[12:])]
                    gene = key_gene_by_name.get(symbol)
                    if not gene:
                        gene = ClinicalReporting.KeyGene(symbol)
                        key_gene_by_name[symbol] = gene
                    gene.cdss.append(cds)

        return key_gene_by_name

# gene_pos = int(get_float_val(fs[1]) + get_float_val(fs[2])) / 2
    def make_key_genes_cov_report(self, ave_depth, key_gene_by_name):
        clinical_cov_metrics = [
            Metric('Gene'),
            Metric('Chr', with_heatmap=False, max_width=20, align='right'),
            Metric('Ave depth', med=ave_depth),
            Metric('% cov at {}x'.format(self.depth_cutoff), unit='%', med=1, low_inner_fence=0.5, low_outer_fence=0.1)]

        clinical_cov_metric_storage = MetricStorage(sections=[ReportSection(metrics=clinical_cov_metrics)])

        key_genes_report = PerRegionSampleReport(sample=self.sample, metric_storage=clinical_cov_metric_storage)

        for gene in sorted(key_gene_by_name.values(), key=lambda g: g.name):
            reg = key_genes_report.add_region()
            reg.add_record('Gene', gene.name)
            reg.add_record('Chr', gene.chrom.replace('chr', ''))
            reg.add_record('Ave depth', gene.ave_depth)
            m = clinical_cov_metric_storage.find_metric('% cov at {}x'.format(self.depth_cutoff))
            reg.add_record(m.name, next((cov for cutoff, cov in gene.cov_by_threshs.items() if cutoff == self.depth_cutoff), None))
            if gene.seq2c:
                reg.add_record('CNV', gene.seq2c)

        key_genes_report.save_tsv(self.sample.clinical_targqc_tsv, human_readable=True)
        info('Saved coverage report to ' + key_genes_report.tsv_fpath)
        info('-' * 70)
        info()
        return key_genes_report


    class Mutation:
        def __init__(self):
            self.gene = None
            self.transcript = None
            self.codon_change = None
            self.chrom = None
            self.pos = None
            self.ref = None
            self.alt = None
            self.depth = None
            self.freq = None
            self.aa_change = None
            self.aa_len = None
            self.type = None
            self.status = None
            self.cosmic_ids = []
            self.dbsnp_ids = []

    def make_mutations_report(self, mutations):
        max_width = '90'
        clinical_mut_metric_storage = MetricStorage(
            sections=[ReportSection(metrics=[
                Metric('Gene'),  # Gene & Transcript
                Metric('Transcript'),  # Gene & Transcript
                Metric('Codon chg', max_width=max_width, class_='long_line'),            # c.244G>A
                Metric('AA chg', max_width=max_width, class_='long_line'),            # p.Glu82Lys
                # Metric('Allele'),             # Het.
                # Metric('Chr', max_width=33, with_heatmap=False),       # chr11
                Metric('Position', sort_by=lambda v: (v.split(':')[0], int(''.join(ch for ch in v.split(':')[1] if ch.isdigit()))),
                       align='left'),       # g.47364249
                Metric('Change', max_width=max_width, class_='long_line'),       # G>A
                Metric('Depth', max_width=48),              # 658
                Metric('Freq', max_width=45, unit='%', with_heatmap=False),          # .19
                Metric('AA len', max_width=50, with_heatmap=False),          # 128
                Metric('DB', max_width=90, class_='long_line'),                 # rs352343, COSM2123
                # Metric('COSMIC', max_width=70, style='', class_='long_line'),                 # rs352343, COSM2123
                Metric('Type', max_width='100', class_='long_line'),               # Frameshift
                Metric('Status'),     # Likely Pathogenic
            ])])
        report = PerRegionSampleReport(sample=self.sample, metric_storage=clinical_mut_metric_storage,
            hide_rows_by=lambda r: r.metric.name == 'Status' and r.value and r.value.lower() == 'unknown',
            highlight_by=lambda r: r.metric.name == 'Status' and r.value and r.value.lower() == 'known')

        for mut in mutations:
            reg = report.add_region()
            reg.add_record('Gene', mut.gene)
            reg.add_record('Transcript', mut.transcript)
            reg.add_record('Codon chg', mut.codon_change)
            reg.add_record('AA chg', 'p.' + mut.aa_change if mut.aa_change else '')
            # reg.add_record('Allele', allele_record)
            # reg.add_record('Chr', chrom.replace('chr', '') if chrom else '')
            c = mut.chrom.replace('chr', '') if mut.chrom else ''
            p = 'g.' + Metric.format_value(mut.start, human_readable=True) if mut.start else ''
            reg.add_record('Position', c + ':' + p)
            reg.add_record('Change', mut.ref + '>' + mut.alt)
            reg.add_record('Depth', mut.depth)
            reg.add_record('Freq', mut.freq)
            reg.add_record('AA len', mut.aa_len)
            db = None
            if mut.dbsnp_ids:
                db += ', '.join(('<a href="http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs=' + rs_id + '">dbSNP</a>') for rs_id in mut.dbsnp_ids)
            if mut.dbsnp_ids and mut.cosmic_ids:
                db += ', '
            if mut.cosmic_ids:
                db += ', '.join('<a href="http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=' + cid + '">COSM</a>' for cid in mut.cosmic_ids)
            reg.add_record('DB', db, parse=False)
            reg.add_record('Type', mut.type)
            reg.add_record('Status', mut.status.lower() if mut.status else mut.status)

        report.save_tsv(self.sample.clinical_mutation_tsv, human_readable=True)
        info('Saved mutations report to ' + report.tsv_fpath)
        info('-' * 70)
        info()
        return report

    def parse_mutations(self, mutations_fpath):
        mutations = []

        info('Preparing mutations stats for key gene tables')
        mutations_fpath = None
        if not verify_file(mutations_fpath, silent=True):
            single_mutations_fpath = add_suffix(mutations_fpath, source.mut_single_suffix)
            paired_mutations_fpath = add_suffix(mutations_fpath, source.mut_paired_suffix)
            if verify_file(single_mutations_fpath, silent=True) and is_sample_presents_in_file(self.sample.name, single_mutations_fpath):
                mutations_fpath = single_mutations_fpath
            elif verify_file(paired_mutations_fpath, silent=True):
                mutations_fpath = paired_mutations_fpath
            else:
                info('Cannot find PASSed mutations fpath')
                return None

        info('Reading mutations from ' + mutations_fpath)
        alts_met_before = set()
        with open(mutations_fpath) as f:
            for i, l in enumerate(f):
                if i == 0:
                    continue
                fs = l.strip().split('\t')
                if len(fs) > 60:
                    sample_name, chrom, start, ids, ref, alt, type_, effect, func, codon_change, aa_change, cdna_change, \
                        aa_len, gene, transcr_biotype, coding, transcript, exon, cosmic_cds_change, cosmic_aa_change, \
                        cosmic_cnt, end, depth, af, bias, pmean, pstd, qual, qstd, sbf, gmaf, vd, clnsif, oddratio, hiaf, \
                        mq, sn, adjaf, nm, shift3, msi, dbsnpbuildid, vtype, status1, paired_pval, paired_oddratiom, \
                        m_depth, m_af, m_vd, m_bias, m_pmean, m_pstd, m_qual, m_qstd, m_hiaf, m_mq, m_sn, m_adjaf, m_nm, \
                        n_sample, n_var, pcnt_sample, ave_af, filter_, var_type, var_class, status = fs[:67]  # 67 of them
                else:
                    sample_name, chrom, start, ids, ref, alt, type_, effect, func, codon_change, aa_change, cdna_change, \
                        aa_len, gene, transcr_biotype, coding, transcript, exon, cosmic_cds_change, cosmic_aa_change, \
                        cosmic_cnt, end, depth, af, bias, pmean, pstd, qual, qstd, sbf, gmaf, vd, clnsif, oddratio, hiaf, \
                        mq, sn, adjaf, nm, shift3, msi, dbsnpbuildid, \
                        n_sample, n_var, pcnt_sample, ave_af, filter_, var_type, var_class, status = fs[:50]  # 50 of them

                if sample_name == self.sample.name and gene in self.key_gene_names:
                    if (chrom, start, ref, alt) in alts_met_before:
                        continue
                    alts_met_before.add((chrom, start, ref, alt))

                    mut = ClinicalReporting.Mutation()
                    mut.gene = gene
                    mut.transcript = transcript
                    mut.codon_change = codon_change
                    mut.aa_change = aa_change
                    mut.aa_len = aa_len
                    mut.chrom = chrom
                    mut.start = start
                    mut.ref = ref
                    mut.alt = alt
                    mut.depth = depth
                    mut.freq = af
                    mut.dbsnp_ids = [''.join(c for c in id_ if c.isdigit()) for id_ in ids.split(';') if id_.startswith('rs')]
                    mut.cosmic_ids = [''.join(c for c in id_ if c.isdigit()) for id_ in ids.split(';') if id_.startswith('COS')]
                    mut.type = (type_[0] + type_[1:].lower().replace('_', ' ')) if type_ else type_
                    mut.status = status

                    mutations.append(mut)
        return mutations


    def proc_actionable_genes(self, mutations_report, seq2c_data_by_genename):
        act_fpath = verify_file(ClinicalReporting.ACTIONABLE_GENES_FPATH, is_critical=False, description='Actionable genes')
        if act_fpath:
            with open(act_fpath) as f:
                actionable_gene_table = [l.split('\t') for l in f.readlines()]
                actionable_gene_dict = dict((l[0], l[1:]) for l in actionable_gene_table)

            return self.make_actionable_genes_report(actionable_gene_dict, mutations_report, seq2c_data_by_genename)
        return None


    def make_actionable_genes_report(cnf, sample, key_gene_names, actionable_genes, mutations_report, seq2c_data_by_genename):
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
        report = PerRegionSampleReport(sample=sample, metric_storage=clinical_action_metric_storage)
        actionable_gene_names = actionable_genes.keys()
        for gene in actionable_gene_names:
            if gene not in key_gene_names:
                continue
            possible_mutations = actionable_genes[gene][1].split('; ')
            skipped_mutations = ['Rearrangement', 'Fusion']
            possible_mutations = [mutation for mutation in possible_mutations if mutation not in skipped_mutations]
            if not possible_mutations:
                continue
            variants = []
            types = []
            amp_del = None
            for region in mutations_report.regions:
                if mutations_report.find_record(region.records, 'Gene').value == gene:
                    variant = mutations_report.find_record(region.records, 'AA chg').value
                    variants.append(variant if variant else '.')
                    types.append(mutations_report.find_record(region.records, 'Type').value)
            if gene in seq2c_data_by_genename:
                amp_del = seq2c_data_by_genename[gene]

            if 'Amplification' in possible_mutations and (not amp_del or 'Amp' not in amp_del):
                possible_mutations.remove('Amplification')
            if 'Deletion' in possible_mutations and (not amp_del or 'Del' not in amp_del):
                possible_mutations.remove('Deletion')
            if not possible_mutations or (not variants and 'Amplification' not in possible_mutations and 'Deletion' not in possible_mutations):
                continue

            reg = report.add_region()
            reg.add_record('Gene', gene)
            reg.add_record('Variant', '\n'.join(variants))
            reg.add_record('Type', '\n'.join(types))
            reg.add_record('Types of recurrent alterations', actionable_genes[gene][1].replace('; ', '\n'))
            reg.add_record('Rationale', actionable_genes[gene][0])
            reg.add_record('Therapeutic Agents', actionable_genes[gene][2])

        report.save_tsv(sample.clinical_target_tsv, human_readable=True)
        info('Saved report for actionable genes to ' + report.tsv_fpath)
        info('-' * 70)
        info()
        return report


    def save_key_genes_cov_json(self, stats_by_genename, mut_info_by_gene, seq2c_data_by_genename, detailed_cov_by_gene):
        for gene in seq2c_data_by_genename.keys():
            if gene not in mut_info_by_gene:
                mut_info_by_gene[gene] = []
            mut_info_by_gene[gene].append([seq2c_data_by_genename[gene]])

        chr_names_lengths = OrderedDict((chr_, l) for chr_, l in (get_chr_lengths(self.cnf)).items()
                                        if '_' not in chr_)  # not drawing extra chromosomes chr1_blablabla
        chr_names = chr_names_lengths.keys()
        chr_short_names = [chrom[3:] for chrom in chr_names_lengths.keys()]
        chr_lengths = [chrom for chrom in chr_names_lengths.values()]
        chr_cum_lens = [sum(chr_lengths[:i]) for i in range(len(chr_lengths)+1)]
        ticks_x = [[(chr_cum_lens[i] + chr_cum_lens[i+1])/2, chr_short_names[i]] for i in range(len(chr_names))]

        key_genes = stats_by_genename.keys()
        gene_ave_depths = [stats_by_genename[gene_name][1] for gene_name in key_genes]
        cov_in_thresh = [stats_by_genename[gene_name][2] for gene_name in key_genes]
        coord_x = [chr_cum_lens[chr_names.index(stats_by_genename[gene_name][0])] + stats_by_genename[gene_name][3] for gene_name in key_genes]
        json_txt = '{"coord_x":%s, "coord_y":%s, "cov_in_thresh":%s, "gene_names":%s, "mutations":%s, "ticks_x":%s, ' \
                   '"lines_x":%s, "detailed_cov":%s}'\
                   % (json.dumps(coord_x), json.dumps(gene_ave_depths), json.dumps(cov_in_thresh), json.dumps(key_genes),
                       json.dumps(mut_info_by_gene), json.dumps(ticks_x), json.dumps(chr_cum_lens), json.dumps(detailed_cov_by_gene))
        return json_txt


def get_key_genes(key_genes_fpath):
    key_genes_fpath = verify_file(key_genes_fpath, is_critical=True, description='300 AZ key genes')
    with open(key_genes_fpath) as f:
        key_gene_names = set([l.strip() for l in f.readlines() if l.strip() != ''])

    return key_gene_names


def get_target_fraction(sample, targqc_json_fpath):
    with open(targqc_json_fpath) as f:
        data = load(f, object_pairs_hook=OrderedDict)
    sr = SampleReport.load(data, sample, None)
    r = sr.find_record(sr.records, 'Percentage of target covered by at least 1 read')
    if not r:
        r = sr.find_record(sr.records, 'Percentage of genome covered by at least 1 read')
    return r.value if r else None


def get_gender(sample, targqc_json_fpath):
    with open(targqc_json_fpath) as f:
        data = load(f, object_pairs_hook=OrderedDict)
    sr = SampleReport.load(data, sample, None)
    r = sr.find_record(sr.records, 'Gender')
    return r.value if r else None


def get_ave_coverage(sample, targqc_json_fpath):
    with open(targqc_json_fpath) as f:
        data = load(f, object_pairs_hook=OrderedDict)
    sr = SampleReport.load(data, sample, None)
    r = sr.find_record(sr.records, 'Average target coverage depth')
    if not r:
        r = sr.find_record(sr.records, 'Average genome coverage depth')
    return r.value if r else None


def get_min_coverage(sample, targqc_json_fpath):
    with open(targqc_json_fpath) as f:
        data = load(f, object_pairs_hook=OrderedDict)
    sr = SampleReport.load(data, sample, None)
    r = sr.find_record(sr.records, 'Minimum target coverage depth')
    if not r:
        r = sr.find_record(sr.records, 'Minimum genome coverage depth')
    return r.value if r else None


def get_total_variants_number(sample, varqc_json_fpath):
    with open(varqc_json_fpath) as f:
        data = load(f, object_pairs_hook=OrderedDict)
    sr = SampleReport.load(data, sample, None)
    r = sr.find_record(sr.records, 'Total with rejected')
    return r.value if r else None


def is_sample_presents_in_file(sample_name, mutations_fpath):
    with open(mutations_fpath) as f:
        for i, l in enumerate(f):
            if i == 0:
                continue
            fs = l.strip().split('\t')
            if fs[0] == sample_name:
                return True
    return False


def tooltip_long(string, max_len=30):
    if len(string) < max_len:
        return string
    else:
        return '<a class="tooltip-link" rel="tooltip" title="' + string + '">' + string[:max_len - 2] + '...</a>'



