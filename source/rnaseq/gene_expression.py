from collections import defaultdict
from os.path import join, dirname, abspath, basename


import source
from source import info, verify_file
from source.file_utils import file_transaction, adjust_system_path
from source.logger import err
from source.reporting.reporting import MetricStorage, Metric, PerRegionSampleReport, ReportSection, BaseReport

from ngs_reporting.utils import get_key_genes


class Counts:
    def __init__(self, name, chrom=None, gene_name=None, gene_id=None, transcript_id=None, counts=None, is_hidden_row=False):
        self.name = name
        self.chrom = chrom
        self.gene_name = gene_name
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.counts = counts
        self.is_hidden_row = is_hidden_row

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash(self.name)


def make_gene_expression_heatmaps(cnf, bcbio_structure, counts_fpath, genes_dict, report_fpath, report_name=None,
                                  not_rename_genes=False):
    key_gene_names = get_key_genes(verify_file(adjust_system_path(cnf.key_genes), 'key genes'))
    if not verify_file(counts_fpath):
        annotate_gene_counts(cnf, counts_fpath, genes_dict, report_name)

    gene_counts, samples = parse_gene_counts(counts_fpath, key_gene_names, report_name[0].lower() + report_name[1:], not_rename_genes)
    report = make_expression_heatmap(bcbio_structure, gene_counts)
    counts_fname = basename(counts_fpath)
    counts_link = '<a href="{counts_fpath}" target="_blank">{counts_fname}</a>'.format(**locals())
    data_dict = {
        'file_link': counts_link
    }
    BaseReport.save_html(report, cnf, report_fpath, caption=report_name,
                         extra_js_fpaths=[join(dirname(abspath(__file__)), 'static', 'rnaseq_heatmaps.js')],
                         extra_css_fpaths=[join(dirname(abspath(__file__)), 'static', 'rnaseq.css')],
                         tmpl_fpath=join(dirname(abspath(__file__)), 'template.html'), data_dict=data_dict)
    info(report_name + ' heatmap saved in ' + report_fpath)
    return


def annotate_gene_counts(cnf, counts_fpath, genes_dict, report_name):
    info('Annotating ' + report_name + ' from ' + counts_fpath)
    unannotated_fpath = counts_fpath.replace('annotated_', '')
    with file_transaction(cnf.work_dir, counts_fpath) as tx:
        with open(tx, 'w') as annotated_f:
            with open(unannotated_fpath) as f:
                for i, l in enumerate(f):
                    if i == 0:
                        header = l.replace('\n', '').split('\t')
                        l = '\t'.join(header + ['symbol'])
                        annotated_f.write(l + '\n')
                        continue
                    fs = l.replace('\n', '').split('\t')
                    gene_and_exon = fs[0].split(':')
                    gene_id = gene_and_exon[0]
                    if gene_id not in genes_dict:
                        continue
                    gene_symbol = genes_dict[gene_id]
                    l = '\t'.join(fs + [gene_symbol])
                    annotated_f.write(l + '\n')


def parse_gene_counts(counts_fpath, key_gene_names, report_name, not_rename_genes):
    gene_counts = defaultdict(list)
    info('Preparing ' + report_name + ' stats for expression heatmaps')
    info('Checking ' + counts_fpath)
    if not verify_file(counts_fpath):
        err('Cannot find ' + report_name + ' fpath')
        return []

    info('Reading ' + report_name + ' from ' + counts_fpath)
    samples_cols = dict()
    samples = []
    gene_col = None

    with open(counts_fpath) as f:
        for i, l in enumerate(f):
            if i == 0:
                header = l.strip().split('\t')
                gene_col = header.index('symbol')
                samples = header[1:gene_col]
                samples_cols = {sample: col + 1 for col, sample in enumerate(samples)}
                continue
            fs = l.replace('\n', '').split('\t')
            gene_name = fs[gene_col]
            if gene_name not in key_gene_names:
                continue
            is_hidden_row = False
            name = gene_name
            if ':' in fs[0]:  ## exon number
                is_hidden_row = True
                exon_number = fs[0].split(':')[1]
                name += ':' + exon_number
            gene_expression_dict = {sample: int(float(fs[col])) if float(fs[col]).is_integer() else float(fs[col])
                                    for sample, col in samples_cols.iteritems()}
            if not_rename_genes:
                is_hidden_row = True
                name = fs[0]  # use id
            gene = Counts(name, gene_name=gene_name, counts=gene_expression_dict, is_hidden_row=is_hidden_row)
            gene_counts[gene_name].append(gene)

    return gene_counts, samples


def make_expression_heatmap(bcbio_structure, gene_counts):
    samples_names = [sample.name for sample in bcbio_structure.samples]
    metrics = [Metric('Gene')] + [Metric(sample_name, max_width=40) for sample_name in samples_names]
    metric_storage = MetricStorage(sections=[ReportSection(metrics=metrics, name='samples')])
    report = PerRegionSampleReport(metric_storage=metric_storage, expandable=True, unique=True, heatmap_by_rows=True,
                                   keep_order=True, large_table=True, vertical_sample_names=True)
    printed_genes = set()

    # Writing records
    for gene in sorted(gene_counts.keys()):
        first_record = gene_counts[gene][0]
        if first_record.is_hidden_row:
            printed_genes.add(first_record.gene_name)
            row = report.add_row()
            row.add_record('Gene', first_record.gene_name)
            for sample in samples_names:
                row.add_record(sample, sum([record.counts[sample] for record in gene_counts[gene]]))
            row.class_ = ' expandable_gene_row collapsed'
        for record in gene_counts[gene]:
            gene_expression = record.counts
            row = report.add_row()
            if record.is_hidden_row:
                row_class = ' row_to_hide row_hidden'
            else:
                row_class = ' expandable_gene_row collapsed'
            row.add_record('Gene', record.name)
            for sample, count in gene_expression.iteritems():
                if sample in samples_names:
                    row.add_record(sample, count)
            row.class_ = row_class
    return report

