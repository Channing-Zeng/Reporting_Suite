from collections import OrderedDict
from json import load, dump, JSONEncoder, dumps

from source import verify_file, info
from source.reporting import MetricStorage, Metric, PerRegionSampleReport, ReportSection, SampleReport
from source.targetcov.flag_regions import get_depth_cutoff
from source.targetcov.summarize_targetcov import get_float_val, get_val


def make_key_genes_reports(cnf, sample, key_gene_names):
    info('Preparing coverage stats key gene tables')

    ave_depth = get_ave_coverage(sample, sample.targetcov_json_fpath)

    depth_cutoff = get_depth_cutoff(ave_depth, cnf.coverage_reports.depth_thresholds)

    stats_by_genename = dict()
    with open(sample.targetcov_detailed_tsv) as f_inp:
        for l in f_inp:
            if not l.startswith('#') and ('Whole-Gene' in l or 'Gene-Exon' in l):
                fs = l.split('\t')
                gene_name = get_val(fs[4])
                if gene_name in key_gene_names:
                    for t, field in zip(cnf.coverage_reports.depth_thresholds, fs[12:]):
                        if int(t) == depth_cutoff:
                            stats_by_genename[gene_name] = get_float_val(fs[9]), get_float_val(field)
                            continue

    clinical_cov_metric_storage = MetricStorage(
        sections=[ReportSection(metrics=[
            Metric('Gene'),
            Metric('Mean coverage'),
            Metric('% covered at {}x'.format(depth_cutoff), unit='%')])])
    key_genes_report = PerRegionSampleReport(sample=sample, metric_storage=clinical_cov_metric_storage)
    for gene_name in key_gene_names:
        ave_cov, depth_in_thresh = stats_by_genename.get(gene_name, (None, None))
        reg = key_genes_report.add_region()
        reg.add_record('Gene', gene_name)
        reg.add_record('Mean coverage', ave_cov)
        reg.add_record('% covered at {}x'.format(depth_cutoff), depth_in_thresh)

    key_genes_report.save_tsv(sample.clinical_targqc_tsv, human_readable=True)
    info('Saved coverage report to ' + key_genes_report.tsv_fpath)
    info('-' * 70)
    info()
    return key_genes_report


def get_ave_coverage(sample, targqc_json_fpath):
    with open(targqc_json_fpath) as f:
        data = load(f, object_pairs_hook=OrderedDict)
    sr = SampleReport.load(data, sample, None)
    r = sr.find_record(sr.records, 'Average target coverage depth')
    if not r:
        r = sr.find_record(sr.records, 'Average genome coverage depth')
    return r.value


def make_mutations_report(cnf, sample, key_gene_names, mutations_fpath):
    info('Preparing mutations stats for key gene tables')

    clinical_mut_metric_storage = MetricStorage(
        sections=[ReportSection(metrics=[
            Metric('Gene & Transcript'),  # Gene & Transcript
            Metric('Variant'),            # c.244G>A, p.Glu82Lys
            Metric('Allele'),             # Het.
            Metric('Genomic Pos.'),       # hg19 chr11:g.47364249G>A
            Metric('Depth'),              # 658
            Metric('Frequency'),          # .19
            Metric('AA length'),          # 128
            Metric('ID'),                 # rs352343, COSM2123
            Metric('Type'),               # Frameshift
            Metric('Classification'),     # Likely Pathogenic
        ])])
    report = PerRegionSampleReport(sample=sample, metric_storage=clinical_mut_metric_storage)

    info('Reading mutations from ' + mutations_fpath)
    with open(mutations_fpath) as f:
        for i, l in enumerate(f):
            if i == 0:
                continue
            fs = l.strip().split('\t')
            if len(fs) > 60:
                sample_name, chrom, start, id, ref, alt, type_, effect, func, codon_change, aa_change, cdna_change, \
                    aa_len, gene, transcr_biotype, coding, transcript, exon, cosmic_cds_change, cosmic_aa_change, \
                    cosmic_cnt, end, depth, af, bias, pmean, pstd, qual, qstd, sbf, gmaf, vd, clnsif, oddratio, hiaf, \
                    mq, sn, adjaf, nm, shift3, msi, dbsnpbuildid, vtype, status, paired_pval, paired_oddratiom, \
                    m_depth, m_af, m_vd, m_bias, m_pmean, m_pstd, m_qual, m_qstd, m_hiaf, m_mq, m_sn, m_adjaf, m_nm, \
                    n_sample, n_var, pcnt_sample, ave_af, filter_, var_type, var_class, status = fs[:67]  # 67 of them
            else:
                sample_name, chrom, start, id, ref, alt, type_, effect, func, codon_change, aa_change, cdna_change, \
                    aa_len, gene, transcr_biotype, coding, transcript, exon, cosmic_cds_change, cosmic_aa_change, \
                    cosmic_cnt, end, depth, af, bias, pmean, pstd, qual, qstd, sbf, gmaf, vd, clnsif, oddratio, hiaf, \
                    mq, sn, adjaf, nm, shift3, msi, dbsnpbuildid, \
                    n_sample, n_var, pcnt_sample, ave_af, filter_, var_type, var_class, status = fs[:50]  # 50 of them

            if sample_name == sample.name:  # and gene in key_gene_names:
                reg = report.add_region()
                reg.add_record('Gene & Transcript', gene + ' ' + transcript)
                reg.add_record('Variant', codon_change + (' p.' + aa_change if aa_change else ''))
                reg.add_record('Allele', None)
                reg.add_record('Genomic Pos.', str(cnf.genome) + ' ' + chrom + ':g.' +
                                   (Metric.format_value(int(start), human_readable=True) if start else '') +
                                   ' ' + ref + '>' + alt)
                reg.add_record('Depth', depth)
                reg.add_record('Frequency', af)
                reg.add_record('AA length', aa_len)
                reg.add_record('ID', id)
                reg.add_record('Type', type_)
                if status == 'likely':
                    status += ' pathogenic'
                reg.add_record('Classification', status)

    report.save_tsv(sample.clinical_mutation_tsv, human_readable=True)
    info('Saved mutations report to ' + report.tsv_fpath)
    info('-' * 70)
    info()
    return report
