from collections import OrderedDict
from json import load, dump, JSONEncoder, dumps

from source import verify_file, info
from source.reporting import MetricStorage, Metric, PerRegionSampleReport, ReportSection, SampleReport
from source.targetcov.flag_regions import get_depth_cutoff
from source.targetcov.summarize_targetcov import get_float_val, get_val


def make_key_genes_reports(cnf, sample):
    cnf.key_genes = verify_file(cnf.key_genes, is_critical=True)
    info('Generating tables for AZ300')

    with open(cnf.key_genes) as f:
        key_gene_names = set([l.strip() for l in f.readlines()])

    ave_depth = get_ave_coverage(sample, sample.targetcov_json_fpath)

    depth_cutoff = get_depth_cutoff(ave_depth, cnf.coverage_reports.depth_thresholds)

    ms = MetricStorage(
        sections=[ReportSection(metrics=[
            Metric('Gene'),
            Metric('Mean coverage'),
            Metric('% covered at {}x'.format(depth_cutoff), unit='%')])])
    key_genes_report = PerRegionSampleReport(sample=sample, metric_storage=ms)

    with open(sample.targetcov_detailed_tsv) as f_inp:
        for l in f_inp:
            if not l.startswith('#') and ('Whole-Gene' in l or 'Gene-Exon' in l):
                fs = l.split('\t')
                gene_name = get_val(fs[4])
                if gene_name in key_gene_names:
                    reg = key_genes_report.add_region()
                    reg.add_record('Gene', gene_name)
                    reg.add_record('Mean coverage', get_float_val(fs[9]))
                    for t, field in zip(cnf.coverage_reports.depth_thresholds, fs[12:]):
                        if int(t) == depth_cutoff:
                            value = get_float_val(field)
                            reg.add_record('% covered at {}x'.format(depth_cutoff), value)
                            continue

    key_genes_report.save_tsv(sample.clinical_targqc_tsv, human_readable=True)
    return key_genes_report


def get_ave_coverage(sample, targqc_json_fpath):
    with open(targqc_json_fpath) as f:
        data = load(f, object_pairs_hook=OrderedDict)
    sr = SampleReport.load(data, sample, None)
    r = sr.find_record(sr.records, 'Average target coverage depth')
    if not r:
        r = sr.find_record(sr.records, 'Average genome coverage depth')
    return r.value


def make_mutations_report(cnf, sample, gene_names, mutations_fpath):
    if not verify_file(cnf.key_genes, silent=True):
        return None

    info('Generating tables for AZ 300')
    with open(cnf.key_genes) as f:
        key_gene_names = set([l.strip() for l in f.readlines()])

    # grep sample.name in vardict_caller.single_mut_res_fpath
    # grep sample.name in vardict_caller.paired_mut_res_fpath
    # select green fields
    # make PerRegionSampleReport and save_tsv()