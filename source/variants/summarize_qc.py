from collections import OrderedDict
from source.logger import info
from source.reporting import summarize, write_summary_reports, get_sample_report_fpaths_for_bcbio_final_dir, Metric, \
    Record
from source.utils import OrderedDefaultDict

main_novelty = 'all'
metrics_header = 'Metric'
novelty_header = 'Novelty'
average_header = 'Average'


class VariantCaller:
    def __init__(self, suf):
        self.name = suf
        self.suf = suf
        self.single_qc_rep_fpaths = []
        self.summary_qc_report = None
        self.summary_qc_rep_fpaths = []


def make_summary_reports(cnf, sample_names):
    varqc_dir = cnf['base_name']

    vcf_sufs = cnf['vcf_suf'].split(',')
    callers = [VariantCaller(suf) for suf in vcf_sufs]

    for caller in callers:
        fpaths, sample_names = get_sample_report_fpaths_for_bcbio_final_dir(
            cnf['bcbio_final_dir'], sample_names, varqc_dir,
            '-' + caller.suf + '.varqc.txt')
        if fpaths:
            caller.single_qc_rep_fpaths = fpaths

    if len(callers) > 1:
        _make_for_multiple_variant_callers(callers, cnf, sample_names)

    else:
        _make_for_single_variant_caller(callers, cnf, sample_names)


def _make_for_single_variant_caller(callers, cnf, sample_names):
    full_report = summarize(sample_names, callers[0].single_qc_rep_fpaths, get_parse_qc_sample_report(cnf))

    full_summary_fpaths = write_summary_reports(
        cnf['output_dir'], cnf['work_dir'], full_report, 'varQC', 'Variant QC')

    info()
    info('*' * 70)
    for fpath in full_summary_fpaths:
        info(fpath)


def _make_for_multiple_variant_callers(callers, cnf, sample_names):
    for caller in callers:
        caller.summary_qc_report = summarize(
            sample_names, caller.single_qc_rep_fpaths, get_parse_qc_sample_report(cnf))

        caller.summary_qc_rep_fpaths = write_summary_reports(
            cnf['output_dir'], cnf['work_dir'], caller.summary_qc_report,
            caller.suf + '.varQC', 'Variant QC for ' + caller.name)

    all_single_reports = [r for c in callers for r in c.single_qc_rep_fpaths]
    all_sample_names = [sample_name + '-' + c.suf for sample_name in sample_names for c in callers]

    full_summary_report = summarize(all_sample_names, all_single_reports, get_parse_qc_sample_report(cnf))

    full_summary_fpaths = write_summary_reports(
        cnf['output_dir'], cnf['work_dir'], full_summary_report, 'varQC', 'Variant QC')

    info()
    info('*' * 70)

    for caller in callers:
        info(caller.name)
        for fpath in caller.summary_qc_rep_fpaths:
            info('  ' + fpath)
        info()

    info('Total')
    for fpath in full_summary_fpaths:
        info('  ' + fpath)


def __to_dict(metrics):
    return {m.name: m for m in metrics}

METRICS = __to_dict([
    Metric('nEvalVariants',   'total',       'Total variants evaluated'),
    Metric('nSNPs',           'SNP',         'SNPs'),
    Metric('nInsertions',     'ins',         'Insertions'),
    Metric('nDeletions',      'del',         'Deletions'),
    Metric('nVariantsAtComp', 'at comp',     'Number of eval sites at comp sites (that is, sharing the same locus as a variant in the comp track, regardless of whether the alternate allele is the same)'),
    Metric('compRate',        'comp rate',   'Percentage of eval sites at comp sites'),
    Metric('nConcordant',     'concord',     'Number of concordant sites (that is, for the sites that share the same locus as a variant in the comp track, those that have the same alternate allele)'),
    Metric('concordantRate',  'conc rate',   'Concordance rate'),
    Metric('variantRate',     'var/loci',    'Variants per loci rate'),
    Metric('basesPerVariant', 'bp/var',      'Bases per variant rate'),
    Metric('hetHomRatio',     'het/hom',     'Heterozygosity to homozygosity ratio'),
    Metric('tiTvRatio',       'ti/tv',       'Transition to transversion ratio'),
])


def get_parse_qc_sample_report(cnf):
    def _parse_qc_sample_report(report_fpath):
        """ returns row_per_sample =
                dict(metricName=None, value=None,
                isMain=True, quality='More is better')
        """

        records = OrderedDefaultDict(Record)
        rest_headers = []
        # metrics[metric_name]['meta']
        with open(report_fpath) as f:
            # parsing Sample name and Database columns
            main_value_col_id = None
            novelty_col_id = None
            for line in f:
                if not line.strip():
                    continue

                elif line.startswith(metrics_header):
                    if cnf.quality_control.db_for_summary in line:
                        main_value_col_id = line.split().index(cnf.quality_control.db_for_summary)

                    if novelty_header in line:
                        novelty_col_id = line.split().index(novelty_header)

                    rest_headers = line.split()[2:]

                elif novelty_col_id:
                    # parsing rest of the report
                    metric_name = line.split()[0]
                    novelty = line.split()[novelty_col_id]

                    records[metric_name].metric = METRICS[metric_name]
                    records[metric_name].meta[novelty] = dict(zip(rest_headers, line.split()[2:]))

                    if novelty == main_novelty:
                        val = line.split()[main_value_col_id].replace(' ', '').replace(',', '')

                        num_chars = []
                        unit_chars = []

                        i = 0
                        while i < len(val) and (val[i].isdigit() or val[i] == '.'):
                            num_chars += val[i]
                            i += 1
                        while i < len(val):
                            unit_chars += val[i]
                            i += 1

                        val_num = ''.join(num_chars)
                        val_unit = ''.join(unit_chars)

                        if val_unit:
                            records[metric_name].metric.unit = val_unit

                        try:
                            val = int(val_num)
                        except ValueError:
                            try:
                                val = float(val_num)
                            except ValueError:
                                val = val_num

                        records[metric_name].value = val

        return records

    return _parse_qc_sample_report

