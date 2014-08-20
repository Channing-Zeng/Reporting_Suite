import json
from os.path import join
from source.bcbio_structure import BCBioStructure
from source.calling_process import call_subprocess

from source.logger import step_greetings
from source.reporting import Metric, Record, save_json
from source.tools_from_cnf import get_gatk_cmdline
from source.utils import OrderedDefaultDict


final_report_ext = '.txt'


gatk_metrics = Metric.to_dict([
    Metric('nEvalVariants',   'Total',       'Total variants evaluated'),
    Metric('nSNPs',           'SNP',         'SNPs'),
    Metric('nInsertions',     'Ins',         'Insertions'),
    Metric('nDeletions',      'Del',         'Deletions'),
    Metric('nVariantsAtComp', 'At comp',     'Number of eval sites at comp sites (that is, sharing the same locus as a variant in the comp track, regardless of whether the alternate allele is the same)'),
    Metric('compRate',        'Comp rate',   'Percentage of eval sites at comp sites'),
    Metric('nConcordant',     'Concord',     'Number of concordant sites (that is, for the sites that share the same locus as a variant in the comp track, those that have the same alternate allele)'),
    Metric('concordantRate',  'Conc rate',   'Concordance rate'),
    Metric('variantRate',     'Var/loci',    'Variants per loci rate'),
    Metric('basesPerVariant', 'Bp/var',      'Bases per variant rate'),
    Metric('hetHomRatio',     'Het/hom',     'Heterozygosity to homozygosity ratio'),
    Metric('tiTvRatio',       'Ti/Tv',       'Transition to transversion ratio'),
])


def _dict_to_objects(report_dict, cnf_databases, cnf_novelty, main_db):
    records = []

    for m_name, cur_metric in gatk_metrics.items():
        rec = Record()
        records.append(rec)
        rec.metric = cur_metric

        for cur_novelty in cnf_novelty:
            rec.meta[cur_novelty] = dict()

            sum_for_all_novelties = 0.0
            for cur_database in cnf_databases:
                value = report_dict[m_name][cur_database][cur_novelty]
                try:
                    value = int(value)
                except ValueError:
                    value = float(value)
                rec.meta[cur_novelty][cur_database] = value
                sum_for_all_novelties += value

            average = sum_for_all_novelties / len(cnf_databases)
            rec.meta[cur_novelty]['average'] = average

        rec.value = rec.meta['all'][main_db]

    return records


def gatk_qc(cnf, vcf_fpath):
    step_greetings('Quality control reports')

    qc_cnf = cnf.quality_control
    cnf_databases = qc_cnf.database_vcfs
    cnf_novelty = qc_cnf.novelty
    cnf_metrics = qc_cnf.metrics

    gatk = get_gatk_cmdline(cnf)

    ref_fpath = cnf.genome['seq']
    report_fpath = join(cnf.work_dir, cnf.name + '_gatk.report')

    cmdline = ('{gatk} -R {ref_fpath} -T VariantEval '
               '--eval:tmp {vcf_fpath} -o {report_fpath}'
    ).format(**locals())

    if 'dbsnp' in cnf_databases:
        cmdline += ' -D ' + cnf_databases['dbsnp']
    for db_name, db_path in cnf_databases.items():
        if not db_name == 'dbsnp':
            cmdline += ' -comp:' + db_name + ' ' + db_path

    call_subprocess(cnf, cmdline, None, report_fpath, stdout_to_outputfile=False,
         to_remove=[vcf_fpath + '.idx'])

    report_dict = _parse_gatk_report(report_fpath, cnf_databases.keys(), cnf_novelty)
    records = _dict_to_objects(report_dict, cnf_databases.keys(), cnf_novelty, cnf.quality_control.db_for_summary)

    f_basename = join(cnf.output_dir, cnf.name + '-' + cnf.caller + '.' + BCBioStructure.varqc_name)
    save_json(records, f_basename + '.json')
    final_report_fpath = f_basename + final_report_ext

    _make_final_report(records, final_report_fpath, cnf_databases.keys(), cnf_novelty)
    return final_report_fpath, records


def _parse_gatk_report(report_filename, cnf_databases, cnf_novelty):
    database_col_name = 'CompRod'
    database_col_id = None
    novelty_col_name = 'Novelty'
    novelty_col_id = None

    report = dict()
    comments_section = False
    cur_header = []
    cur_metrics_ids = []
    for line in open(report_filename):
        if not line.strip():
            continue
        if line.startswith('#'):  # comment line
            comments_section = True
            continue
        elif comments_section:
            comments_section = False
            line = line.replace('variantRatePerBp', 'basesPerVariant')
            cur_header = line.split()
            cur_metrics_ids = []
            database_col_id = cur_header.index(database_col_name)
            novelty_col_id = cur_header.index(novelty_col_name)
            for metric in gatk_metrics.keys():
                if metric in cur_header:
                    cur_metrics_ids.append(cur_header.index(metric))
                    if metric not in report:
                        report[metric] = dict()
        elif cur_metrics_ids:  # process lines only if there are metrics in current section
            values = line.split()
            if len(values) <= max(cur_metrics_ids + [novelty_col_id, database_col_id]):
                continue

            cur_database = values[database_col_id]
            cur_novelty = values[novelty_col_id]
            if (cur_database not in cnf_databases) or (cur_novelty not in cnf_novelty):
                continue
            for metric_id in cur_metrics_ids:
                if cur_database not in report[cur_header[metric_id]]:
                    report[cur_header[metric_id]][cur_database] = dict()
                report[cur_header[metric_id]][cur_database][cur_novelty] = values[metric_id]
    return report


def _make_final_report(records, report_fpath, cnf_databases, cnf_novelties):
    header = ['Metric', 'Novelty'] + cnf_databases + ['Average']
    rows = [header]

    for rec in records:
        for novelty in cnf_novelties:
            row = [rec.metric.name, novelty]
            for db in cnf_databases + ['average']:
                row.append(rec.metric.format(rec.meta[novelty][db]))

            rows.append(row)

    col_widths = [len(h) for h in header]
    for row in rows:
        for i, value in enumerate(row):
            col_widths[i] = max(len(value), col_widths[i])

    with open(report_fpath, 'w') as out:
        for row in rows:
            out.write('  '.join('%-*s' % (col_width, value)
                      for col_width, value in zip(col_widths, row)) + "\r\n")

