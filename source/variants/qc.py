import gzip
import re
from os.path import join

from source.variants import vcf_parser

import source
from source.logger import step_greetings, warn
from source.file_utils import open_gzipsafe
from source.reporting.reporting import Metric, MetricStorage, ReportSection, SampleReport
from source.utils import get_db_path
from source.variants.vcf_processing import get_sample_column_index
import source.variants.vcf_processing as vcf_processing


metric_storage = MetricStorage(
    sections=[
        ReportSection('basic', '', [
            Metric('Total variants',      'Total',               'Total number of passed variants with'),
            Metric('SNPs',                'SNP',                 'SNPs'),
            Metric('Insertions',          'Ins',                 'Insertions'),
            Metric('Deletions',           'Del',                 'Deletions'),
            Metric('Novel',               'Novel',               'Novel (not in dbSNP or Cosmic'),
            Metric('Novel, %',            '%',                   '% novel varinats', unit='%'),
            # Metric('dbsnp_loci',          'Loci in dnSNP',       'Loci in dbSNP (just CHROM:POS matches, regardless if allele is the same)'),
            # Metric('dbsnp_loci_percent',  '%',                   '% loci in dbSNP (just CHROM:POS matches, regardless if allele is the same)', unit='%'),
            Metric('In dbSNP',            'dbSNP',               'Variants in dbSNP'),
            Metric('In dbSNP, %',         '%',                   '% variants in dbSNP', unit='%'),
            # Metric('cosmic_loci',         'Loci in Cosmic',      'Loci in Cosmic (just CHROM:POS matches, regardless if allele is the same)'),
            # Metric('cosmic_loci_percent', '%',                   '% loci in Cosmic (just CHROM:POS matches, regardless if allele is the same)', unit='%'),
            Metric('In Cosmic',           'Cosmic',              'Variants in Cosmic'),
            Metric('In Cosmic, %',        '%',                   '% variants in Cosmic', unit='%'),
            # Metric('bases_per_variant',   'Bp/var',              'Reference bases per variant', quality='Less is better'),
            Metric('Het/hom',             'Het/hom',             'Heterozygosity to homozygosity ratio'),
            Metric('Ti/tv',               'Ti/tv',               'Transition (T<->C, A<->G) to transversion (A<->C, C<->G, G<->T, T<->A) ratio. Should be 2 to 3 or higher (depending on the species and region)'),
            Metric('Total with rejected', 'Total with rejected', 'Total number of records in VCF, regardless FILTER column'),
        ])
    ]
)


def set_db_versions(cnf):
    m = metric_storage.find_metric('In dbSNP')
    m.description = 'Variants in dbSNP' + get_db_version(cnf, 'dbsnp')
    m = metric_storage.find_metric('In Cosmic')
    m.description = 'Variants in Cosmic' + get_db_version(cnf, 'cosmic')


def get_db_version(cnf, dbname):
    db_fpath = get_db_path(cnf, cnf.annotation[dbname], dbname)
    if db_fpath:
        if db_fpath.endswith('.gz'):
            opener = gzip.open
        else:
            opener = open
        with opener(db_fpath, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    break
                if (dbname == 'dbsnp' and line.startswith('##dbSNP_BUILD_ID=')) or \
                        (dbname == 'cosmic' and line.startswith('##source=')):
                    db_version = re.findall('\d+', line.split('=')[1])[0]
                    return ' v' + str(db_version)
    return ''


def make_report(cnf, vcf_fpath, sample):
    set_db_versions(cnf)
    step_greetings('Quality control reports')

    total_with_rejected = 0
    total = 0
    snps = 0
    inss = 0
    dels = 0
    dbsnps = 0
    cosmics = 0
    novels = 0
    hets = 0
    homs = 0
    transitions = 0
    transversions = 0

    with open_gzipsafe(vcf_fpath) as f:
        reader = vcf_parser.Reader(f)
        for rec in (vcf_processing.Record(rec, vcf_fpath, i) for i, rec in enumerate(reader)):
            total_with_rejected += 1

            if not rec.FILTER or rec.FILTER == 'PASS':
                if rec.FILTER:
                    warn('Warn: ' + rec.get_variant() + ' FILTER=' + str(rec.FILTER))

                total += 1

                if rec.is_snp:
                    snps += 1
                    if rec.is_transition:
                        transitions += 1
                    elif len(rec.ALT) == 1:
                        transversions += 1
                elif rec.is_indel:
                    if rec.is_deletion:
                        dels += 1
                    elif len(rec.ALT) == 1:
                        inss += 1

                if not rec.ID:
                    novels += 1
                else:
                    ids = rec.ID
                    if isinstance(ids, basestring):
                        ids = [ids]
                    if any(id.startswith('COS') for id in ids):
                        cosmics += 1
                    if any(id.startswith('rs') for id in ids):
                        dbsnps += 1

                call = rec.samples[0]
                if call.called:
                    if call.gt_type == 1:
                        hets += 1
                    elif call.gt_type == 2:
                        homs += 1

    report = SampleReport(sample, metric_storage=metric_storage)
    report.add_record('Total variants',      total)
    report.add_record('SNPs',                snps)
    report.add_record('Insertions',          inss)
    report.add_record('Deletions',           dels)
    report.add_record('Novel',               novels)
    report.add_record('Novel, %',            1.0 * novels / total if total else None)
    report.add_record('In dbSNP',            dbsnps)
    report.add_record('In dbSNP, %',         1.0 * dbsnps / total if total else None)
    report.add_record('In Cosmic',           cosmics)
    report.add_record('In Cosmic, %',        1.0 * cosmics / total if total else None)
    report.add_record('Het/hom',             float(hets) / homs if homs != 0 else None)
    report.add_record('Ti/tv',               float(transitions) / transversions if transversions != 0 else None)
    report.add_record('Total with rejected', total_with_rejected)

    return report


def save_report(cnf, report, sample, callername, output_dir, proc_name=source.varqc_name):
    f_basepath = join(output_dir, sample.name + ('-' + callername if callername else '') + '.' + proc_name)
    report.save_txt(join(output_dir, f_basepath + '.txt'))
    report.save_json(join(output_dir, f_basepath + '.json'))
    report.save_html(cnf,
        join(output_dir, sample.name + (('-' + callername) if callername else '') + '.' + proc_name + '.html'),
        caption='Variant QC for ' + sample.name + ((' (caller: ' + callername + ')') if callername else ''))
    return report
