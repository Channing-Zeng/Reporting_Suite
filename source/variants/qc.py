from os.path import join
from source.calling_process import call_subprocess

from source.logger import step_greetings, warn
from source.file_utils import open_gzipsafe
from source.reporting import Metric, Record, MetricStorage, ReportSection, SampleReport
from ext_modules import vcf_parser
from source.variants.vcf_processing import get_sample_column_index, iterate_vcf
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
            Metric('In dbSNP',            'dnSNP',               'Variants in dbSNP'),
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


def make_report(cnf, vcf_fpath, sample):
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

    main_sample_index = get_sample_column_index(vcf_fpath, sample.name)

    with open_gzipsafe(vcf_fpath) as f:
        reader = vcf_parser.Reader(f)
        for rec in (vcf_processing.Record(rec, vcf_fpath, i) for i, rec in enumerate(reader)):
            total_with_rejected += 1

            if rec.FILTER == [] or rec.FILTER == 'PASS':
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

                call = rec.get_main_sample(main_sample_index)
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

    save_report(cnf, report)
    return report


def save_report(cnf, report):
    f_basepath = join(cnf.output_dir, cnf.name + ('-' + cnf.caller if cnf.caller else '') + '.' + cnf.proc_name)

    report.save_txt(cnf.output_dir, f_basepath)
    report.save_json(cnf.output_dir, f_basepath)


