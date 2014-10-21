#!/usr/bin/env python
import math
from collections import defaultdict, OrderedDict
import os
from os.path import join
import sys
from ext_modules.simplejson import load
from source.bcbio_structure import BCBioStructure
from source.calling_process import call_subprocess, call_pipe
from source.file_utils import verify_file
from source.logger import info, err, step_greetings, critical, send_email
from source.reporting import write_tsv_rows, Record, SampleReport
from source.targetcov.cov import make_and_save_general_report, make_targetseq_reports
from source.tools_from_cnf import get_script_cmdline

from source.utils import OrderedDefaultDict, get_chr_len_fpath
from source.utils import median, mean


def cnv_reports(cnf, bcbio_structure):
    step_greetings('Coverage statistics for each gene for all samples')

    info('Collecting sample reports...')
    summary_report_fpath_by_sample = bcbio_structure.get_targetcov_report_fpaths_by_sample('json')
    gene_report_fpaths_by_sample = bcbio_structure.get_gene_report_fpaths_by_sample()

    for sample in bcbio_structure.samples:
        if (not verify_file(summary_report_fpath_by_sample.get(sample.name)) or
            not verify_file(gene_report_fpaths_by_sample.get(sample.name))):
            # TargetSeq was not run but needed, thus running.
            proc_name, name, output_dir = cnf.proc_name, cnf.name, cnf.output_dir
            cnf.proc_name, cnf.name, cnf.output_dir = \
                BCBioStructure.targetseq_name, sample.name, join(sample.dirpath, BCBioStructure.targetseq_name)
            make_targetseq_reports(cnf, sample)
            cnf.proc_name, cnf.name, cnf.output_dir = proc_name, name, output_dir

    info('Calculating normalized coverages for CNV...')
    amplicon_cnv_rows, gene_cnv_rows = _summarize_copy_number(cnf, bcbio_structure,
        gene_report_fpaths_by_sample, summary_report_fpath_by_sample)

    cnv_ampl_report_fpath, cnv_gene_ampl_report_fpath = None, None
    if amplicon_cnv_rows:
        cnv_ampl_report_fpath = write_tsv_rows(amplicon_cnv_rows, cnf.output_dir, BCBioStructure.seq2c_name + '_amplicons')
    if gene_cnv_rows:
        cnv_gene_ampl_report_fpath = write_tsv_rows(gene_cnv_rows, cnf.output_dir, BCBioStructure.seq2c_name)

    info()
    info('*' * 70)
    msg = 'Seq2C was not generated.'
    if cnv_ampl_report_fpath or cnv_gene_ampl_report_fpath:
        info('Seq2C:')
        msg = 'Seq2C:'
        if cnv_ampl_report_fpath:
            info('  Amplicon level:      ' + cnv_ampl_report_fpath)
            msg += '\n  Amplicon level:      ' + cnv_ampl_report_fpath
        if cnv_gene_ampl_report_fpath:
            info('  Gene-Amplicon level: ' + cnv_gene_ampl_report_fpath)
            msg += '\n  Gene-Amplicon level: ' + cnv_gene_ampl_report_fpath

    # if msg:
        # send_email(msg)

    return [cnv_ampl_report_fpath, cnv_gene_ampl_report_fpath]


def _get_lines_by_region_type(report_fpath, region_type):
    gene_summary_lines = []

    with open(report_fpath, 'r') as f:
        for line in f:
            if region_type in line:
                ts = line.split()
                gene_summary_lines.append(ts[:10])

    if not gene_summary_lines:
        critical('Regions of type ' + region_type + ' not found in ' + report_fpath)

    return gene_summary_lines


def _summarize_copy_number(cnf, bcbio_structure, gene_reports_by_sample, report_fpath_by_sample):
    amplicon_summary_lines = []
    gene_amplicon_summary_lines = []
    mapped_reads_by_sample = OrderedDict()

    for sample_name, gene_report_fpath in gene_reports_by_sample.items():
        json_fpath = report_fpath_by_sample[sample_name]

        # amplicon_summary_lines += _get_lines_by_region_type(gene_report_fpath, 'Amplicon')
        gene_amplicon_summary_lines += _get_lines_by_region_type(gene_report_fpath, 'Gene-Amplicon')

        with open(json_fpath) as f:
            data = load(f, object_pairs_hook=OrderedDict)
        sample = next((sample for sample in bcbio_structure.samples
                       if sample.name == sample_name), None)
        if not sample:
            continue
        cov_report = SampleReport.load(data, sample, bcbio_structure)

        mapped_reads_by_sample[sample_name] = int(next(
            rec.value for rec in cov_report.records
            if rec.metric.name == 'Mapped reads'))

    # results = run_copy_number(mapped_reads_by_sample, gene_summary_lines)

    results_amplicon = None
    # results_amplicon = run_copy_number__cov2cnv2(cnf, mapped_reads_by_sample, amplicon_summary_lines)
    results_gene_amplicon = run_copy_number__cov2cnv2(cnf, mapped_reads_by_sample, gene_amplicon_summary_lines)

    # save_results_separate_for_samples(results)

    return [results_amplicon, results_gene_amplicon]


def run_copy_number__cov2cnv2(cnf, mapped_reads_by_sample, gene_summary_lines):
    """
    Normalize the coverage from targeted sequencing to CNV log2 ratio. The algorithm assumes the medium
    is diploid, thus not suitable for homogeneous samples (e.g. parent-child).
    """

    mapped_read_fpath = join(cnf.work_dir, 'mapped_reads_by_sample.txt')
    with open(mapped_read_fpath, 'w') as f:
        for sample_name, mapped_reads in mapped_reads_by_sample.items():
            f.write(sample_name + '\t' + str(mapped_reads) + '\n')

    gene_depths_fpaths = join(cnf.work_dir, 'gene_depths.txt')
    with open(gene_depths_fpaths, 'w') as f:
        for tokens in gene_summary_lines:
            sample, chrom, s, e, gene, exon, strand, tag, size, cov = tokens
            s, e, size, cov = [''.join(c for c in l if c.isdigit()) for l in [s, e, size, cov]]
            reordered = sample, gene, chrom, s, e, tag, size, cov
            f.write('\t'.join(reordered) + '\n')

    os.environ['PERL5LIB'] = '/group/cancer_informatics/tools_resources/NGS/lib/perl5' + \
        (':' + os.environ['PERL5LIB'] if 'PERL5LIB' in os.environ else '')
    cov2cnv2 = get_script_cmdline(cnf, 'perl', 'cov2cnv2')
    if not cov2cnv2: sys.exit(1)
    cmdline = '{cov2cnv2} {mapped_read_fpath} {gene_depths_fpaths}'.format(**locals())

    proc = call_pipe(cnf, cmdline)
    results = [l.split() for l in proc.stdout]
    return results


# def save_results_separate_for_samples(results):
#     header = results[0]
#
#     results_per_sample = OrderedDict()
#
#     for fields in results[1:]:
#         sample_name = fields[0]
#         if sample_name not in results_per_sample:
#             results_per_sample[sample_name] = [header]
#
#         results_per_sample[sample_name].append(fields)
#
#     for sample_name, fields in results_per_sample.items():


def __proc_sample(sample, norm_depths_by_sample, factors_by_gene, med_depth, median_depth_by_sample):
    info('  Processing ' + sample + '...')
    norm_depths_by_gene = OrderedDefaultDict(dict)
    norm2 = OrderedDefaultDict(dict)
    norm3 = OrderedDefaultDict(dict)

    i = 0
    for gene, norm_depth_by_sample in norm_depths_by_sample.items():
        if i % 1000 == 0:
            info('    ' + sample + ': processing gene #{0:,}: {1}'.format(i, str(gene)))
        i += 1

        if sample not in norm_depth_by_sample:
            continue

        gene_norm_depth = norm_depth_by_sample[sample] * factors_by_gene[gene] + 0.1  # norm1b

        norm_depths_by_gene[gene][sample] = gene_norm_depth

        norm2[gene][sample] = math.log(gene_norm_depth / med_depth, 2) if med_depth else 0

        norm3[gene][sample] = math.log(gene_norm_depth / median_depth_by_sample[sample], 2) if \
            median_depth_by_sample[sample] else 0

    info('  ' + sample + ': Done. Processed {0:,} genes.'.format(i))

    return norm_depths_by_gene, norm2, norm3


def run_copy_number(mapped_reads_by_sample, gene_depth):
    mapped_reads_by_sample = {k: v for k, v in mapped_reads_by_sample.items() if 'Undetermined' not in k}

    info('Parsing rows...')
    records = _report_row_to_objects(gene_depth)

    info('Calculating depths normalized by samples...')
    norm_depths_by_sample = _get_norm_depths_by_sample(mapped_reads_by_sample, records)

    med_depth = median([depth for gn, vs in norm_depths_by_sample.items() for sn, depth in vs.items()])

    info('Getting factors by genes...')
    factors_by_gene = _get_factors_by_gene(norm_depths_by_sample, med_depth)

    info('Getting median norm...')
    median_depth_by_sample = _get_med_norm_depths(mapped_reads_by_sample, norm_depths_by_sample)

    info()
    info('Norm depth by gene, norm2, norm3...')

    samples = set(rec.sample_name for rec in records)

    norm_depths_by_gene = OrderedDefaultDict(dict)
    norm2 = OrderedDefaultDict(dict)
    norm3 = OrderedDefaultDict(dict)
    for sample in samples:
        info('  Processing ' + sample + '...')
        i = 0
        for gene, norm_depth_by_sample in norm_depths_by_sample.items():
            if i % 1000 == 0:
                info('    ' + sample + ': processing gene #{0:,}: {1}'.format(i, str(gene)))
            i += 1

            if sample not in norm_depth_by_sample:
                continue

            gene_norm_depth = norm_depth_by_sample[sample] * factors_by_gene[gene] + 0.1  # norm1b
            norm_depths_by_gene[gene][sample] = gene_norm_depth

            norm2[gene][sample] = math.log(gene_norm_depth / med_depth, 2) if med_depth else 0
            norm3[gene][sample] = math.log(gene_norm_depth / median_depth_by_sample[sample], 2) if median_depth_by_sample[sample] else 0

        info('  ' + sample + ': Done. Processed {0:,} genes.'.format(i))

    return _make_report_data(records, norm2, norm3, norm_depths_by_gene, norm_depths_by_sample)


def _get_med_norm_depths(mapped_reads_by_sample, norm_depths_by_sample):
    mean_depths_by_sample = defaultdict(list)

    for sample_name in mapped_reads_by_sample:
        for gn, v_by_sample in norm_depths_by_sample.items():
            if sample_name not in v_by_sample:
                info(sample_name + ' - no ' + gn)
            mean_depths_by_sample[sample_name].append(v_by_sample.get(sample_name) or 0)

    median_depth_by_sample = OrderedDict()
    for sample_name in mapped_reads_by_sample:
        median_depth_by_sample[sample_name] = median(mean_depths_by_sample[sample_name])

    return median_depth_by_sample


# mean_reads =  mean for all the samples mapped reads
# return factor_by_sample = sample mapped read/mean for all the samples mapped reads
def _get_factors_by_sample(mapped_reads_by_sample):
    factor_by_sample = dict()

    mean_reads = mean(mapped_reads_by_sample.values())
    for sample, mapped_reads in mapped_reads_by_sample.items():
        factor_by_sample[sample] = mean_reads / mapped_reads if mapped_reads else 0

    return factor_by_sample


#med_depth = median for all gene
#return factor_by_gene = median gene for all samples/ med_depth
def _get_factors_by_gene(norm_depths_1, med_depth):
    # mean_depth_by_genes = defaultdict(list)
    # for gene_name, values in norm_depths_by_gene:
    #     mean_depth_by_genes[gene_name].append(rec.mean_depth)

    factors_by_gene = dict()
    for gene_name, values in norm_depths_1.items():
        med = median(values.values())
        factors_by_gene[gene_name] = med_depth / med if med else 0

    return factors_by_gene


#MeanDepth_Norm1
# gene -> { sample -> [] }
def _get_norm_depths_by_sample(mapped_reads_by_sample, records):
    norm_depths = defaultdict(dict)

    factor_by_sample = _get_factors_by_sample(mapped_reads_by_sample)

    # Initialization
    for sample, factor in factor_by_sample.items():
        for rec in records:
            norm_depths[rec.gene_name][sample] = 0

    # Filling the dict with available balues
    for sample, factor in factor_by_sample.items():
        info('  ' + sample)

        for rec in [rec for rec in records if rec.sample_name == sample]:
            norm_depths[rec.gene_name][sample] = rec.mean_depth * factor

    return norm_depths


def _make_report_data(records, norm2, norm3, norm_depths_by_gene, norm_depths_by_seq_distr):
    header = ["Sample", "Gene", "Chr", "Start", "Stop", "Length", "MeanDepth", "MeanDepth_Norm1",
              "MeanDepth_Norm2", "log2Ratio_norm1", "log2Ratio_norm2"]
    report_data = []

    for rec in records:
        gene_name = rec.gene_name
        sample = rec.sample_name
        report_data.append(map(str,
           (sample, gene_name, rec.chrom, rec.start_position,
            rec.end_position, rec.size,
            '{0:.3f}'.format(rec.mean_depth, norm_depths_by_seq_distr[gene_name][sample]),
            '{0:.3f}'.format(norm_depths_by_gene[gene_name][sample]),
            '{0:.3f}'.format(norm2[gene_name][sample]),
            '{0:.3f}'.format(norm3[gene_name][sample]))))

    # report_data = sorted(report_data, key=itemgetter(2, 1))
    report_data = [header] + report_data

    return report_data


# sample_order = '''PC9_IRLR PC9_9291-L0B_1 PC9_IR-GM PC9_9291_6 PC9_9291-5 PC9_IR-4
# PC9_IR-6 PC9_9291-2 PC9 PC9_9291-3 gDNA'''.split()
#
# gene_order = 'PIK3CA KRAS NRAS HRAS BRAF EGFR'.split()

def _report_row_to_objects(gene_depth):
    inputs = OrderedDefaultDict(OrderedDict)
    # for gene_name in gene_order:
    #     for sample_name in sample_order:
    #         inputs[gene_name][sample_name] = None

    for read in gene_depth:
        rec = CovRec(*read)
        if 'Undetermined' not in rec.sample_name:
            inputs[rec.gene_name][rec.sample_name] = rec

    details_list = []
    for gene_dict in inputs.values():
        for rec in gene_dict.values():
            details_list.append(rec)

    return details_list


class CovRec:
    def __init__(self, sample_name=None, chrom=None, start_position=None, end_position=None, gene_name=None,
                 type="Gene-Amplicon", size=None, mean_depth=None):
        self.sample_name = sample_name
        self.chrom = chrom
        self.start_position = int(start_position)
        self.end_position = int(end_position)
        self.gene_name = gene_name
        self.type = type
        self.size = int(size)
        self.mean_depth = float(mean_depth)

    def __str__(self):
        values = [self.sample_name, self.chrom, self.start_position, self.end_position, self.gene_name, self.type, self.size,
                  self.mean_depth]
        return '"' + '\t'.join(map(str, values)) + '"'

    def __repr__(self):
        return repr(
            (self.gene_name, self.sample_name))

