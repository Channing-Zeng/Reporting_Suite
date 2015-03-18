#!/usr/bin/env python
from genericpath import isfile
import math
from collections import defaultdict, OrderedDict
import os
from os.path import join, splitext, basename
import sys
from ext_modules.simplejson import load
from source.bcbio_structure import BCBioStructure
from source.calling_process import call_subprocess, call_pipe, call
from source.file_utils import verify_file, adjust_path
from source.logger import info, err, step_greetings, critical, send_email, warn
from source.reporting import write_tsv_rows, Record, SampleReport
from source.targetcov.cov import make_and_save_general_report, make_targetseq_reports
from source.tools_from_cnf import get_script_cmdline, get_system_path
from source.utils import OrderedDefaultDict, get_chr_len_fpath
from source.utils import median, mean
import source


def _get_exons_and_genes(cnf, sample):
    exons_bed_fpath = adjust_path(cnf.exons) if cnf.exons else adjust_path(cnf.genomes[sample.genome].exons)
    info('Exons: ' + exons_bed_fpath)

    genes_fpath = None
    if cnf.genes:
        genes_fpath = adjust_path(cnf.genes)
        info('Custom genes list: ' + genes_fpath)

    return exons_bed_fpath, genes_fpath


def cnv_reports(cnf, bcbio_structure):
    step_greetings('Coverage statistics for each gene for all samples')

    for sample in bcbio_structure.samples:
        if not verify_file(sample.targetcov_json_fpath) or not verify_file(sample.targetcov_detailed_tsv):
            exons_bed_fpath, genes_fpath = _get_exons_and_genes(cnf, sample)

            # TargetSeq was not run but needed, thus running.
            proc_name, name, output_dir = cnf.proc_name, cnf.name, cnf.output_dir
            cnf.proc_name, cnf.name, cnf.output_dir = source.targetseq_name, sample.name, join(sample.dirpath, source.targetseq_name)

            make_targetseq_reports(cnf, sample, exons_bed_fpath, genes_fpath)  # cnf.vcfs_by_callername

            cnf.proc_name, cnf.name, cnf.output_dir = proc_name, name, output_dir

    info('Calculating normalized coverages for CNV...')
    cnv_report_fpath, cnv_report_fpath__mine = _seq2c(cnf, bcbio_structure)

    info()
    info('*' * 70)
    msg = 'Seq2C was not generated.'
    if cnv_report_fpath or cnv_report_fpath__mine:
        info('Seq2C:')
        msg = 'Seq2C:'
        if cnv_report_fpath:
            info('   ' + cnv_report_fpath)
            msg += '\n   ' + cnv_report_fpath
        if cnv_report_fpath__mine:
            info('   Mine: ' + cnv_report_fpath__mine)
            msg += '\n   Mine: ' + cnv_report_fpath__mine

    # if msg:
        # send_email(msg)

    return [cnv_report_fpath, cnv_report_fpath__mine]


def _get_whole_genes_and_amlicons(report_fpath):
    gene_summary_lines = []

    info('Reading from ' + report_fpath)

    with open(report_fpath, 'r') as f:
        for i, line in enumerate(f):
            if 'Capture' in line or 'Whole-Gene' in line:
                ts = line.split('\t')
                # #Sample  Chr  Start  End  Size  Gene  Strand  Feature  Biotype  Min depth  Ave depth  Std dev.  W/n 20% of ave
                ts = ts[:4] +               ts[5:6] + ts[7:8] + ts[4:5] + ts[10:11]
                # sample_name, chrom, s, e, gene,     tag,      size,     cov
                gene_summary_lines.append(ts)

    if not gene_summary_lines:
        critical('No Capture or Whole-Gene is not found in ' + report_fpath)

    return gene_summary_lines


def _seq2c(cnf, bcbio_structure):
    """
    Normalize the coverage from targeted sequencing to CNV log2 ratio. The algorithm assumes the medium
    is diploid, thus not suitable for homogeneous samples (e.g. parent-child).
    """

    info('Getting reads and cov stats with seq2cov.pl and bam2readsl.pl')
    read_stats_fpath, combined_gene_depths_fpath = __get_mapped_reads_and_cov_by_seq2c_itself(cnf, bcbio_structure.samples)
    info()

    if not read_stats_fpath or not combined_gene_depths_fpath:
        err('Error: no read_stats_fpath or combined_gene_depths_fpath by Seq2C, making ours...')
        return None, None

    info('Getting old way reads and cov stats, but with amplicons')
    read_stats_fpath__mine, combined_gene_depths_fpath__mine = __get_mapped_reads_and_cov(cnf.work_dir, bcbio_structure)
    info()

    cnv_gene_ampl_report_fpath = join(cnf.output_dir, BCBioStructure.seq2c_name + '.tsv')
    cnv_gene_ampl_report_fpath__mine = join(cnf.work_dir, BCBioStructure.seq2c_name + '.mine.tsv')

    # cnv_gene_ampl_report_fpath = __cov2cnv2(cnf, read_stats_fpath, combined_gene_depths_fpath, cnv_gene_ampl_report_fpath)
    cnv_gene_ampl_report_fpath = __new_seq2c(cnf, read_stats_fpath, combined_gene_depths_fpath, cnv_gene_ampl_report_fpath)
    cnv_gene_ampl_report_fpath__mine = __new_seq2c(cnf, read_stats_fpath__mine, combined_gene_depths_fpath__mine, cnv_gene_ampl_report_fpath__mine)

    # run_copy_number__cov2cnv2(cnf, mapped_reads_by_sample, amplicon_summary_lines, cnv_ampl_report_fpath)

    # save_results_separate_for_samples(results)

    return [cnv_gene_ampl_report_fpath, cnv_gene_ampl_report_fpath__mine]


def __new_seq2c(cnf, read_stats_fpath, combined_gene_depths_fpath, output_fpath):
    cov2lr = get_script_cmdline(cnf, 'perl', join('Seq2C', 'cov2lr.pl'))
    if not cov2lr: sys.exit(1)
    cov2lr_output = join(cnf.work_dir, splitext(basename(output_fpath))[0] + '.cov2lr.tsv')
    cmdline = '{cov2lr} -a {read_stats_fpath} {combined_gene_depths_fpath}'.format(**locals())
    call(cnf, cmdline, cov2lr_output, exit_on_error=False)
    info()

    if not verify_file(cov2lr_output):
        return None

    lr2gene = get_script_cmdline(cnf, 'perl', join('Seq2C', 'lr2gene.pl'))
    if not lr2gene: sys.exit(1)
    cmdline = lr2gene + ' ' + cov2lr_output
    res = call(cnf, cmdline, output_fpath, exit_on_error=False)
    info()

    if not verify_file(output_fpath):
        return None

    return res

# def __cov2cnv2(cnf, mapped_reads_by_sample, cov_info, output_fpath):
#     mapped_read_fpath, gene_depths_fpath = __get_mapped_reads_and_cov(
#         cnf.work_dir, bcbio_structure, report_fpath_by_sample, gene_reports_by_sample)
#
#     cov2cnv2 = get_script_cmdline(cnf, 'perl', 'cov2cnv2')
#     if not cov2cnv2: sys.exit(1)
#     cmdline = '{cov2cnv2} {mapped_read_fpath} {gene_depths_fpath}'.format(**locals())
#     if call(cnf, cmdline, output_fpath):
#         return output_fpath


def __get_mapped_reads_and_cov_by_seq2c_itself(cnf, samples):
    # GENE DEPTHS
    info()
    info('Combining gene depths...')

    combined_gene_depths_fpath = join(cnf.work_dir, 'gene_depths.seq2c.txt')
    with open(combined_gene_depths_fpath, 'w') as out:
        for i, sample in enumerate(samples):
            seq2c_output = join(
                sample.dirpath,
                BCBioStructure.targetseq_dir,
                sample.name + '.' + \
                BCBioStructure.targetseq_name + '_' + \
                BCBioStructure.seq2c_seq2cov_ending)
            if not isfile(seq2c_output):
                return None, None
            with open(seq2c_output) as inp:
                for l in inp:
                    if l.startswith('Sample\tGene\tChr\t'):
                        if i == 0:
                            out.write(l)
                    else:
                        out.write(l)

    info('Saved to ' + combined_gene_depths_fpath)
    info()

    # READ STATS
    info('Getting read counts...')
    bam2reads = get_script_cmdline(cnf, 'perl', join('Seq2C', 'bam2reads.pl'))
    if not bam2reads: sys.exit(1)

    bam2reads_list_of_bams_fpath = join(cnf.work_dir, 'seq2c_list_of_bams.txt')
    with open(bam2reads_list_of_bams_fpath, 'w') as f:
        for sample in samples:
            f.write(sample.name + '\t' + sample.bam + '\n')

    samtools = get_system_path(cnf, 'samtools')
    read_stats_fpath = join(cnf.work_dir, 'seq2c_read_stats.txt')
    cmdline = '{bam2reads} {bam2reads_list_of_bams_fpath} -m {samtools}'.format(**locals())
    if not call(cnf, cmdline, read_stats_fpath): return None, None

    return read_stats_fpath, combined_gene_depths_fpath


def __get_mapped_reads_and_cov(work_dir, bcbio_structure):
    coverage_info = []
    mapped_reads_by_sample = OrderedDict()

    for sample in bcbio_structure.samples:
        for tokens in _get_whole_genes_and_amlicons(sample.targetcov_detailed_tsv):
            sample_name, chrom, s, e, size, gene, tag, cov = tokens
            s, e, size, cov = [''.join(c for c in l if c != ',') for l in [s, e, size, cov]]
            if cov != '.' and float(cov) != 0:
                reordered = sample_name, gene, chrom, s, e, tag, size, cov
                coverage_info.append(reordered)

        with open(sample.targetcov_json_fpath) as f:
            data = load(f, object_pairs_hook=OrderedDict)
        sample = next((s for s in bcbio_structure.samples if s.name == sample.name), None)
        if not sample: continue
        cov_report = SampleReport.load(data, sample, bcbio_structure)

        mapped_reads_by_sample[sample.name] = int(next(
            rec.value for rec in cov_report.records
            if rec.metric.name == 'Mapped reads'))

    # Writing results
    mapped_read_fpath = join(work_dir, 'mapped_reads_by_sample.txt')
    with open(mapped_read_fpath, 'w') as f:
        for sample_name, mapped_reads in mapped_reads_by_sample.items():
            f.write(sample_name + '\t' + str(mapped_reads) + '\n')

    gene_depths_fpath = join(work_dir, 'gene_depths.txt')
    with open(gene_depths_fpath, 'w') as f:
        for tokens in coverage_info:
            f.write('\t'.join(tokens) + '\n')

    return mapped_read_fpath, gene_depths_fpath



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


def my_copy_number(mapped_reads_by_sample, gene_depth):
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
                 type="Gene-Capture", size=None, mean_depth=None):
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

