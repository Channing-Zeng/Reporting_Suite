#!/usr/bin/env python
from genericpath import isfile
import math
from collections import defaultdict, OrderedDict
import os
from os.path import join, splitext, basename
from joblib import Parallel, delayed
import sys
from ext_modules.simplejson import load
from source.bcbio_structure import BCBioStructure
from source.calling_process import call_subprocess, call_pipe, call
from source.config import CallCnf
from source.file_utils import verify_file, adjust_path, iterate_file
from source.logger import info, err, step_greetings, critical, send_email, warn
from source.ngscat.bed_file import verify_bed
from source.reporting import write_tsv_rows, Record, SampleReport
from source.targetcov.bam_and_bed_utils import sort_bed, count_bed_cols, annotate_amplicons, cut, \
    group_and_merge_regions_by_gene
from source.targetcov.cov import make_and_save_general_report, make_targetseq_reports, remove_dups
from source.tools_from_cnf import get_script_cmdline, get_system_path
from source.utils import OrderedDefaultDict, get_chr_len_fpath
from source.utils import median, mean
import source


def cnv_reports(cnf, bcbio_structure):
    step_greetings('Coverage statistics for each gene for all samples')

    info('Calculating normalized coverages for CNV...')
    cnv_report_fpath, cnv_report_dups_fpath = _seq2c(cnf, bcbio_structure)

    info()
    info('*' * 70)
    msg = 'Seq2C was not generated.'
    if cnv_report_fpath or cnv_report_dups_fpath:
        info('Seq2C:')
        msg = 'Seq2C:'
        if cnv_report_fpath:
            info('   ' + cnv_report_fpath)
            msg += '\n   ' + cnv_report_fpath
        if cnv_report_dups_fpath:
            info('   With dups: ' + cnv_report_dups_fpath)
            msg += '\n   With dups: ' + cnv_report_dups_fpath

    # if msg:
        # send_email(msg)

    return [cnv_report_fpath, cnv_report_dups_fpath]


def _get_whole_genes_and_amlicons(report_fpath):
    gene_summary_lines = []

    info('Reading from ' + report_fpath)

    with open(report_fpath, 'r') as f:
        for i, line in enumerate(f):
            if 'Capture' in line or 'Whole-Gene' in line:  # TODO: this is incorect. We need to get Whole-Gene of Capture, not of CDSs!
                ts = line.split('\t')
                # Chr  Start  End  Size  Gene  Strand  Feature  Biotype  Min depth  Ave depth  Std dev.  W/n 20% of ave
                chrom, s, e, size, gene, strand, feature = ts[:7]
                ave_depth = ts[9]
                # chrom, s, e, gene,     feature,  size,     cov
                gene_summary_lines.append([chrom, s, e, size, gene, feature, ave_depth])
    if not gene_summary_lines:
        critical('No Capture or Whole-Gene is not found in ' + report_fpath)

    return gene_summary_lines


def seq2c_seq2cov(cnf, seq2cov, samtools, sample, bam_fpath, amplicons_bed, seq2c_output):
    def fn(l, i): return '\t'.join(l.split('\t')[:4])
    amplicons_bed = iterate_file(cnf, amplicons_bed, fn, suffix='4col')

    sample_name = sample.name

    cmdline = '{seq2cov} -m {samtools} -z -b {bam_fpath} -N {sample_name} {amplicons_bed}'.format(**locals())
    res = call(cnf, cmdline, seq2c_output)
    if not res:
        err('Could not run seq2cov.pl for ' + sample.name)
        return None


def _seq2c(cnf, bcbio_structure):
    """
    Normalize the coverage from targeted sequencing to CNV log2 ratio. The algorithm assumes the medium
    is diploid, thus not suitable for homogeneous samples (e.g. parent-child).
    """
    # Output fpaths:
    cnv_gene_ampl_report_fpath = join(cnf.output_dir, BCBioStructure.seq2c_name + '.tsv')
    cnv_gene_ampl_report_dups_fpath = join(cnf.output_dir, BCBioStructure.seq2c_name + '.dups.tsv')

    info('Getting reads and cov stats with seq2cov.pl and bam2readsl.pl')

    samtools = get_system_path(cnf, 'samtools')
    dedupped_bam_by_sample = dict(zip((s.name for s in bcbio_structure.samples), Parallel(n_jobs=cnf.threads) \
        (delayed(remove_dups)(CallCnf(cnf.__dict__), s.bam, samtools) for s in bcbio_structure.samples)))

    if samtools is None:
        samtools = get_system_path(cnf, 'samtools', is_critical=True)
        
    combined_gene_depths_fpath, combined_gene_depths_dups_fpath = __cov2cnv(cnf, bcbio_structure.samples, dedupped_bam_by_sample)
    mapped_reads_fpath, mapped_reads_dup_fpath = __get_mapped_reads(cnf, bcbio_structure, dedupped_bam_by_sample)
    info()
    if not mapped_reads_fpath or not combined_gene_depths_fpath:
        err('Error: no mapped_reads_fpath or combined_gene_depths_fpath by Seq2C')
        cnv_gene_ampl_report_fpath = None
    else:
        cnv_gene_ampl_report_fpath = __new_seq2c(cnf, mapped_reads_fpath, combined_gene_depths_fpath, cnv_gene_ampl_report_fpath)

    if not mapped_reads_dup_fpath or not combined_gene_depths_dups_fpath:
        warn('Warning: no mapped_reads_dup_fpath or combined_gene_depths_dups_fpath by Seq2C (with dups)')
        cnv_gene_ampl_report_dups_fpath = None
    else:
        cnv_gene_ampl_report_dups_fpath = __new_seq2c(cnf, mapped_reads_dup_fpath, combined_gene_depths_dups_fpath, cnv_gene_ampl_report_dups_fpath)

    # info('Getting old way reads and cov stats, but with amplicons')
    # info()
    # cnv_gene_ampl_report_fpath = __cov2cnv2(cnf, read_stats_fpath, combined_gene_depths_fpath, cnv_gene_ampl_report_fpath)
    # cnv_gene_ampl_report_fpath__mine = __new_seq2c(cnf, read_stats_fpath, combined_gene_depths_fpath__mine, cnv_gene_ampl_report_fpath__mine)

    # run_copy_number__cov2cnv2(cnf, mapped_reads_by_sample, amplicon_summary_lines, cnv_ampl_report_fpath)

    # save_results_separate_for_samples(results)

    return [cnv_gene_ampl_report_fpath, cnv_gene_ampl_report_dups_fpath]


def __new_seq2c(cnf, read_stats_fpath, combined_gene_depths_fpath, output_fpath):
    cov2lr = get_script_cmdline(cnf, 'perl', join('Seq2C', 'cov2lr.pl'), is_critical=True)
    cov2lr_output = join(cnf.work_dir, splitext(basename(output_fpath))[0] + '.cov2lr.tsv')
    cmdline = '{cov2lr} -a {read_stats_fpath} {combined_gene_depths_fpath}'.format(**locals())
    call(cnf, cmdline, cov2lr_output, exit_on_error=False)
    info()

    if not verify_file(cov2lr_output):
        return None

    lr2gene = get_script_cmdline(cnf, 'perl', join('Seq2C', 'lr2gene.pl'), is_critical=True)
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


def __prep_bed(cnf, bed_fpath, exons_bed):
    info()
    info('Preparing BED file for seq2c...')

    info()
    info('Sorting exons by (chrom, gene name, start); and merging regions within genes...')
    exons_bed = group_and_merge_regions_by_gene(cnf, exons_bed, keep_genes=True)

    info()
    info('bedtools-sotring BED file...')
    bed_fpath = sort_bed(cnf, bed_fpath)
    cols = count_bed_cols(bed_fpath)

    if cols < 4:
        info('Annotating amplicons with gene names from Ensembl...')
        bed_fpath = annotate_amplicons(cnf, bed_fpath, exons_bed)

    elif 8 > cols > 4:
        bed_fpath = cut(cnf, bed_fpath, 4)

    elif cols > 8:
        bed_fpath = cut(cnf, bed_fpath, 8)

    # removing regions with no gene annotation
    def f(l, i):
        if l.split('\t')[3].strip() == '.': return None
        else: return l
    bed_fpath = iterate_file(cnf, bed_fpath, f, 'filt')

    info('Done: ' + bed_fpath)
    return bed_fpath


def _run_cov2cnv(cnf, seq2cov, samtools, sample, bed_fpath, dedupped_bam_by_sample):
    if not verify_file(sample.seq2cov_output_fpath):
        seq2c_seq2cov(cnf, seq2cov, samtools, sample, sample.bam, bed_fpath, sample.seq2cov_output_fpath)

    if not verify_file(sample.seq2cov_output_dup_fpath):
        seq2c_seq2cov(cnf, seq2cov, samtools, sample, dedupped_bam_by_sample[sample.name], bed_fpath, sample.seq2cov_output_dup_fpath)


def __cov2cnv(cnf, samples, dedupped_bam_by_sample):
    info()
    info('Combining gene depths...')

    result = []
    bed_fpath = next((adjust_path(s.bed) for s in samples if s.bed), cnf.genome.az_exome)
    for s in samples:
        if not verify_file(s.seq2cov_output_fpath, silent=True):
            exons_bed_fpath = adjust_path(cnf.exons) if cnf.exons else adjust_path(cnf.genome.exons)
            verify_bed(bed_fpath, is_critical=True)
            bed_fpath = __prep_bed(cnf, bed_fpath, exons_bed_fpath)
            break

    info('Running first for the de-dupped version, then for the original version.')
    seq2cov = get_script_cmdline(cnf, 'perl', join('Seq2C', 'seq2cov.pl'), is_critical=True)
    samtools = get_system_path(cnf, 'samtools')
    Parallel(n_jobs=cnf.threads) \
        (delayed(_run_cov2cnv)(CallCnf(cnf.__dict__), seq2cov, samtools, s, bed_fpath, dedupped_bam_by_sample)
            for s in samples)

    for suf in '', '_dup':
        combined_gene_depths_fpath = join(cnf.work_dir, 'gene_depths.seq2c' + suf + '.txt')

        with open(combined_gene_depths_fpath, 'w') as out:
            for i, sample in enumerate(samples):
                seq2cov_fpath = sample.seq2cov_output_dup_fpath if suf else sample.seq2cov_output_fpath

                if not verify_file(seq2cov_fpath):
                    combined_gene_depths_fpath = None
                    break
                with open(seq2cov_fpath) as inp:
                    for l in inp:
                        if l.startswith('Sample\tGene\tChr\t'):
                            if i == 0:
                                out.write(l)
                        else:
                            out.write(l)

        result.append(combined_gene_depths_fpath)
        if combined_gene_depths_fpath:
            info('Saved to ' + combined_gene_depths_fpath)
            info()

    return result

    # # READ STATS
    # info('Getting read counts...')
    # bam2reads = get_script_cmdline(cnf, 'perl', join('Seq2C', 'bam2reads.pl'))
    # if not bam2reads: sys.exit(1)
    #
    # bam2reads_list_of_bams_fpath = join(cnf.work_dir, 'seq2c_list_of_bams.txt')
    # with open(bam2reads_list_of_bams_fpath, 'w') as f:
    #     for sample in samples:
    #         dedup_bam_fpath = remove_dups(cnf, sample.bam)
    #         f.write(sample.name + '\t' + dedup_bam_fpath + '\n')
    #
    # samtools = get_system_path(cnf, 'samtools')
    # read_stats_fpath = join(cnf.work_dir, 'seq2c_read_stats.txt')
    # cmdline = '{bam2reads} -m {samtools} {bam2reads_list_of_bams_fpath}'.format(**locals())
    # if not call(cnf, cmdline, read_stats_fpath): return None, None
    #
    # return read_stats_fpath, combined_gene_depths_fpath


def __get_mapped_reads(cnf, bcbio_structure, dedupped_bam_by_sample):
    # coverage_info = []
    mapped_reads_by_sample = OrderedDict()
    mapped_reads_by_sample_dup = OrderedDict()

    for sample in bcbio_structure.samples:
        # for tokens in _get_whole_genes_and_amlicons(sample.targetcov_detailed_tsv):
        #     chrom, s, e, size, gene, feature, ave_depth = tokens
        #     s, e, size, cov = [''.join(c for c in l if c != ',') for l in [s, e, size, ave_depth]]
        #     if cov != '.' and float(cov) != 0:
        #         coverage_info.append([sample.name, gene, chrom, s, e, feature, size, cov])

        if verify_file(sample.targetcov_json_fpath):
            with open(sample.targetcov_json_fpath) as f:
                data = load(f, object_pairs_hook=OrderedDict)
            sample = next((s for s in bcbio_structure.samples if s.name == sample.name), None)
            if not sample: continue
            cov_report = SampleReport.load(data, sample, bcbio_structure)

            mapped_reads = int(next(
                rec.value for rec in cov_report.records
                if rec.metric.name == 'Mapped reads'))

            try:
                dup_rate = float(next(
                    rec.value for rec in cov_report.records
                    if rec.metric.name == 'Duplication rate'))
            except StopIteration:
                dup_rate = float(next(
                    rec.value for rec in cov_report.records
                    if rec.metric.name == 'Duplication rate (picard)'))

            info(sample.name + ': ')
            info('  Mapped reads: ' + str(mapped_reads))
            info('  Dup rate: ' + str(dup_rate))
            mapped_reads_by_sample[sample.name] = int(float(mapped_reads) * (1 - dup_rate))
            mapped_reads_by_sample_dup[sample.name] = mapped_reads
            info('  Mapped not dup reads: ' + str(mapped_reads_by_sample[sample.name]))

    mapped_read_fpath = join(cnf.work_dir, 'mapped_reads_by_sample.txt')
    mapped_read_dup_fpath = join(cnf.work_dir, 'mapped_reads_by_sample_dups.txt')

    if len(mapped_reads_by_sample.keys()) == len(bcbio_structure.samples):
    # TargetSeq was run correctly; writing results
        with open(mapped_read_fpath, 'w') as f:
            for sample_name, mapped_reads in mapped_reads_by_sample.items():
                f.write(sample_name + '\t' + str(mapped_reads) + '\n')

        with open(mapped_read_dup_fpath, 'w') as f:
            for sample_name, mapped_reads in mapped_reads_by_sample_dup.items():
                f.write(sample_name + '\t' + str(mapped_reads) + '\n')

    else:
        info('Getting read counts...')
        bam2reads = get_script_cmdline(cnf, 'perl', join('Seq2C', 'bam2reads.pl'), is_critical=True)

        bam2reads_list_of_bams_fpath = join(cnf.work_dir, 'seq2c_list_of_bams.txt')
        bam2reads_list_of_dup_bams_fpath = join(cnf.work_dir, 'seq2c_list_of_dup_bams.txt')
        with open(bam2reads_list_of_bams_fpath, 'w') as f, open(bam2reads_list_of_dup_bams_fpath, 'w') as fd:
            for sample in bcbio_structure.samples:
                fd.write(sample.name + '\t' + sample.bam + '\n')
                f.write(sample.name + '\t' + dedupped_bam_by_sample[sample.name] + '\n')

        samtools = get_system_path(cnf, 'samtools')
        cmdline = '{bam2reads} -m {samtools} {bam2reads_list_of_bams_fpath}'.format(**locals())
        if not call(cnf, cmdline, mapped_read_fpath): return None, None
        cmdline = '{bam2reads} -m {samtools} {bam2reads_list_of_dup_bams_fpath}'.format(**locals())
        if not call(cnf, cmdline, mapped_read_dup_fpath): return None, None

    return mapped_read_fpath, mapped_read_dup_fpath



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

