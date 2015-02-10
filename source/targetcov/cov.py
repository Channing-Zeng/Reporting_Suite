# coding=utf-8

import sys
import traceback
# import pysam
from collections import OrderedDict, defaultdict
from os.path import join, basename, isfile, abspath, realpath

import source
import source.targetcov
from source.calling_process import call, call_pipe
from source.file_utils import intermediate_fname, splitext_plus, verify_file, file_exists, iterate_file, tmpfile
from source import logger
from source.logger import step_greetings, critical, info, err, warn
from source.reporting import Metric, SampleReport, MetricStorage, ReportSection, PerRegionSampleReport, write_txt_rows, write_tsv_rows
from source.targetcov.Region import Region, save_regions_to_bed, GeneInfo
from source.targetcov.bam_file import index_bam
from source.tools_from_cnf import get_system_path, get_script_cmdline
from source.utils import get_chr_len_fpath


header_metric_storage = MetricStorage(
    general_section=ReportSection('general_section', '', [
        Metric('Target', short_name='Target', common=True),
        Metric('Regions in target', short_name='Regions in target', common=True),
        Metric('Bases in target', short_name='Target bp', unit='bp', common=True),
        Metric('Genes', short_name='Genes', common=True),
        Metric('Genes in target', short_name='Genes in target', common=True),
    ]),
    sections=[
        ReportSection('reads', 'Reads', [
            Metric('Reads'),

            Metric('Mapped reads', short_name='Mapped'),
            Metric('Percentage of mapped reads', short_name='%', unit='%'),

            Metric('Unmapped reads', short_name='Unmapped', quality='Less is better'),
            Metric('Percentage of unmapped reads', short_name='%', unit='%', quality='Less is better'),

            Metric('Duplication rate (picard)', short_name='Dup rate', description='Percent duplication, reported by Picard', quality='Less is better', unit='%'),
        ]),

        ReportSection('target_metrics', 'Target', [
            Metric('Covered bases in target', short_name='Covered in trg', unit='bp'),
            Metric('Percentage of target covered by at least 1 read', short_name='%', unit='%'),

            Metric('Reads mapped on target', short_name='Reads on trg'),
            Metric('Percentage of reads mapped on target', short_name='%', unit='%'),

            Metric('Reads mapped on padded target', 'On padded trg'),
            Metric('Percentage of reads mapped on padded target', short_name='%', unit='%'),

            Metric('Read bases mapped on target', short_name='Read bp on trg', unit='bp'),
        ]),

        ReportSection('depth_metrics', 'Target coverage depth', [
            Metric('Average target coverage depth', short_name='Avg'),
            Metric('Std. dev. of target coverage depth', short_name='Std dev', quality='Less is better'),
            Metric('Maximum target coverage depth', short_name='Max'),
            Metric('Percentage of target within 20% of mean depth', short_name='&#177;20% avg', unit='%', quality='Less is better')
        ]),
    ]
)


def _prep_bed_for_seq2c(cnf, seq2c_bed, amplicons_bed):
    if cnt_fields_in_bed(seq2c_bed) < 4:
        seq2c_bed = amplicons_bed

    elif cnt_fields_in_bed(seq2c_bed) > 4:
        cmdline = 'cut -f1,2,3,4 ' + amplicons_bed
        seq2c_bed = intermediate_fname(cnf, amplicons_bed, 'cut')
        call(cnf, cmdline, amplicons_bed)

    # removing regions with no gene annotation
    def f(l, i):
        if l.split('\t')[3].strip() == '.':
            return None
        else:
            return l
    seq2c_bed = iterate_file(cnf, seq2c_bed, f)

    return seq2c_bed


def _prep_files(cnf, sample, exons_bed):
    if not sample.bam:
        critical(sample.name + ': BAM file is required.')
    if not isfile(sample.bam + '.bai'):
        info('Indexing bam ' + sample.bam)
        index_bam(cnf, sample.bam)

    seq2c_bed = amplicons_bed = sample.bed
    if not sample.bed:
        err(sample.name + ': BED file was not provided. Using AZ-Exome as default: ' + cnf.genome.az_exome)
        seq2c_bed = amplicons_bed = cnf.genome.az_exome

    if not exons_bed:
        critical('Error: no exons specified for the genome in system config.')
    elif abspath(exons_bed) == abspath(amplicons_bed):
        warn('Same file used for exons and amplicons: ' + exons_bed)

    # Exons
    info()
    info('Sorting exons by (chrom, gene name, start); and merging regions within genes...')
    exons_bed = _merge_bed(cnf, exons_bed)

    info()
    info('bedtools-sotring amplicons...')
    amplicons_bed = sort_bed(cnf, amplicons_bed)

    if cnf.reannotate or cnt_fields_in_bed(seq2c_bed) < 4:
        info()
        info('Annotating amplicons with gene names from Ensembl...')
        amplicons_bed = _annotate_amplicons(cnf, amplicons_bed, exons_bed)

    seq2c_bed = _prep_bed_for_seq2c(cnf, seq2c_bed, amplicons_bed)

    info()
    info('Merging amplicons...')
    amplicons_bed = _merge_bed(cnf, amplicons_bed)

    return exons_bed, amplicons_bed, seq2c_bed


def cnt_fields_in_bed(bed_fpath):
    with open(bed_fpath) as f:
        for l in f:
            if l and l.strip() and not l.startswith('#'):
                return len(l.split('\t'))
    critical('Empty bed file: ' + bed_fpath)


def _annotate_amplicons(cnf, amplicons_bed, exons_bed):
    output_fpath = intermediate_fname(cnf, amplicons_bed, 'ann')

    annotate_bed_py = get_system_path(cnf, 'python', join('tools', 'annotate_bed.py'))
    bedtools = get_system_path(cnf, 'bedtools')

    cmdline = '{annotate_bed_py} {amplicons_bed} {cnf.work_dir} {exons_bed} {bedtools}'.format(**locals())
    call(cnf, cmdline, output_fpath)

    return output_fpath


def _filter_bed_with_gene_set(cnf, bed_fpath, gene_names_set):
    def fn(l, i):
        if l:
            fs = l.split('\t')
            new_gns = []
            for g in fs[3].split(','):
                if g in gene_names_set:
                    new_gns.append(g)
            if new_gns:
                return l.replace(fs[3], ','.join(new_gns))

    return iterate_file(cnf, bed_fpath, fn, suffix='key')


def _get_genes_and_filter(cnf, amplicons_bed, exons_bed, genes_fpath):
    gene_names_set = set()
    gene_names_list = list()

    info()
    if genes_fpath:
        with open(genes_fpath) as f:
            gene_names_list = [g.strip() for g in f.read().split('\n') if g]
            gene_names_set = set(gene_names_list)
        info('Using genes from ' + genes_fpath + ', filtering exons and amplicons with this genes.')
        amplicons_bed = _filter_bed_with_gene_set(cnf, amplicons_bed, gene_names_set)
        exons_bed = _filter_bed_with_gene_set(cnf, exons_bed, gene_names_set)
    else:
        gene_names_set, gene_names_list = _get_gene_names(amplicons_bed)
        info('Using genes from amplicons list, filtering exons with this genes.')
        exons_bed = _filter_bed_with_gene_set(cnf, exons_bed, gene_names_set)

    fixed_gene_names_list = []
    added_gene_names_set = set()
    for i in range(len(gene_names_list)):
        gene_name = gene_names_list[i]
        if gene_name not in added_gene_names_set:
            fixed_gene_names_list.append(gene_name)
            added_gene_names_set.add(gene_name)

    return exons_bed, amplicons_bed, gene_names_set, fixed_gene_names_list


class TargetInfo:
    def __init__(self, fpath=None, regions_num=None, bases_num=None, genes_fpath=None, genes_num=None):
        self.fpath = realpath(fpath) if fpath else None
        self.regions_num = regions_num
        self.bases_num = bases_num
        self.genes_fpath = realpath(genes_fpath) if genes_fpath else None
        self.genes_num = genes_num


def make_targetseq_reports(cnf, sample, exons_bed, genes_fpath=None):
    exons_bed, amplicons_bed, seq2c_bed = _prep_files(cnf, sample, exons_bed)

    exons_bed, amplicons_bed, gene_names_set, gene_names_list = \
        _get_genes_and_filter(cnf, amplicons_bed, exons_bed, genes_fpath)

    info()
    info('Calculation of coverage statistics for the regions in the input BED file...')
    amplicons, combined_region, max_depth, total_bed_size = bedcoverage_hist_stats(
        cnf, sample.name, sample.bam, amplicons_bed)

    target_info = TargetInfo(
        fpath=cnf.bed, regions_num=len(amplicons), bases_num=total_bed_size,
        genes_fpath=genes_fpath, genes_num=len(gene_names_list))

    info()
    general_rep_fpath = make_and_save_general_report(cnf, sample, combined_region, max_depth, target_info)

    # Building the gene list
    genes_by_name = OrderedDict()
    for gn in gene_names_list:
        genes_by_name[gn] = GeneInfo(sample_name=sample.name, gene_name=gn)

    un_annotated_amplicons = []

    info('Calculating coverage statistics for exons...')
    exons_with_optional_genes, _, _, _ = bedcoverage_hist_stats(cnf, sample.name, sample.bam, exons_bed)
    for exon_or_gene in exons_with_optional_genes:
        exon_or_gene.sample_name = sample.name

        ef = exon_or_gene.extra_fields
        if ef:
            exon_or_gene.exon_num = ef[0]
            if len(ef) >= 2:
                exon_or_gene.strand = ef[1]
            if len(ef) >= 3:
                exon_or_gene.feature = ef[2]
            else:
                exon_or_gene.feature = 'CDS'
            if len(ef) >= 4:
                exon_or_gene.biotype = ef[3]

        if exon_or_gene.gene_name in genes_by_name:
            gene = genes_by_name[exon_or_gene.gene_name]
            if exon_or_gene.feature == 'gene':
                gene.start = exon_or_gene.start
                gene.end = exon_or_gene.end
                gene.biotype = exon_or_gene.biotype
            else:
                gene.chrom = exon_or_gene.chrom
                gene.strand = exon_or_gene.strand
                gene.add_exon(exon_or_gene)

    for ampl in amplicons:
        ampl.feature = 'Amplicon'
        ampl.sample_name = sample.name

        if ampl.gene_name != '.':
            gene = genes_by_name[ampl.gene_name]
            gene.add_amplicon(ampl)
        else:
            un_annotated_amplicons.append(ampl)

    # per_gene_rep_fpath, genes_by_name = make_and_save_region_report(
    #     cnf, exons_bed_fpath, sample, amplicons, amplicons_bed, gene_names)

    per_gene_rep_fpath = _generate_region_cov_report(cnf, sample, cnf.output_dir, sample.name,
         genes_by_name.values(), un_annotated_amplicons)

    # info()
    # info('Sorting amplicons')
    # amplicons = sorted((a for a in amplicons if a.gene_name), key=Region.get_order_key)

    # info()
    # info('Saving only amplicons overlapped with exons, updated with a gene name')
    # amplicons_bed = save_regions_to_bed(cnf, amplicons, 'targeted_amplicons_with_gene_names',
    #                                     save_original_fields=True)

    info()
    info('Running seq2cov.pl for ' + sample.name)
    seq2c_seq2cov(cnf, sample, seq2c_bed)

    return general_rep_fpath, per_gene_rep_fpath


def seq2c_seq2cov(cnf, sample, amplicons_bed):
    seq2cov = get_script_cmdline(cnf, 'perl', join('Seq2C', 'seq2cov.pl'))
    if not seq2cov: sys.exit(1)

    def fn(l, i): return '\t'.join(l.split('\t')[:4])
    amplicons_bed = iterate_file(cnf, amplicons_bed, fn, suffix='4col')

    seq2c_output = join(
        cnf.output_dir,
        sample.name + '.' + \
        source.targetseq_name + '_' + \
        source.seq2c_seq2cov_ending)
    sample_name = sample.name
    bam = sample.bam
    cmdline = '{seq2cov} -z -b {bam} -N {sample_name} {amplicons_bed}'.format(**locals())
    res = call(cnf, cmdline, seq2c_output)
    if not res:
        err('Could not run seq2cov.pl for ' + sample.name)
        return None


def make_and_save_general_report(cnf, sample, combined_region, max_depth, target_info):
    step_greetings('Target coverage summary report')

    chr_len_fpath = get_chr_len_fpath(cnf)
    ref_fapth = cnf.genome.seq

    summary_report = generate_summary_report(cnf, sample, chr_len_fpath, ref_fapth,
        cnf.coverage_reports.depth_thresholds, cnf.padding, combined_region, max_depth, target_info)

    summary_report_json_fpath = summary_report.save_json(cnf.output_dir, sample.name + '.' + source.targetseq_name)
    summary_report_txt_fpath  = summary_report.save_txt (cnf.output_dir, sample.name + '.' + source.targetseq_name)
    summary_report_html_fpath = summary_report.save_html(cnf.output_dir, sample.name + '.' + source.targetseq_name,
        caption='Target coverage statistics for ' + sample.name)
    info()
    info('Saved to ')
    info('  ' + summary_report_txt_fpath)
    return summary_report_txt_fpath


def _get_gene_names(exons_bed, gene_index=3):
    gene_names_set = set()
    gene_names_list = list()
    # getting gene names for all exons overlapped with amplicons
    with open(exons_bed) as f:
        for line in f:
            if not line or not line.strip() or line.startswith('#'):
                continue

            tokens = line.split()
            if len(tokens) <= gene_index:
                continue

            gene_names_set |= set(tokens[gene_index].split(','))
            gene_names_list.extend(tokens[gene_index].split(','))

    return gene_names_set, gene_names_list


# def make_and_save_region_report(cnf, exons_bed, sample, amplicons, gene_names):
#     step_greetings('Analysing regions.')
#
#     ampli_bed_fpath = save_regions_to_bed(cnf, exons, sample.name + '_exons_by_gene_name')
#     exons_bed_fpath = save_regions_to_bed(cnf, amplicons, sample.name + '_amplicons_by_gene_name')
#
#     info('Groupping exons by gene, getting GeneInfo instances, adding exons to genes...')
#
#     # gene_infos_by_name = _get_exons_combined_by_genes(exons, gene_names)
#
#     info()
#     info('Finding amplicons overlap with exons, adding gene names to amplicons and adding amplicons to genes...')
#
#     amplicons_by_gene_name = _combine_amplicons_by_genes(cnf, sample, ampli_bed_fpath, exons_bed_fpath, gene_names)
#     for gene_name, ampls in amplicons_by_gene_name.items():
#         gene = genes_by_name[gene_name]
#         for a in ampls:
#             gene.add_amplicon(a)
#
#     non_overlapping_exons = [e for g in gene_infos_by_name.values() for e in g.non_overlapping_exons]
#
#     info()
#     non_overlapping_exons_bed_fpath = save_regions_to_bed(cnf, non_overlapping_exons, 'non_overlapping_exons')
#     info()
#     info('Calculating coverage statistics for whole genes, getting Region instances for Genes...')
#     non_overlapping_exons_2, _, _, _ = bedcoverage_hist_stats(cnf, sample.name, sample.bam, non_overlapping_exons_bed_fpath)
#     genes = []
#     for non_overlapping_exon in non_overlapping_exons_2:
#         gene_info = gene_infos_by_name[non_overlapping_exon.gene_name]
#         for d, bs in non_overlapping_exon.bases_by_depth.items():
#             gene_info.bases_by_depth[d] += bs
#         genes.append(gene_info)
#
#     info()
#     info('Building region coverage report.')
#     gene_report_fpath = _generate_region_cov_report(cnf, sample, cnf.output_dir, sample.name, genes)
#
#     return gene_report_fpath, gene_infos_by_name


# def _add_other_exon_of_genes(cnf, gene_names, exons_bed, overlapped_exons_bed):
#     new_overlp_exons_bed = intermediate_fname(cnf, overlapped_exons_bed, 'by_genes')
#     # finding exons for our genes that did not overlapped with amplicons, and adding them too
#     with open(exons_bed) as exons_f, \
#          open(new_overlp_exons_bed, 'w') as new_overl_f:
#         for line in exons_f:
#             if not line or not line.strip() or line.startswith('#') or len(line.split()) < 4:
#                 new_overl_f.write(line)
#             elif line.split()[3] in gene_names:
#                 new_overl_f.write(line)
#
#     return sort_bed(cnf, new_overlp_exons_bed)


def get_records_by_metrics(records, metrics):
    _records = []
    for rec in records:
        if rec.metric.name in metrics:
            rec.metric = metrics[rec.metric.name]
            _records.append(rec)
    return _records


def generate_summary_report(
        cnf, sample, chr_len_fpath, ref_fapth,
        depth_thresholds, padding, combined_region, max_depth, target_info):

    for depth in depth_thresholds:
        name = 'Part of target covered at least by ' + str(depth) + 'x'
        header_metric_storage.add_metric(
            Metric(name, short_name=str(depth) + 'x', description=name, unit='%'),
            'depth_metrics')

    report = SampleReport(sample, metric_storage=header_metric_storage)

    info('* General coverage statistics *')
    info('Getting number of reads...')
    v_number_of_reads = number_of_reads(cnf, sample.bam)
    report.add_record('Reads', v_number_of_reads)

    info('Getting number of mapped reads...')
    v_mapped_reads = number_of_mapped_reads(cnf, sample.bam)
    v_percent_mapped = 100.0 * v_mapped_reads / v_number_of_reads if v_number_of_reads else None
    v_percent_unmapped = 100.0 * (v_number_of_reads - v_mapped_reads) / v_number_of_reads if v_number_of_reads else None
    report.add_record('Mapped reads', v_mapped_reads)
    report.add_record('Percentage of mapped reads', v_percent_mapped)
    report.add_record('Unmapped reads', v_number_of_reads - v_mapped_reads)
    report.add_record('Percentage of unmapped reads', v_percent_unmapped)
    info('')

    info('* Target coverage statistics *')
    report.add_record('Target', target_info.fpath)
    report.add_record('Regions in target', target_info.regions_num)
    report.add_record('Bases in target', target_info.bases_num)
    if target_info.genes_fpath:
        report.add_record('Genes', target_info.genes_fpath)
    report.add_record('Genes in target', target_info.genes_num)

    combined_region.sum_up(depth_thresholds)

    v_covered_bases_in_targ = combined_region.bases_within_threshs.items()[0][1]
    report.add_record('Covered bases in target', v_covered_bases_in_targ)

    v_percent_covered_bases_in_targ = 100.0 * v_covered_bases_in_targ / target_info.bases_num \
        if target_info.bases_num else None
    report.add_record('Percentage of target covered by at least 1 read', v_percent_covered_bases_in_targ)

    info('Getting number of mapped reads on target...')
    v_mapped_reads_on_target = number_mapped_reads_on_target(cnf, sample.bed, sample.bam)
    report.add_record('Reads mapped on target', v_mapped_reads_on_target)

    v_percent_mapped_on_target = 100.0 * v_mapped_reads_on_target / v_mapped_reads if v_mapped_reads else None
    report.add_record('Percentage of reads mapped on target ', v_percent_mapped_on_target)

    info('Making bed file for padded regions...')
    padded_bed = get_padded_bed_file(cnf, sample.bed, chr_len_fpath, padding)

    info('Getting number of mapped reads on padded target...')
    v_reads_on_padded_targ = number_mapped_reads_on_target(cnf, padded_bed, sample.bam)
    report.add_record('Reads mapped on padded target', v_reads_on_padded_targ)

    v_percent_mapped_on_padded_target = 100.0 * v_reads_on_padded_targ / v_mapped_reads if v_mapped_reads else None
    report.add_record('Percentage of reads mapped on padded target', v_percent_mapped_on_padded_target)

    v_read_bases_on_targ = int(target_info.bases_num * combined_region.avg_depth)  # sum of all coverages
    report.add_record('Read bases mapped on target', v_read_bases_on_targ)

    info('')
    report.add_record('Average target coverage depth', combined_region.avg_depth)
    report.add_record('Std. dev. of target coverage depth', combined_region.std_dev)
    report.add_record('Maximum target coverage depth', max_depth)
    report.add_record('Percentage of target within 20% of mean depth', combined_region.percent_within_normal)

    for depth, bases in combined_region.bases_within_threshs.items():
        percent_val = 100.0 * bases / target_info.bases_num if target_info.bases_num else 0
        report.add_record('Part of target covered at least by ' + str(depth) + 'x', percent_val)

    info()

    picard = get_system_path(cnf, 'java', 'picard')
    if picard:
        info('Picard duplication metrics for "' + basename(sample.bam) + '"')
        bam_fpath = sample.bam
        # bam_fpath = _prepare_bam_for_picard(cnf, bam_fpath)
        dup_metrics_txt = join(cnf.work_dir, 'picard_dup_metrics.txt')
        cmdline = '{picard} MarkDuplicates' \
                  ' I={bam_fpath}' \
                  ' O=/dev/null' \
                  ' METRICS_FILE={dup_metrics_txt}' \
                  ' VALIDATION_STRINGENCY=LENIENT'
        cmdline = cmdline.format(**locals())
        call(cnf, cmdline, output_fpath=dup_metrics_txt, stdout_to_outputfile=False, exit_on_error=False)

        if verify_file(dup_metrics_txt, silent=True):
            _parse_picard_dup_report(report, dup_metrics_txt)

        info('Picard ins size hist for "' + basename(sample.bam) + '"')
        picard_ins_size_hist_pdf = join(cnf.output_dir, 'picard_ins_size_hist.pdf')
        picard_ins_size_hist_txt = join(cnf.output_dir, 'picard_ins_size_hist.txt')
        cmdline = '{picard} CollectInsertSizeMetrics' \
                  ' I={bam_fpath}' \
                  ' O={picard_ins_size_hist_txt}' \
                  ' H={picard_ins_size_hist_pdf}' \
                  ' VALIDATION_STRINGENCY=LENIENT'
        cmdline = cmdline.format(**locals())
        call(cnf, cmdline, output_fpath=picard_ins_size_hist_txt, stdout_to_outputfile=False, exit_on_error=False)

    return report


def _prepare_bam_for_picard(cnf, bam_fpath):
    def __process_problem_read_aligns(read_aligns):
        # each alignment: 0:NAME 1:FLAG 2:CHR 3:COORD 4:MAPQUAL 5:CIGAR 6:MATE_CHR 7:MATE_COORD TLEN SEQ ...
        def __get_key(align):
            return align.split('\t')[2] + '@' + align.split('\t')[3]

        def __get_mate_key(align):
            return (align.split('\t')[6] if align.split('\t')[2] != '=' else align.split('\t')[2]) \
                   + '@' + align.split('\t')[7]

        chr_coord = OrderedDict()
        for align in read_aligns:
            key = __get_key(align)
            if key not in chr_coord:
                chr_coord[key] = []
            chr_coord[key].append(align)
        correct_pairs = []
        for align in read_aligns:
            mate_key = __get_mate_key(align)
            if mate_key in chr_coord:
                for pair_align in chr_coord[mate_key]:
                    if read_aligns.index(pair_align) <= read_aligns.index(align):
                        continue
                    if __get_mate_key(pair_align) == __get_key(align):
                        correct_pairs.append((align, pair_align))
        if not correct_pairs:
            return []
        if len(correct_pairs) > 1:
            # sort by sum of mapping quality of both alignments
            correct_pairs.sort(key=lambda pair: pair[0].split('\t')[4] + pair[1].split('\t')[4], reverse=True)
        return [correct_pairs[0][0], correct_pairs[0][1]]

    samtools = get_system_path(cnf, 'samtools')
    ## TODO if REFERENCE_SEQUENCE is needed for Picard (both in CollectInsertSizeMetrics and MarkDuplicates)
    ## problem: header should contain exactly the same chromosomes as in REFERENCE_SEQUENCE
    ## solution:
    ## header_lines = pysam.view("-H", bam_fpath)
    ## modify @SQ entries according to cnf.genome.seq + ".fai" or cnf.genome.chr_lengths
    ## save modified header into COR_HEADER_SAM (plain text file)
    ## samtools reheader COR_HEADER_SAM INPUT_BAM > CORRECTED_BAM

    #TODO: fix alignments for MarkDuplicates (problem more than 2 alignments per read-pair)
    # version 1: sort BAM by queryname, process reads which are present more than twice
    qname_sorted_bam_fpath = intermediate_fname(cnf, bam_fpath, 'qname_sorted')
    cmdline = '{samtools} sort -n -o {bam_fpath} prefix'.format(**locals())  # queryname sorting, output to stdout, 'prefix' is not used
    call(cnf, cmdline, qname_sorted_bam_fpath)
    problem_reads = dict()
    cur_read_aligns = []
    for cur_align in pysam.Samfile(qname_sorted_bam_fpath, 'rb'):
        cur_align = str(cur_align)
        if cur_read_aligns:
            if cur_align.split('\t')[0] != cur_read_aligns[0].split('\t')[0]:
                if len(cur_read_aligns) > 2:
                    problem_reads[cur_read_aligns[0].split('\t')[0]] = cur_read_aligns
                cur_read_aligns = []
        flag = int(cur_align.split('\t')[1])
        cur_read_aligns.append(cur_align)
    if len(cur_read_aligns) > 2:
        problem_reads[cur_read_aligns[0].split('\t')[0]] = cur_read_aligns
    #TODO if speed up is important and disk space is critical
    # version 2: read initial BAM and count number of alignments for each read, process the ones with n>2

    for read_id, read_aligns in problem_reads.items():
        problem_reads[read_id] = __process_problem_read_aligns(read_aligns)
    # correct input BAM
    bam_file = pysam.Samfile(bam_fpath, 'rb')
    fixed_bam_fpath = intermediate_fname(cnf, bam_fpath, 'fixed_for_picard')
    fixed_bam_file = pysam.Samfile(fixed_bam_fpath, "wb", template=bam_file)
    for cur_align in bam_file:
        read_name = str(cur_align).split('\t')[0]
        if read_name not in problem_reads or str(cur_align) in problem_reads[read_name]:
            fixed_bam_file.write(cur_align)
    fixed_bam_file.close()
    return fixed_bam_fpath


def _parse_picard_dup_report(report, dup_report_fpath):
    records = []

    with open(dup_report_fpath) as f:
        for l in f:
            if l.startswith('## METRICS CLASS'):
                try:
                    l_LIBRARY = next(f)
                    l_EMPTY = next(f)
                    l_UNKNOWN = next(f)
                except StopIteration:
                    pass
                else:
                    if l_UNKNOWN:
                        ts = l_UNKNOWN.split()
                        if len(ts) >= 9:
                            dup_rate = 100.0 * float(ts[8])
                            report.add_record('Duplication rate (picard)', dup_rate)
                            return records
    err('Error: cannot read duplication rate from ' + dup_report_fpath)


def _generate_region_cov_report(cnf, sample, output_dir, sample_name, genes, un_annotated_amplicons):
    final_regions = []
    info('Combining all regions for final report...')
    i = 0
    for gene in genes:
        if i and i % 100000 == 0:
            info('Processed {0:,} genes, current gene {1}'.format(i, gene.gene_name))
        i += 1

        final_regions.extend(gene.get_amplicons())
        final_regions.extend(gene.get_exons())
        final_regions.append(gene)

    for ampl in un_annotated_amplicons:
        final_regions.append(ampl)

    info('Processed {0:,} genes.'.format(i))

    info()
    info('Summing up region stats...')
    i = 0
    for region in final_regions:
        i += 1
        if i % 10000 == 0:
            info('Processed {0:,} regions.'.format(i))
        region.sum_up(cnf.coverage_reports.depth_thresholds)

    # info('Sorting...')
    # sorted_regions = sorted(regions, lambda r: (r.start, r.end))
    # info('Saving sorted version...')
    # sorted_regions_fpath = join(cnf.work_dir, (sample.name + '.' + BCBioStructure.detail_sorted_gene_report_baseending))
    # with open(sorted_regions_fpath, 'w') as f:
    #     for r in sorted_regions:
    #         f.write('\t'.join(map(str, [
    #             r.chr, r.start, r.end, r.gene, r.feature,
    #             r.size, r.avg_depth, r.std_dev, r.percent_within_normal])) + '\n')

    info('Saving report...')
    report = _make_flat_region_report(sample, final_regions, cnf.coverage_reports.depth_thresholds)

    gene_report_basename = sample.name + '.' + source.targetseq_name + source.targetcov.detail_gene_report_baseending
    txt_rep_fpath = report.save_txt(output_dir, gene_report_basename)
    tsv_rep_fpath = report.save_tsv(output_dir, gene_report_basename)
    info('')
    info('Regions (total ' + str(len(final_regions)) + ') saved into:')
    info('  ' + txt_rep_fpath)

    return txt_rep_fpath


# def _combine_amplicons_by_genes(cnf, sample, exons_bed_fpath, ampli_bed_fpath, gene_names):
#     exons_ovelaps_with_amplicons_fpath = join(cnf.work_dir, sample.name + '_exons_amplicons_ovelaps.bed')
#     bedtools = get_system_path(cnf, 'bedtools')
#     cmdline = '{bedtools} intersect -wao -a {ampli_bed_fpath} -b {exons_bed_fpath}'.format(**locals())
#     call(cnf, cmdline, output_fpath=exons_ovelaps_with_amplicons_fpath)
#     info()
#
#     amplicons_by_site = defaultdict(list)
#     for a in amplicons:
#         amplicons_by_site[(a.chrom, a.get_start(), a.get_end())].append(a)
#
#     amplicons_by_gene_name = defaultdict(list)
#
#     info()
#     info('Groupping amplicons by genes...')
#     i = 0
#     with open(exons_ovelaps_with_amplicons_fpath) as f:
#         for line in f:
#             a_chrom, a_start, a_end, a_gene_name, a_feature, \
#             e_chrom, e_start, e_end, e_gene_name, e_feature, \
#             overlap_size = line.split('\t')
#
#             if e_gene_name != '.' or a_gene_name != '.':
#                 gene_name = None
#
#                 if e_gene_name != '.':  # hit
#                     if a_gene_name != '.' and a_gene_name != e_gene_name:
#                         err('Amplicon gene name != exon gene name for line: ' + line.strip())
#                     if e_gene_name not in gene_names:
#                         err(e_gene_name + ' from exons not in gene_names from exons and amplicons')
#                         print str(gene_names)
#                         continue
#                     gene_name = e_gene_name
#
#                 elif a_gene_name != '.':  # not hit, but a_gene_name != '.', so amplicons gene names provided
#                     if a_gene_name not in gene_names:
#                         err(a_gene_name + ' from amplicons not in gene_names from exons and amplicons')
#                         continue
#                     gene_name = a_gene_name
#
#                 for a in amplicons_by_site[(a_chrom, int(a_start), int(a_end))]:
#                     amplicons_by_gene_name.append(a)
#
#             if i and i % 10000 == 0:
#                 info('  Processed {0:,} regions, current gene {1}'.format(i, e_gene_name))
#             i += 1
#     info('  Processed {0:,} regions.'.format(i))
#     return amplicons_by_gene_name

    # if exon_gene_summary.intersect(amplicon):
    #     # if amplicon.gene_name and exon_gene.gene_name != amplicon.gene_name:
    #     #     warn('Warning: for ' + str(amplicon) )
    #     amplicon_gene_summary = amplicon_genes_by_name.get(exon_gene_summary.gene_name)
    #     if amplicon_gene_summary is None:
    #         amplicon_gene_summary = SummaryRegion(
    #             sample_name=amplicon.sample_name, chrom=amplicon.chrom,
    #             start=amplicon.start, end=amplicon.end, size=0,
    #             feature='Gene-' + amplicon.feature,
    #             gene_name=exon_gene_summary.gene_name)
    #         amplicon_genes_by_name[amplicon_gene_summary.gene_name] = amplicon_gene_summary
    #     # amplicon_copy = copy.copy(amplicon)
    #     amplicon_gene_summary.add_subregion(amplicon)
    #     # amplicon_copy.gene_name = amplicon_gene_summary.gene_name


def _get_exons_combined_by_genes(exons, ampl_gene_names):
    genes_by_name = OrderedDict()
    # for gn in ampl_gene_names:
    #     genes_by_name[gn] = GeneInfo(sample_name=exon.sample_name, gene_name=exon.gene_name,
    #         chrom=exon.chrom, strand=exon.strand)  # TODO: create genes from amplicons that are not in exons

    i = 0
    for exon in exons:
        if not exon.gene_name:
            info()
            err('  No gene name info in the record: ' + str(exon) + '. Skipping.')
            continue

        if i and i % 10000 == 0:
            info('  Processed {0:,} exons, current gene {1}'.format(i, exon.gene_name))
        i += 1

        gene = genes_by_name.get(exon.gene_name)
        if gene is None:
            gene = GeneInfo(sample_name=exon.sample_name, gene_name=exon.gene_name,
                            chrom=exon.chrom, strand=exon.strand)
            genes_by_name[exon.gene_name] = gene
        gene.add_exon(exon)

    info('  Processed {0:,} exons.'.format(i))
    return genes_by_name


# def _make_region_report(sample, regions, depth_threshs):
#     metrics = [
#         Metric('Sample'),
#         Metric('Chr'),
#         Metric('Start'),
#         Metric('End'),
#         Metric('Gene'),
#         Metric('Feature'),
#         Metric('Size'),
#         Metric('Mean Depth'),
#         Metric('Std Dev', description='Coverage depth standard deviation'),
#         Metric('Wn 20% of Mean Depth', description='Persentage of the region that lies within 20% of an avarage depth.'),
#         Metric('Depth/Median', description='Average depth of coveraege of the region devided by median coverage of the sample'),
#     ]
#     for thres in depth_threshs:
#         metrics.append(Metric('{}x'.format(thres), description='Bases covered by at least {} reads'.format(thres)))
#
#     region_metric_storage = MetricStorage(sections=[ReportSection(metrics=metrics)])
#
#     report = SquareSampleReport(sample, metric_storage=region_metric_storage,
#         report_name='Amplicons and exons coverage statistics.')
#
#     i = 0
#     for region in regions:
#         i += 1
#         if i % 10000 == 0:
#             info('Processed {0:,} regions.'.format(i))
#
#         region.sum_up(depth_threshs)
#         report.add_row('Sample', region.sample)
#
#         # row = map(str, [region.sample,
#         # row = map(str, [region.sample,
#
#         # report
#         row = map(str, [region.sample,
#                         region.chrom,
#                         region.start,
#                         region.end,
#                         region.gene_name])
#         row += [region.feature]
#         row += [str(region.get_size())]
#         row += ['{0:.2f}'.format(region.avg_depth)]
#         row += ['{0:.2f}'.format(region.std_dev)]
#         row += ['{0:.2f}%'.format(region.percent_within_normal) if region.percent_within_normal is not None else '-']
#         if print_cov_factor:
#             row += ['{0:.2f}%'.format(region.cov_factor)]
#
#         for depth_thres, bases in region.bases_within_threshs.items():
#             if int(region.get_size()) == 0:
#                 percent_str = '-'
#             else:
#                 percent = 100.0 * bases / region.get_size()
#                 percent_str = '{0:.2f}%'.format(percent)
#             row.append(percent_str)
#
#         all_rows.append(row)
#         # max_lengths = map(max, izip(max_lengths, chain(map(len, line_fields), repeat(0))))
#     info('Processed {0:,} regions.'.format(i))
#     return all_rows


def _make_flat_region_report(sample, regions, depth_threshs):
    metric_storage = MetricStorage(sections=[ReportSection(metrics=[
        Metric('Sample'),
        Metric('Chr'),
        Metric('Start'),
        Metric('End'),
        Metric('Size'),
        Metric('Gene'),
        Metric('Strand'),
        Metric('Feature'),
        Metric('Biotype'),
        Metric('Min depth'),
        Metric('Ave depth'),
        Metric('Std dev', description='Coverage depth standard deviation'),
        Metric('W/n 20% of ave depth', description='Persentage of the region that lies within 20% of an avarage depth.', unit='%'),
        # Metric('Norm depth', description='Ave region depth devided by median depth of sample'),
    ] + [
        Metric('{}x'.format(thresh), description='Bases covered by at least {} reads'.format(thresh), unit='%')
        for thresh in depth_threshs
    ])])

    # header_fields = ['#Sample', 'Chr', 'Start', 'End', 'Size', 'Gene', 'Strand', 'Feature',
    #                  'Biotype', 'Min Depth', 'Avg Depth', 'Std Dev.', 'Within 20% of Mean']
    # for thres in depth_threshs:
    #     header_fields.append('{}x'.format(thres))
    #
    # all_rows = [header_fields]

    report = PerRegionSampleReport(sample=sample, metric_storage=metric_storage)

    i = 0
    for region in regions:
        i += 1
        if i % 10000 == 0:
            info('Processed {0:,} regions.'.format(i))

        rep_region = report.add_region()
        rep_region.add_record('Sample', region.sample_name)
        rep_region.add_record('Chr', region.chrom)
        rep_region.add_record('Start', region.start)
        rep_region.add_record('End', region.end)
        rep_region.add_record('Size', region.get_size())
        rep_region.add_record('Gene', region.gene_name)
        rep_region.add_record('Strand', region.strand)
        rep_region.add_record('Feature', region.feature)
        rep_region.add_record('Biotype', region.biotype)
        rep_region.add_record('Min depth', region.min_depth)
        rep_region.add_record('Ave depth', region.avg_depth)
        rep_region.add_record('Std dev', region.std_dev)
        rep_region.add_record('W/n 20% of ave depth', region.percent_within_normal)

        for thresh in depth_threshs:
            percent = None
            if region.bases_within_threshs is None:
                err('Error: no bases_within_threshs for ' + str(region))
            else:
                bases = region.bases_within_threshs.get(thresh)
                if bases is not None and region.get_size() > 0:
                    percent = 100.0 * bases / region.get_size()
                    if percent > 100:
                        critical(
                            'Error: percent = ' + str(percent) + ', bases = ' + str(bases) +
                            ', size = ' + str(region.get_size()) +
                            ', start = ' + str(region.start) + ', end = ' + str(region.end))
            rep_region.add_record('{}x'.format(thresh), percent)

        # except:
        #     err('Err in region ' + ' '.join(map(str, [region.sample_name, region.chrom, region.start,
        #                                               region.end, region.gene_name, region.exon_num])))
        #     traceback.print_exc()
        # else:

        # max_lengths = map(max, izip(max_lengths, chain(map(len, line_fields), repeat(0))))

    info('Processed {0:,} regions.'.format(i))
    return report


def bedcoverage_hist_stats(cnf, sample_name, bam, bed, reuse=False):
    if not bam or not bed:
        err()
        if not bam: err('BAM file is required.')
        if not bed: err('BED file is required.')
        sys.exit(1)

    regions, max_depth, total_bed_size = [], 0, 0

    bedtools = get_system_path(cnf, 'bedtools')
    cmdline = '{bedtools} coverage -abam {bam} -b {bed} -hist'.format(**locals())
    bedcov_output = join(cnf.work_dir,
        splitext_plus(basename(bed))[0] + '_' +
        splitext_plus(basename(bam))[0] + '_bedcov_output.txt')
    if reuse and file_exists(bedcov_output) and verify_file(bedcov_output):
        pass
    else:
        call(cnf, cmdline, bedcov_output)

    _total_regions_count = 0

    info('Anylising bedcoverage output...')
    with open(bedcov_output) as f:
        for next_line in f:
            if not next_line.strip() or next_line.startswith('#'):
                continue

            line_tokens = next_line.strip().split()
            chrom = line_tokens[0]
            start, end, gene_name = None, None, None
            try:
                depth, bases, region_size = map(int, line_tokens[-4:-1])
            except:
                critical('Undexpected error: incorrect line in the coverageBed output:\n' + next_line)
                return

            if next_line.startswith('all'):
                max_depth = max(max_depth, depth)
                total_bed_size += bases
                extra_fields = ()
            else:
                start, end = map(int, line_tokens[1:3])
                gene_name = line_tokens[3]
                extra_fields = tuple(line_tokens[4:-4])

            line_region_key_tokens = (sample_name, chrom, start, end, gene_name, extra_fields)

            if regions == [] or hash(line_region_key_tokens) != \
                    hash((regions[-1].sample_name, regions[-1].chrom,
                          regions[-1].start, regions[-1].end,
                          regions[-1].gene_name, regions[-1].extra_fields)):
                region = Region(
                    sample_name=sample_name, chrom=chrom,
                    start=start, end=end, size=region_size,
                    gene_name=gene_name, extra_fields=extra_fields)
                regions.append(region)

                _total_regions_count += 1

                if _total_regions_count > 0 and _total_regions_count % 100000 == 0:
                    info('  Processed {0:,} regions'.format(_total_regions_count))

            regions[-1].add_bases_for_depth(depth, bases)

            if regions[-1].min_depth is None:
                regions[-1].min_depth = depth

    if _total_regions_count % 100000 != 0:
        info('  Processed {0:,} regions'.format(_total_regions_count))

    total_region = regions[-1]
    regions = regions[:-1]
    # info('Sorting genes...')
    # regions = sorted(regions[:-1], key=Region.get_order_key)

    return regions, total_region, max_depth, total_bed_size


def _merge_bed(cnf, bed_fpath):
    output_fpath = intermediate_fname(cnf, bed_fpath, 'merge')

    merge_bed_py = get_system_path(cnf, 'python', join('tools', 'merge_bed.py'))

    cmdline = '{merge_bed_py} {bed_fpath}'.format(**locals())

    call(cnf, cmdline, output_fpath)

    return output_fpath


# def _fix_amplicons_gene_names(cnf, amplicons_fpath):
#     output_fpath = intermediate_fname(cnf, amplicons_fpath, 'fixedgenenames')
#
#     with open(amplicons_fpath) as f, open(output_fpath, 'w') as out:
#         prev_end = None
#         prev_chr = None
#
#         for line in f:
#             if not line.strip() or line.startswith('#'):
#                 out.write(line)
#                 continue
#
#             ts = line.split()
#             assert len(ts) >= 3
#
#             # cur_chr = ts[0]
#             # if prev_chr is not None and prev_chr == cur_chr and prev_end is not None:
#             #     cur_start = int(ts[1])
#                 # if prev_end > cur_start:
#                 #     err(line + ': prev region end ' + str(prev_end) + ' is more then current start ' + str(cur_start))
#             # prev_end = int(ts[2])
#             # prev_chr = cur_chr
#
#             if len(ts) >= 4:
#                 if ':' in ts[3]:
#                     ts[3] = ts[3].split(':')[-1]
#                 if len(ts) < 8:
#                     ts = ts[:4]
#                 if len(ts) > 8:
#                     ts = ts[:8]
#             out.write('\t'.join(ts) + '\n')
#
#     info('Saved to ' + output_fpath)
#     return output_fpath


def sort_bed(cnf, bed_fpath):
    bedtools = get_system_path(cnf, 'bedtools')
    cmdline = '{bedtools} sort -i {bed_fpath}'.format(**locals())
    output_fpath = intermediate_fname(cnf, bed_fpath, 'sorted')
    call(cnf, cmdline, output_fpath)
    return output_fpath


def _unique_longest_exons(cnf, exons_bed_fpath):
    unique_exons_dict = OrderedDict()

    with open(exons_bed_fpath) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue

            ts = line.split()

            if len(ts) < 4:
                pass

            elif len(ts) < 5:
                chrom, start, end, gene = ts
                unique_exons_dict[(gene, '')] = ts

            else:
                chrom, start, end, gene, exon_num = ts[:5]
                prev_ts = unique_exons_dict.get((gene, exon_num))
                if not prev_ts:
                    unique_exons_dict[(gene, exon_num)] = ts
                else:
                    size = int(ts[2]) - int(ts[1])
                    prev_size = int(prev_ts[2]) - int(prev_ts[1])
                    if size > prev_size:
                        unique_exons_dict[(gene, exon_num)] = ts

    unique_bed_fpath = intermediate_fname(cnf, exons_bed_fpath, 'uniq')
    with open(unique_bed_fpath, 'w') as f:
        for ts in unique_exons_dict.values():
            f.write('\t'.join(ts) + '\n')

    info('Saved to ' + unique_bed_fpath)
    return unique_bed_fpath


def intersect_bed(cnf, bed1, bed2):
    bed1_fname, _ = splitext_plus(basename(bed1))
    bed2_fname, _ = splitext_plus(basename(bed2))
    output_fpath = join(cnf['work_dir'], bed1_fname + '__' + bed2_fname + '.bed')
    bedtools = get_system_path(cnf, 'bedtools')
    cmdline = '{bedtools} intersect -u -a {bed1} -b {bed2}'.format(**locals())
    call(cnf, cmdline, output_fpath)
    return output_fpath


#TODO: works slow.
def number_of_reads(cnf, bam):
    samtools = get_system_path(cnf, 'samtools')
    output_fpath = join(cnf.work_dir, 'num_reads')
    if not isfile(output_fpath):  # TODO: tmp
        cmdline = '{samtools} view -c {bam}'.format(**locals())
        call(cnf, cmdline, output_fpath)
    with open(output_fpath) as f:
        return int(f.read().strip())


#TODO: works slow.
def number_of_mapped_reads(cnf, bam):
    samtools = get_system_path(cnf, 'samtools')
    output_fpath = join(cnf.work_dir, 'num_mapped_reads')
    if not isfile(output_fpath):  # TODO: tmp
        cmdline = '{samtools} view -c -F 4 {bam}'.format(**locals())
        call(cnf, cmdline, output_fpath)
    with open(output_fpath) as f:
        return int(f.read().strip())


#TODO: works slow.
def number_of_unmapped_reads(cnf, bam):
    samtools = get_system_path(cnf, 'samtools')
    output_fpath = join(cnf.work_dir, 'num_unmapped_reads')
    if not isfile(output_fpath):  # TODO: tmp
        cmdline = '{samtools} view -c -f 4 {bam}'.format(**locals())
        call(cnf, cmdline, output_fpath)
    with open(output_fpath) as f:
        return int(f.read().strip())


#TODO: works slow.
def number_mapped_reads_on_target(cnf, bed, bam):
    samtools = get_system_path(cnf, 'samtools')
    output_fpath = join(cnf.work_dir, 'num_mapped_reads_target')
    if not isfile(output_fpath):  # TODO: tmp
        cmdline = '{samtools} view -c -F 4 -L {bed} {bam}'.format(**locals())
        call(cnf, cmdline, output_fpath)
    with open(output_fpath) as f:
        return int(f.read().strip())


#TODO: works slow.
def number_bases_in_aligned_reads(cnf, bam):
    samtools = get_system_path(cnf, 'samtools')
    cmdline = '{samtools} depth {bam}'.format(**locals())
    proc = call_pipe(cnf, cmdline)
    count = 0
    while True:
        coverage_line = proc.stdout.readline()
        if coverage_line:
            values = coverage_line.strip().split('\t')
            count += int(values[2])
    return count


def get_padded_bed_file(cnf, bed, genome, padding):
    bedtools = get_system_path(cnf, 'bedtools')
    cmdline = '{bedtools} slop -i {bed} -g {genome} -b {padding}'.format(**locals())
    output_fpath = intermediate_fname(cnf, bed, 'padded')
    call(cnf, cmdline, output_fpath)
    return output_fpath


def _bases_by_depth(depth_vals, depth_thresholds):
    bases_by_min_depth = {depth: 0 for depth in depth_thresholds}

    for depth_value in depth_vals:
        for threshold in depth_thresholds:
            if depth_value >= threshold:
                bases_by_min_depth[threshold] += 1

        return [100.0 * bases_by_min_depth[thres] / len(depth_vals) if depth_vals else 0
                for thres in depth_thresholds]


