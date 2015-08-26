# coding=utf-8

from collections import OrderedDict, defaultdict
from os.path import join, basename, isfile, abspath, realpath, splitext, normpath
import shutil

import source
from source.bcbio_structure import BCBioStructure
from source.qsub_utils import submit_job, wait_for_jobs
from source.qualimap.report_parser import parse_qualimap_sample_report
import source.targetcov
from source.calling_process import call, call_pipe
from source.file_utils import intermediate_fname, splitext_plus, verify_file, file_exists, iterate_file, tmpfile, \
    safe_mkdir, add_suffix
from source.logger import step_greetings, critical, info, err, warn
from source.reporting import Metric, SampleReport, MetricStorage, ReportSection, PerRegionSampleReport, write_txt_rows, write_tsv_rows, \
    load_records
from source.targetcov.Region import Region, save_regions_to_bed, GeneInfo, calc_bases_within_threshs
from source.targetcov.bam_and_bed_utils import index_bam, prepare_beds, filter_bed_with_gene_set, get_total_bed_size, \
    total_merge_bed, count_bed_cols, annotate_amplicons, group_and_merge_regions_by_gene, sort_bed, fix_bed_for_qualimap, \
    remove_dups
from source.targetcov.coverage_hist import bedcoverage_hist_stats
from source.tools_from_cnf import get_system_path, get_script_cmdline
from source.utils import get_chr_len_fpath


def get_header_metric_storage(depth_thresholds):
    ms = MetricStorage(
        general_section=ReportSection('general_section', '', [
            Metric('Target', short_name='Target', common=True),
            # Metric('Target ready', short_name='Target ready', common=True),
            Metric('Regions in target', short_name='Regions in target', common=True),
            Metric('Bases in target', short_name='Target bp', unit='bp', common=True),
            Metric('Genes', short_name='Genes', common=True),
            Metric('Genes in target', short_name='Genes in target', common=True),
        ]),
        sections=[
            ReportSection('reads', 'Reads', [
                Metric('Reads'),

                Metric('Mapped reads', short_name='Mapped', description='samtools view -c -F 4'),
                Metric('Percentage of mapped reads', short_name='%', unit='%'),

                Metric('Unmapped reads', short_name='Unmapped', quality='Less is better', description='samtools view -c -f 4'),
                Metric('Percentage of unmapped reads', short_name='%', unit='%', quality='Less is better'),

                Metric('Properly paired reads percent', short_name='Paired %', unit='%', description='Pecent of properly paired mapped reads (-f 2).'),

                Metric('Duplication rate', short_name='Dup rate', description='Percent of mapped reads (-F 4), marked as duplicates (-f 1024)', quality='Less is better', unit='%'),
                Metric('Dedupped mapped reads', short_name='Dedupped', description='Mapped reads (-F 4), not makred as duplicates (-f 1024)'),
            ]),

            ReportSection('target_metrics', 'Target (duplicate reads are not counted)', [
                Metric('Covered bases in target', short_name='Trg covered', unit='bp'),
                Metric('Percentage of target covered by at least 1 read', short_name='%', unit='%'),

                Metric('Reads mapped on target', short_name='Reads on trg'),
                Metric('Percentage of reads mapped on target', short_name='%', unit='%'),

                Metric('Percentage of reads mapped off target', short_name='% off trg', unit='%', quality='Less is better'),

                Metric('Reads mapped on padded target', 'On padded trg'),
                Metric('Percentage of reads mapped on padded target', short_name='%', unit='%'),

                Metric('Read bases mapped on target', short_name='Read bp on trg', unit='bp'),
            ]),

            ReportSection('depth_metrics', 'Target coverage depth (duplicate reads are not counted)', [
                Metric('Average target coverage depth', short_name='Avg'),
                Metric('Std. dev. of target coverage depth', short_name='Std dev', quality='Less is better'),
                Metric('Maximum target coverage depth', short_name='Max'),
                Metric('Percentage of target within 20% of mean depth', short_name='&#177;20% avg', unit='%')
            ]),
        ])

    for depth in depth_thresholds:
        name = 'Part of target covered at least by ' + str(depth) + 'x'
        ms.add_metric(
            Metric(name, short_name=str(depth) + 'x', description=name, unit='%'),
            'depth_metrics')

    return ms


def _extract_gene_names_and_filter_exons(cnf, sample_name, target_bed, exons_bed, exons_no_genes_bed, genes_fpath):
    gene_by_name = OrderedDict()

    gene_names_set = set()
    gene_names_list = []

    info()
    info('Getting gene list')
    if genes_fpath:
        with open(genes_fpath) as f:
            gene_names_list = [g.strip() for g in f.read().split('\n') if g]
            gene_names_set = set(gene_names_list)
        info('Using genes from ' + genes_fpath + ', filtering exons and amplicons with this genes.')
        if target_bed:
            target_bed = filter_bed_with_gene_set(cnf, target_bed, gene_names_set)
        if exons_bed:
            exons_bed = filter_bed_with_gene_set(cnf, exons_bed, gene_names_set)
            exons_no_genes_bed = filter_bed_with_gene_set(cnf, exons_no_genes_bed, gene_names_set)

    elif target_bed or exons_bed:
        info()
        gene_names_set, gene_names_list = _get_gene_names(target_bed or exons_bed)
        info('Using genes from amplicons list, trying filtering exons with this genes.')
        if target_bed and exons_bed:
            exons_anno_bed = filter_bed_with_gene_set(cnf, exons_bed, gene_names_set)
            if not verify_file(exons_anno_bed):
                info()
                warn('No gene symbols from the capture bed file was found in Ensemble. Re-annotating target...')
                target_bed = annotate_amplicons(cnf, target_bed, exons_bed)
                info('Merging regions within genes...')
                target_bed = group_and_merge_regions_by_gene(cnf, target_bed, keep_genes=False)
                info('Sorting amplicons_bed by (chrom, gene name, start)')
                target_bed = sort_bed(cnf, target_bed, cnf.genome.name)
                info('Getting gene names again...')
                gene_names_set, gene_names_list = _get_gene_names(target_bed)
                info()
                info('Using genes from amplicons list, filtering exons with this genes.')
                exons_anno_bed = filter_bed_with_gene_set(cnf, exons_bed, gene_names_set)
                if not verify_file(exons_anno_bed):
                    critical('No gene symbols from the capture bed file was found in Ensemble.')
            exons_bed = exons_anno_bed
            exons_no_genes_bed = filter_bed_with_gene_set(cnf, exons_no_genes_bed, gene_names_set)

    info()

    info('Making unique gene list without affecting the order')
    fixed_gene_names_list = []
    added_gene_names_set = set()
    for i in range(len(gene_names_list)):
        gene_name = gene_names_list[i]
        if gene_name not in added_gene_names_set:
            fixed_gene_names_list.append(gene_name)
            added_gene_names_set.add(gene_name)
    gene_names_list = fixed_gene_names_list
    info('Uniq gene list contains ' + str(len(gene_names_list)) + ' genes')

    info('Building the gene list')
    for gn in gene_names_list:
        gene_by_name[gn] = GeneInfo(sample_name=sample_name, gene_name=gn)
    info('Processed ' + str(len(gene_names_list)) + ' gene records -> ' + str(len(gene_by_name)) + ' uniq gene sybmols')

    if exons_bed:
        info()
        info('Filtering exon bed file to have only gene records...')
        exons_only_genes_bed = intermediate_fname(cnf, exons_bed, 'only_genes')
        call(cnf, 'grep -w Gene ' + exons_bed, output_fpath=exons_only_genes_bed)
        info('Saved genes to ' + exons_only_genes_bed)

        info()
        info('Setting start and end for the genes')
        i = 0
        with open(exons_only_genes_bed) as f:
            for l in f:
                l = l.strip()
                if l and not l.startswith('#'):
                    fs = l.split('\t')
                    chrom, start, end, symbol = fs[:4]
                    gene_by_name[symbol].start = int(start)
                    gene_by_name[symbol].end = int(end)
                    if len(fs) >= 8:
                        gene_by_name[symbol].biotype = fs[7]
                    i += 1
        info('Processed ' + str(i) + ' genes')
        info()

    return target_bed, exons_bed, exons_no_genes_bed, gene_by_name


class TargetInfo:
    def __init__(self, fpath=None, bed=None, original_target_bed=None, regions_num=None,
                 bases_num=None, genes_fpath=None, genes_num=None):
        self.fpath = realpath(fpath) if fpath else None  # raw source file - to demonstrate where we took one
        self.bed = bed                                   # processed (sorted, merged...), to do real calculations
        self.original_target_bed = original_target_bed
        self.regions_num = regions_num
        self.bases_num = bases_num
        self.genes_fpath = realpath(genes_fpath) if genes_fpath else None
        self.genes_num = genes_num


def _run_qualimap(cnf, sample, bam_fpath, bed_fpath=None):
    safe_mkdir(sample.qualimap_dirpath)

    bed = ''
    if bed_fpath:
        qualimap_bed_fpath = join(cnf.work_dir, 'tmp_qualimap.bed')
        fix_bed_for_qualimap(bed_fpath, qualimap_bed_fpath)
        bed = '--bed ' + qualimap_bed_fpath

    qm = get_system_path(cnf, 'python', join('scripts', 'post', 'qualimap.py'))
    cmdl = '{qm} --bam {bam_fpath} {bed} -o {output_dirpath} -t {cnf.threads}'.format(**locals())
    call(cnf, cmdl, sample.qualimap_dirpath, output_is_dir=True)
    return sample.qualimap_dirpath

    # join(self.qualimap, sample_name, bam=sample.bam, sample=sample.name,
    #     genome=sample.genome, bed=qualimap_gff_cmdl, threads=self.threads_per_sample)
 # hist = qualimap_cov_hist_fpath
 #
 #
 #                self.qualimap = Step(cnf, run_id,
 #            name=BCBioStructure.qualimap_name, short_name='qm',
 #            interpreter='python',
 #            script=join('scripts', 'post', 'qualimap.py'),
 #            dir_name=BCBioStructure.qualimap_dir,
 #            paramln=params_for_one_sample + ' --bam {bam} {bed} -o {output_dir}',
 #        )


def _dedup_and_flag_stat(cnf, bam_fpath):
    bam_stats = samtools_flag_stat(cnf, bam_fpath)
    info('Total reads: ' + Metric.format_value(bam_stats['total']))
    info('Total mapped reads: ' + Metric.format_value(bam_stats['mapped']))
    info('Total dup reads: ' + Metric.format_value(bam_stats['duplicates']))
    info('Total properly paired reads: ' + Metric.format_value(bam_stats['properly paired']))

    dedup_bam_dirpath = join(cnf.work_dir, source.dedup_bam)
    safe_mkdir(dedup_bam_dirpath)
    dedup_bam_fpath = join(dedup_bam_dirpath, add_suffix(basename(bam_fpath), source.dedup_bam))
    remove_dups(cnf, bam_fpath, output_fpath=dedup_bam_fpath, use_grid=False)
    dedup_bam_stats = samtools_flag_stat(cnf, dedup_bam_fpath)
    info('Total reads after dedup (samtools view -F 1024): ' + Metric.format_value(dedup_bam_stats['total']))
    info('Total mapped reads after dedup (samtools view -F 1024): ' + Metric.format_value(dedup_bam_stats['mapped']))
    return dedup_bam_fpath, bam_stats, dedup_bam_stats


def _parse_qualimap_results(qualimap_html_fpath, qualimap_cov_hist_fpath, depth_thresholds):
    if not verify_file(qualimap_cov_hist_fpath):
        err('Qualimap hist fpath is not found, cannot build histogram')
    else:
        bases_by_depth = OrderedDict()
        with open(qualimap_cov_hist_fpath) as f:
            for l in f:
                if l.startswith('#'):
                    pass
                else:
                    cov, bases = map(int, l.strip().split())
                    bases_by_depth[cov] = bases

        max_depth = bases_by_depth.items()[-1][0]
        n = len(bases_by_depth)

    if not verify_file(qualimap_html_fpath):
        err('Qualimap report was not found')
    else:
        records = parse_qualimap_sample_report(qualimap_html_fpath)
        ave_depth = next(r for r in records if r.metric.name == 'Coverage Mean')
        stddev_depth = next(r for r in records if r.metric.name == 'Coverage Standard Deviation')

    return bases_by_depth, max_depth, ave_depth, stddev_depth, within20percent_ofmean


def make_targetseq_reports(cnf, output_dir, sample, bam_fpath, exons_bed, exons_no_genes_bed, target_bed, genes_fpath=None):
    info('Starting targeqSeq for ' + sample.name + ', saving into ' + output_dir)
    original_target_bed = cnf.original_target_bed or target_bed
    original_exons_bed = cnf.original_exons_bed or exons_no_genes_bed

    bam_fpath, bam_stats, dedup_bam_stats = _dedup_and_flag_stat(cnf, bam_fpath)

    target_info = None
    targets = []
    if target_bed or exons_bed:
        target_bed, exons_bed, exons_no_genes_bed, gene_by_name = \
            _extract_gene_names_and_filter_exons(cnf, sample.name, target_bed, exons_bed, exons_no_genes_bed, genes_fpath)

        # info()
        # info('Calculation of coverage statistics for the regions in the input BED file...')
        # targets, _, max_depth = bedcoverage_hist_stats(cnf, sample.name, bam_fpath, target_bed or exons_no_genes_bed)
        # info()
        # info('Merging capture BED file to get total target cov statistics...')
        total_merged_target_bed = total_merge_bed(cnf, target_bed or exons_no_genes_bed)
        # info()
        # info('Calculation of coverage statistics for total target...')
        # _, combined_region, _ = bedcoverage_hist_stats(cnf, sample.name, bam_fpath, total_merged_target_bed)

        total_bed_size = get_total_bed_size(cnf, target_bed or exons_no_genes_bed)

        ready_target_bed = join(output_dir, 'target.bed')
        shutil.copy(target_bed or exons_bed, ready_target_bed)

        target_info = TargetInfo(
            fpath=target_bed, bed=target_bed, original_target_bed=original_target_bed,
            regions_num=len(targets), bases_num=total_bed_size,
            genes_fpath=genes_fpath, genes_num=len(gene_by_name))

    info()
    summary_report = make_and_save_general_report(cnf, output_dir, sample, bam_fpath, bam_stats, dedup_bam_stats, target_info)

    un_annotated_amplicons = []

    if exons_bed:
        if not target_bed:
            exons_with_optional_genes = targets
        else:
            info()
            info('Calculating coverage statistics for exons...')
            exons_with_optional_genes, _, _ = bedcoverage_hist_stats(cnf, sample.name, bam_fpath, exons_no_genes_bed)

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

            if exon_or_gene.gene_name in gene_by_name:
                gene = gene_by_name[exon_or_gene.gene_name]
                gene.chrom = exon_or_gene.chrom
                gene.strand = exon_or_gene.strand
                gene.add_exon(exon_or_gene)

        if target_bed:
            for ampl in targets:
                ampl.feature = 'Capture'
                ampl.sample_name = sample.name

                if ampl.gene_name != '.':
                    gene = gene_by_name[ampl.gene_name]
                    gene.add_amplicon(ampl)
                else:
                    un_annotated_amplicons.append(ampl)

            un_annotated_amplicons = sorted(un_annotated_amplicons, key=lambda r: (r.start, r.end))

    per_gene_report = None
    if exons_bed or target_bed:
        per_gene_report = _generate_region_cov_report(cnf, sample, output_dir,
            gene_by_name.values(), un_annotated_amplicons)

    info()
    return combined_region.avg_depth, gene_by_name, [summary_report, per_gene_report]


def make_and_save_general_report(cnf, output_dir, sample, bam_fpath, bam_stats, dedup_bam_stats, target_info=None):
    step_greetings('Target coverage summary report')

    chr_len_fpath = get_chr_len_fpath(cnf)
    ref_fapth = cnf.genome.seq

    summary_report = generate_summary_report(cnf, output_dir, sample, bam_fpath, chr_len_fpath, ref_fapth,
        bam_stats, dedup_bam_stats,
        cnf.coverage_reports.depth_thresholds, cnf.padding, target_info=target_info)

    summary_report.save_json(output_dir, sample.name + '.' + source.targetseq_name)
    summary_report.save_txt (output_dir, sample.name + '.' + source.targetseq_name)
    summary_report.save_html(output_dir, sample.name + '.' + source.targetseq_name, caption='Target coverage statistics for ' + sample.name)
    info()
    info('Saved to ')
    info('  ' + summary_report.txt_fpath)
    return summary_report


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

            for gn in tokens[gene_index].split(','):
                if gn != '.':
                    gene_names_set.add(gn)
                    gene_names_list.append(gn)

    return gene_names_set, gene_names_list


def get_records_by_metrics(records, metrics):
    _records = []
    for rec in records:
        if rec.metric.name in metrics:
            rec.metric = metrics[rec.metric.name]
            _records.append(rec)
    return _records


def generate_summary_report(
        cnf, output_dir, sample, bam_fpath, chr_len_fpath, ref_fapth,
        target_bed, exons_no_genes_bed,
        bam_stats, dedup_bam_stats,
        depth_thresholds, padding, target_info=None):  # TODO: calculate max_depth independently of target

    report = SampleReport(sample, metric_storage=get_header_metric_storage(depth_thresholds))
    info('* General coverage statistics *')
    report.add_record('Reads', bam_stats['total'])
    report.add_record('Mapped reads', bam_stats['mapped'])
    report.add_record('Unmapped reads', bam_stats['total'] - bam_stats['mapped'])
    percent_mapped = 1.0 * bam_stats['mapped'] / bam_stats['total'] if bam_stats['total'] else None
    assert percent_mapped <= 1.0 or percent_mapped is None, str(percent_mapped)
    report.add_record('Percentage of mapped reads', percent_mapped)
    percent_unmapped = 1.0 * (bam_stats['total'] - bam_stats['mapped']) / bam_stats['total'] if bam_stats['total'] else None
    assert percent_unmapped <= 1.0 or percent_unmapped is None, str(percent_unmapped)
    report.add_record('Percentage of unmapped reads', percent_unmapped)
    total_paired_reads_pecent = 1.0 * bam_stats['properly paired'] / bam_stats['total'] if bam_stats['total'] else None
    assert total_paired_reads_pecent <= 1.0 or total_paired_reads_pecent is None, str(total_paired_reads_pecent)
    report.add_record('Properly paired reads percent', total_paired_reads_pecent)
    info('')

    if dedup_bam_stats:
        dup_rate = 1 - (1.0 * dedup_bam_stats['mapped'] / bam_stats['mapped']) if bam_stats['mapped'] else None
        report.add_record('Duplication rate', dup_rate)
        report.add_record('Dedupped mapped reads', dedup_bam_stats['mapped'])

    _run_qualimap(cnf, sample.name, bam_fpath, target_bed or exons_no_genes_bed)
    bases_by_depth, max_depth, ave_depth, stddev_depth, within20percent_ofmean = \
        _parse_qualimap_results(sample.qualimap_html_fpath, sample.qualimap_cov_hist_fpath, depth_thresholds)
    bases_within_threshs, rates_within_threshs = calc_bases_within_threshs(bases_by_depth, depth_thresholds)

    if target_info:
        info('* Target coverage statistics *')
        if target_info.original_target_bed:
            report.add_record('Target', target_info.original_target_bed)
            # report.add_record('Target ready', target_info.fpath)
        else:
            report.add_record('Target', target_info.fpath)

        report.add_record('Regions in target', target_info.regions_num)
        report.add_record('Bases in target', target_info.bases_num)
        if target_info.genes_fpath: report.add_record('Genes', target_info.genes_fpath)
        report.add_record('Genes in target', target_info.genes_num)

        v_covered_bases_in_targ = bases_within_threshs.items()[0][1]
        report.add_record('Covered bases in target', v_covered_bases_in_targ)
        v_percent_covered_bases_in_targ = 1.0 * v_covered_bases_in_targ / target_info.bases_num if target_info.bases_num else None
        report.add_record('Percentage of target covered by at least 1 read', v_percent_covered_bases_in_targ)
        assert v_percent_covered_bases_in_targ <= 1.0 or v_percent_covered_bases_in_targ is None, str(v_percent_covered_bases_in_targ)

        info('Getting number of mapped reads on target...')
        mapped_reads_on_target = number_mapped_reads_on_target(cnf, target_info.bed, bam_fpath)
        report.add_record('Reads mapped on target', mapped_reads_on_target)
        if (dedup_bam_stats or bam_stats)['mapped']:
            percent_mapped_on_target = 1.0 * mapped_reads_on_target / (dedup_bam_stats or bam_stats)['mapped']
            report.add_record('Percentage of reads mapped on target', percent_mapped_on_target)
            assert percent_mapped_on_target <= 1.0 or percent_mapped_on_target is None, str(percent_mapped_on_target)
            percent_mapped_off_target = 1.0 - percent_mapped_on_target
            report.add_record('Percentage of reads mapped off target ', percent_mapped_off_target)

        info('Making bed file for padded regions...')
        padded_bed = get_padded_bed_file(cnf, target_info.bed, chr_len_fpath, padding)
        info('Getting number of dedupped mapped reads on padded target...')
        reads_on_padded_targ = number_mapped_reads_on_target(cnf, padded_bed, bam_fpath)
        report.add_record('Reads mapped on padded target', reads_on_padded_targ)
        percent_mapped_on_padded_target = 1.0 * reads_on_padded_targ / (dedup_bam_stats or bam_stats)['mapped'] if (dedup_bam_stats or bam_stats)['mapped'] else None
        report.add_record('Percentage of reads mapped on padded target', percent_mapped_on_padded_target)
        assert percent_mapped_on_padded_target <= 1.0 or percent_mapped_on_padded_target is None, str(percent_mapped_on_padded_target)

        read_bases_on_targ = int(target_info.bases_num * combined_region.avg_depth)  # sum of all coverages
        report.add_record('Read bases mapped on target', read_bases_on_targ)

        info('')
        report.add_record('Average target coverage depth', combined_region.avg_depth)
        report.add_record('Std. dev. of target coverage depth', combined_region.std_dev)
        report.add_record('Maximum target coverage depth', max_depth)
        report.add_record('Percentage of target within 20% of mean depth', combined_region.rate_within_normal)
        assert combined_region.rate_within_normal <= 1.0 or combined_region.rate_within_normal is None, str(combined_region.rate_within_normal)

        for depth, bases in combined_region.bases_within_threshs.items():
            percent_val = 1.0 * bases / target_info.bases_num if target_info.bases_num else 0
            report.add_record('Part of target covered at least by ' + str(depth) + 'x', percent_val)
            assert percent_val <= 1.0 or percent_val is None, str(percent_val)

        info()

    picard = get_system_path(cnf, 'java', 'picard')
    if picard:
        info()
        info('Picard ins size hist for "' + basename(bam_fpath) + '"')
        picard_ins_size_hist_pdf = join(output_dir, 'picard_ins_size_hist.pdf')
        picard_ins_size_hist_txt = join(output_dir, 'picard_ins_size_hist.txt')
        cmdline = '{picard} CollectInsertSizeMetrics' \
                  ' I={bam_fpath}' \
                  ' O={picard_ins_size_hist_txt}' \
                  ' H={picard_ins_size_hist_pdf}' \
                  ' VALIDATION_STRINGENCY=LENIENT'

        cmdline = cmdline.format(**locals())
        call(cnf, cmdline, output_fpath=picard_ins_size_hist_txt, stdout_to_outputfile=False, exit_on_error=False)

    return report


def _generate_region_cov_report(cnf, sample, output_dir, genes, un_annotated_amplicons):
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

    info('Saving report...')
    report = make_flat_region_report(sample, final_regions, cnf.coverage_reports.depth_thresholds)

    gene_report_basename = sample.name + '.' + source.targetseq_name + source.detail_gene_report_baseending
    report.save_txt(output_dir, gene_report_basename)
    report.save_tsv(output_dir, gene_report_basename)
    info('')
    info('Regions (total ' + str(len(final_regions)) + ') saved into:')
    info('  ' + report.txt_fpath)

    return report


def get_detailed_metric_storage(depth_threshs):
    return MetricStorage(
        general_section=ReportSection(metrics=[
            Metric('Sample'),
        ]),
        sections=[ReportSection(metrics=[
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


def make_flat_region_report(sample, regions, depth_threshs):
    report = PerRegionSampleReport(sample=sample, metric_storage=get_detailed_metric_storage(depth_threshs))
    report.add_record('Sample', sample.name)

    i = 0
    for region in regions:
        i += 1
        if i % 10000 == 0:
            info('Processed {0:,} regions.'.format(i))

        rep_region = report.add_region()
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
        rep_region.add_record('W/n 20% of ave depth', region.rate_within_normal)

        for thresh in depth_threshs:
            if region.rates_within_threshs is None:
                err('Error: no rates_within_threshs for ' + str(region))
            else:
                # bases = region.bases_within_threshs.get(thresh)
                # if bases is not None and region.get_size() > 0:
                #     rate = 1.0 * bases / region.get_size()
                #     if rate > 1:
                #         critical(
                #             'Error: rate = ' + str(rate) + ', bases = ' + str(bases) +
                #             ', size = ' + str(region.get_size()) +
                #             ', start = ' + str(region.start) + ', end = ' + str(region.end))
                rep_region.add_record('{}x'.format(thresh), region.rates_within_threshs[thresh])

    info('Processed {0:,} regions.'.format(i))
    return report


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


def samtools_flag_stat(cnf, bam):
    samtools = get_system_path(cnf, 'samtools')
    output_fpath = join(cnf.work_dir, basename(bam) + '_flag_stats')
    cmdline = '{samtools} flagstat {bam}'.format(**locals())
    call(cnf, cmdline, output_fpath)
    stats = dict()
    with open(output_fpath) as f:
        lines = f.readlines()
        for stat, fun in [('total', number_of_reads),
                           ('duplicates', number_of_dup_reads),  # '-f 1024'
                           ('mapped', number_of_mapped_reads),   # '-F 4'
                           ('properly paired', number_of_properly_paired_reads)]:  # '-f 2'
            try:
                val = next(l.split()[0] for l in lines if stat in l)
            except StopIteration:
                warn('Cannot extract ' + stat + ' from flagstat output ' + output_fpath + '. Trying samtools view -c...')
                val = None
            else:
                try:
                    val = int(val)
                except ValueError:
                    warn('Cannot parse value ' + str(val) + ' from ' + stat + ' from flagstat output ' + output_fpath + '. Trying samtools view -c...')
                    val = None
            if val is not None:
                stats[stat] = val
            else:
                stats[stat] = fun(cnf, bam)

    return stats


def number_of_reads(cnf, bam, suf=''):
    samtools = get_system_path(cnf, 'samtools')
    output_fpath = join(cnf.work_dir, basename(bam) + '_' + suf + 'num_reads')
    cmdline = '{samtools} view -c {bam}'.format(**locals())
    call(cnf, cmdline, output_fpath)
    with open(output_fpath) as f:
        return int(f.read().strip())


def number_of_mapped_reads(cnf, bam, suf=''):
    samtools = get_system_path(cnf, 'samtools')
    output_fpath = join(cnf.work_dir, basename(bam) + '_' + suf + 'num_mapped_reads')
    cmdline = '{samtools} view -c -F 4 {bam}'.format(**locals())
    call(cnf, cmdline, output_fpath)
    with open(output_fpath) as f:
        return int(f.read().strip())


def number_of_properly_paired_reads(cnf, bam):
    samtools = get_system_path(cnf, 'samtools')
    output_fpath = join(cnf.work_dir, basename(bam) + '_num_paired_reads')
    cmdline = '{samtools} view -c -f 2 {bam}'.format(**locals())
    call(cnf, cmdline, output_fpath)
    with open(output_fpath) as f:
        return int(f.read().strip())


def number_of_dup_reads(cnf, bam):
    samtools = get_system_path(cnf, 'samtools')
    output_fpath = join(cnf.work_dir, basename(bam) + '_num_dup_reads')
    cmdline = '{samtools} view -c -f 1024 {bam}'.format(**locals())
    call(cnf, cmdline, output_fpath)
    with open(output_fpath) as f:
        return int(f.read().strip())


def number_of_dup_mapped_reads(cnf, bam):
    samtools = get_system_path(cnf, 'samtools')
    output_fpath = join(cnf.work_dir, basename(bam) + '_num_dup_unmapped_reads')
    cmdline = '{samtools} view -c -F 4 -f 1024 {bam}'.format(**locals())  # 1024 (dup) + 4 (unmpapped)
    call(cnf, cmdline, output_fpath)
    with open(output_fpath) as f:
        return int(f.read().strip())


def number_mapped_reads_on_target(cnf, bed, bam):
    samtools = get_system_path(cnf, 'samtools')
    output_fpath = join(cnf.work_dir, basename(bam) + '_' + basename(bed) + '_num_mapped_reads_target')
    cmdline = '{samtools} view -c -F 4 -L {bed} {bam}'.format(**locals())
    call(cnf, cmdline, output_fpath)
    with open(output_fpath) as f:
        return int(f.read().strip())


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


# def _bases_by_depth(depth_vals, depth_thresholds):
#     bases_by_min_depth = {depth: 0 for depth in depth_thresholds}
#
#     for depth_value in depth_vals:
#         for threshold in depth_thresholds:
#             if depth_value >= threshold:
#                 bases_by_min_depth[threshold] += 1
#
#         return [1.0 * bases_by_min_depth[thres] / len(depth_vals) if depth_vals else 0
#                 for thres in depth_thresholds]


