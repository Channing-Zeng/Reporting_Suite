# coding=utf-8

from collections import OrderedDict, defaultdict
from os.path import join, basename, isfile, abspath, realpath, splitext, normpath, dirname, relpath
import shutil
import traceback

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
from source.targetcov.Region import Region, save_regions_to_bed, GeneInfo, calc_bases_within_threshs, \
    calc_rate_within_normal, build_gene_objects_list
from source.targetcov.bam_and_bed_utils import index_bam, prepare_beds, filter_bed_with_gene_set, get_total_bed_size, \
    total_merge_bed, count_bed_cols, annotate_amplicons, group_and_merge_regions_by_gene, sort_bed, fix_bed_for_qualimap, \
    remove_dups, get_padded_bed_file, number_mapped_reads_on_target, samtools_flag_stat, calc_region_number
from source.targetcov.coverage_hist import bedcoverage_hist_stats
from source.tools_from_cnf import get_system_path, get_script_cmdline
from source.utils import get_chr_len_fpath


chr_y_bed = join(dirname(abspath(__file__)), 'chrY.bed')


def get_header_metric_storage(depth_thresholds, is_wgs=False):
    sections = [
        ReportSection('reads', 'Reads', [
            Metric('Reads'),
            Metric('Mapped reads', short_name='Mapped', description='samtools view -c -F 4', ok_threshold='Percentage of mapped reads', bottom=0),
            Metric('Percentage of mapped reads', short_name='Mapped %', unit='%', ok_threshold=0.98, bottom=0),
            Metric('Properly paired reads percent', short_name='Paired %', unit='%', description='Pecent of properly paired mapped reads (-f 2).', ok_threshold=0.9, bottom=0),
            Metric('Duplication rate', short_name='Dup rate', description='Percent of mapped reads (-F 4), marked as duplicates (-f 1024)', quality='Less is better', unit='%'),
            Metric('Read min length',  'Min len', 'Read min length'),
            Metric('Read max length',  'Max len', 'Read max length'),
            Metric('Read mean length', 'Ave len', 'Read mean length'),
            Metric('Gender', is_hidden=True),
        ]),
    ]
    if not is_wgs:
        sections.extend([
            ReportSection('target_metrics', 'Target coverage', [
                Metric('Covered bases in target', short_name='Trg covered', unit='bp'),
                Metric('Percentage of target covered by at least 1 read', short_name='%', unit='%'),
                Metric('Percentage of reads mapped on target', short_name='% reads on trg', unit='%', description='Percentage of reads overlapping target at least by 1 base'),
                Metric('Percentage of reads mapped off target', short_name='% reads off trg', unit='%', quality='Less is better', description='Percentage of reads that don\'t overlap target even by 1 base'),
                Metric('Percentage of reads mapped on padded target', short_name='% reads on padded trg', unit='%', description='Percentage of reads that overlap target at least by 1 base. Should be 1-2% higher.'),
                Metric('Read bases mapped on target', short_name='Read bp on trg', unit='bp'),
            ]),
        ])
    else:
        sections.extend([
            ReportSection('target_metrics_wgs', 'Genome coverage', [
                Metric('Covered bases in genome', short_name='Genome covered', unit='bp'),
                Metric('Percentage of genome covered by at least 1 read', short_name='Genome cvrd %', unit='%'),
                Metric('Covered bases in exome', short_name='Exome covered', unit='bp'),
                Metric('Percentage of exome covered by at least 1 read', short_name='%', unit='%'),
                Metric('Percentage of reads mapped on exome', short_name='% reads on exome', unit='%'),
                Metric('Percentage of reads mapped off exome', short_name='% off exome', unit='%', quality='Less is better'),
            ]),
        ])

    trg_name = 'target' if not is_wgs else 'genome'
    depth_section = ReportSection('depth_metrics', ('Target' if not is_wgs else 'Genome') + ' coverage depth', [
        Metric('Average ' + trg_name + ' coverage depth', short_name='Avg'),
        Metric('Std. dev. of ' + trg_name + ' coverage depth', short_name='Std dev', quality='Less is better'),
        # Metric('Minimal ' + trg_name + ' coverage depth', short_name='Min', is_hidden=True),
        # Metric('Maximum ' + trg_name + ' coverage depth', short_name='Max', is_hidden=True),
        Metric('Percentage of ' + trg_name + ' within 20% of mean depth', short_name='&#177;20% avg', unit='%')
    ])
    for depth in depth_thresholds:
        name = 'Part of ' + trg_name + ' covered at least by ' + str(depth) + 'x'
        depth_section.add_metric(Metric(name, short_name=str(depth) + 'x', description=name, unit='%'))
    sections.append(depth_section)

    sections.append(
        ReportSection('qualimap', 'Qualimap stats' + ('' if is_wgs else ' within the target'), [
            Metric('Mean Mapping Quality',  'Mean MQ',            'Mean mapping quality, inside of regions'),
            Metric('Mismatches',            'Mismatches',         'Mismatches, inside of regions', quality='Less is better'),  # added in Qualimap v.2.0
            Metric('Insertions',            'Insertions',         'Insertions, inside of regions', quality='Less is better'),
            Metric('Deletions',             'Deletions',          'Deletions, inside of regions', quality='Less is better'),
            Metric('Homopolymer indels',    'Homopolymer indels', 'Percentage of homopolymer indels, inside of regions', quality='Less is better'),
            Metric('Qualimap',              'Qualimap',           'Qualimap report'),
            Metric('ngsCAT',                'ngsCAT',             'ngsCAT report')
        ])
    )

    ms = MetricStorage(
        general_section=ReportSection('general_section', '', [
            Metric('Target', short_name='Target', common=True),
            Metric('Reference size', short_name='Reference bp', unit='bp', common=True),
            # Metric('Target ready', short_name='Target ready', common=True),
            Metric('Regions in target', short_name='Regions in target', common=True),
            Metric('Bases in target', short_name='Target bp', unit='bp', common=True),
            Metric('Percentage of reference', short_name='Percentage of reference', unit='%', common=True),
            # Metric('Genes', short_name='Genes', common=True),
            Metric('Genes in target', short_name='Genes in target', common=True),
        ]),
        sections=sections
    )
    return ms


class TargetInfo:
    def __init__(self, fpath=None, bed=None, original_target_bed=None, regions_num=None,
                 bases_num=None, fraction=None, genes_fpath=None, genes_num=None):
        self.fpath = realpath(fpath) if fpath else None  # raw source file - to demonstrate where we took one
        self.bed = bed                                   # processed (sorted, merged...), to do real calculations
        self.original_target_bed = original_target_bed
        self.regions_num = regions_num
        self.bases_num = bases_num
        self.fraction = fraction
        self.genes_fpath = realpath(genes_fpath) if genes_fpath else None
        self.genes_num = genes_num
        self.type = None  # 'Regional', 'WGS'


def _run_qualimap(cnf, sample, bam_fpath, bed_fpath=None):
    safe_mkdir(dirname(sample.qualimap_dirpath))  # do not create the directory itself!!! it will be check for reuse

    bed = ''
    if bed_fpath:
        qualimap_bed_fpath = join(cnf.work_dir, 'tmp_qualimap.bed')
        fix_bed_for_qualimap(bed_fpath, qualimap_bed_fpath)
        bed = '--bed ' + qualimap_bed_fpath

    qm = get_system_path(cnf, 'python', join('scripts', 'post', 'qualimap.py'))
    cmdl = '{qm} --project-name {cnf.project_name} --sys-cnf {cnf.sys_cnf} --run-cnf {cnf.run_cnf} ' \
           '--bam {bam_fpath} {bed} -o {sample.qualimap_dirpath} -t {cnf.threads}'.format(**locals())
    call(cnf, cmdl, sample.qualimap_html_fpath, stdout_to_outputfile=False)
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

    # dedup_bam_dirpath = join(cnf.work_dir, source.dedup_bam)
    # safe_mkdir(dedup_bam_dirpath)
    # dedup_bam_fpath = join(dedup_bam_dirpath, add_suffix(basename(bam_fpath), source.dedup_bam))
    # remove_dups(cnf, bam_fpath, output_fpath=dedup_bam_fpath, use_grid=False)
    # dedup_bam_stats = samtools_flag_stat(cnf, dedup_bam_fpath)
    # info('Total reads after dedup (samtools view -F 1024): ' + Metric.format_value(dedup_bam_stats['total']))
    # info('Total mapped reads after dedup (samtools view -F 1024): ' + Metric.format_value(dedup_bam_stats['mapped']))
    return bam_stats  # dedup_bam_fpath, bam_stats, dedup_bam_stats


def _parse_qualimap_results(qualimap_html_fpath, qualimap_cov_hist_fpath, depth_thresholds):
    if not verify_file(qualimap_html_fpath):
        critical('Qualimap report was not found')

    depth_stats = dict()

    if not verify_file(qualimap_cov_hist_fpath):
        err('Qualimap hist fpath is not found, cannot build histogram')
    else:
        bases_by_depth = OrderedDict()
        with open(qualimap_cov_hist_fpath) as f:
            for l in f:
                if l.startswith('#'):
                    pass
                else:
                    cov, bases = map(int, map(float, l.strip().split()))
                    bases_by_depth[cov] = bases

        depth_stats['bases_by_depth'] = bases_by_depth

    qualimap_records = parse_qualimap_sample_report(qualimap_html_fpath)

    def find_rec(name):
        return next((r.value for r in qualimap_records if r.metric.name.startswith(name)), None)

    depth_stats['ave_depth'] = find_rec('Coverage Mean')
    # depth_stats['min_depth'] = find_rec('Coverage Min')
    # depth_stats['max_depth'] = find_rec('Coverage Max')
    depth_stats['stddev_depth'] = find_rec('Coverage Standard Deviation')

    target_stats = dict(
        reference_size  = find_rec('Reference size'),
        target_size     = find_rec('Regions size/percentage of reference (on target)'),
        target_fraction = find_rec('Regions size/percentage of reference (on target) %'),
    )

    reads_stats = dict(
        total                    = find_rec('Number of reads'),
        mapped                   = find_rec('Mapped reads'),
        mapped_rate              = find_rec('Mapped reads %'),
        unmapped                 = find_rec('Unmapped reads'),
        unmapped_rate            = find_rec('Unmapped reads %'),
        # mapped_on_target         = find_rec('Mapped reads (on target)'),
        # mapped_rate_on_target    = find_rec('Mapped reads % (on target)'),
        paired                   = find_rec('Paired reads'),
        paired_rate              = find_rec('Paired reads %'),
        dup                      = find_rec('Duplicated reads (flagged)'),
        dup_rate                 = find_rec('Duplicated reads (flagged) %'),
        min_len                  = find_rec('Read min length'),
        max_len                  = find_rec('Read max length'),
        ave_len                  = find_rec('Read mean length'),
    )

    mm_indels_stats = dict(
        mean_mq     = find_rec('Mean Mapping Quality'),
        mismatches  = find_rec('Mismatches'),
        insertions  = find_rec('Insertions'),
        deletions   = find_rec('Deletions'),
        homo_indels = find_rec('Homopolymer indels'),
    )

    return depth_stats, reads_stats, mm_indels_stats, target_stats


def _get_gender(sample, bam_fpath, cnf):
    thres = 100
    info('Detecting gender by chrY key regions coverage. The read number threshold mapped on the regions is ' + str(thres))
    if not bam_fpath:
        critical(sample.name + ': BAM file is required.')
    if not isfile(bam_fpath + '.bai'):
        info('Indexing bam ' + bam_fpath)
        index_bam(cnf, bam_fpath)

    reads_mapped_on_chr_y = number_mapped_reads_on_target(cnf, chr_y_bed, bam_fpath)
    info('Number of reads mapped on chrY key genes is ' + str(reads_mapped_on_chr_y))
    gender = 'F'
    if reads_mapped_on_chr_y > thres:
        gender = 'M'
    info('Gender is ' + gender)
    return gender


def make_targetseq_reports(cnf, output_dir, sample, bam_fpath, exons_bed, exons_no_genes_bed, target_bed, gene_names_list):
    info('Starting targeqSeq for ' + sample.name + ', saving into ' + output_dir)
    gene_by_name = build_gene_objects_list(cnf, sample.name, exons_bed, gene_names_list)

    # ref_fapth = cnf.genome.seq
    original_target_bed = cnf.original_target_bed or target_bed
    target_info = TargetInfo(
        fpath=target_bed, bed=target_bed, original_target_bed=original_target_bed,
        genes_num=len(gene_by_name) if gene_by_name else None)
    if target_bed:
        target_info.regions_num = calc_region_number(target_bed)

    _run_qualimap(cnf, sample, bam_fpath, target_bed)

    depth_stats, reads_stats, mm_indels_stats, target_stats = _parse_qualimap_results(
        sample.qualimap_html_fpath, sample.qualimap_cov_hist_fpath, cnf.coverage_reports.depth_thresholds)
    reads_stats['gender'] = _get_gender(sample, cnf.bam, cnf)

    if 'bases_by_depth' in depth_stats:
        depth_stats['bases_within_threshs'], depth_stats['rates_within_threshs'] = calc_bases_within_threshs(
            depth_stats['bases_by_depth'],
            target_stats['target_size'] or target_stats['reference_size'],
            cnf.coverage_reports.depth_thresholds)

        depth_stats['wn_20_percent'] = calc_rate_within_normal(
            depth_stats['bases_by_depth'],
            depth_stats['ave_depth'],
            target_stats['target_size'] or target_stats['reference_size'])

    if target_stats['target_size']:
        target_info.bases_num = target_stats['target_size']
        target_info.fraction  = target_stats['target_fraction']
    else:
        target_info.bases_num = target_stats['reference_size']

    sample.dedup_bam = add_suffix(bam_fpath, source.dedup_bam)
    remove_dups(cnf, bam_fpath, sample.dedup_bam)
    bam_fpath = sample.dedup_bam

    if target_info.bed:
        reads_stats['mapped_on_target'] = number_mapped_reads_on_target(cnf, target_bed, bam_fpath)

    if target_info.bed:
        padded_bed = get_padded_bed_file(cnf, target_info.bed, get_chr_len_fpath(cnf), cnf.coverage_reports.padding)
        reads_stats['mapped_on_padded_target'] = number_mapped_reads_on_target(cnf, padded_bed, bam_fpath)
    elif exons_no_genes_bed:
        reads_stats['mapped_on_exome'] = number_mapped_reads_on_target(cnf, exons_no_genes_bed, bam_fpath)

    summary_report = make_summary_report(cnf, depth_stats, reads_stats, mm_indels_stats, sample, output_dir, target_info)

    per_gene_report = make_per_gene_report(cnf, sample, bam_fpath, target_bed, exons_bed,
                                           exons_no_genes_bed, output_dir, gene_by_name)

    # key_genes_report = make_key_genes_reports(cnf, sample, gene_by_name, depth_stats['ave_depth'])

    info()
    return depth_stats['ave_depth'], gene_by_name, [summary_report, per_gene_report]


def get_records_by_metrics(records, metrics):
    _records = []
    for rec in records:
        if rec.metric.name in metrics:
            rec.metric = metrics[rec.metric.name]
            _records.append(rec)
    return _records


def make_summary_report(cnf, depth_stats, reads_stats, mm_indels_stats, sample, output_dir, target_info):
    report = SampleReport(sample, metric_storage=get_header_metric_storage(cnf.coverage_reports.depth_thresholds, is_wgs=target_info.bed is None))
    report.add_record('Qualimap', value='Qualimap', html_fpath=relpath(sample.qualimap_html_fpath, output_dir), silent=True)
    report.add_record('Gender', reads_stats['gender'], silent=True)

    info('* General coverage statistics *')
    report.add_record('Reads', reads_stats['total'])
    report.add_record('Mapped reads', reads_stats['mapped'])
    # report.add_record('Unmapped reads', reads_stats['total'] - reads_stats['mapped'])
    percent_mapped = 1.0 * (reads_stats['mapped'] or 0) / reads_stats['total'] if reads_stats['total'] else None
    assert percent_mapped <= 1.0 or percent_mapped is None, str(percent_mapped)
    report.add_record('Percentage of mapped reads', percent_mapped)
    # percent_unmapped = 1.0 * (reads_stats['total'] - reads_stats['mapped']) / reads_stats['total'] if reads_stats['total'] else None
    # assert percent_unmapped <= 1.0 or percent_unmapped is None, str(percent_unmapped)
    # report.add_record('Percentage of unmapped reads', percent_unmapped)
    total_paired_reads_pecent = 1.0 * (reads_stats['paired'] or 0) / reads_stats['total'] if reads_stats['total'] else None
    assert total_paired_reads_pecent <= 1.0 or total_paired_reads_pecent is None, str(total_paired_reads_pecent)
    report.add_record('Properly paired reads percent', total_paired_reads_pecent)
    info('')

    # if dedup_bam_stats:
    # dup_rate = 1 - (1.0 * dedup_bam_stats['mapped'] / bam_stats['mapped']) if bam_stats['mapped'] else None
    report.add_record('Duplication rate', reads_stats['dup_rate'])
    # report.add_record('Dedupped mapped reads', reads_stats['mapped'] - reads_stats[''])

    if target_info.bed:
        info('* Target coverage statistics *')
        if target_info.original_target_bed:
            report.add_record('Target', target_info.original_target_bed)
            # report.add_record('Ready target', target_info.fpath)
        else:
            report.add_record('Target', target_info.fpath)
        report.add_record('Bases in target', target_info.bases_num)
        report.add_record('Percentage of reference', target_info.fraction)
        report.add_record('Regions in target', target_info.regions_num)
    else:
        info('* Genome coverage statistics *')
        report.add_record('Target', 'whole genome')
        report.add_record('Reference size', target_info.bases_num)

    report.add_record('Genes in target', target_info.genes_num)

    trg_type = 'target' if target_info.bed else 'genome'

    if 'bases_within_threshs' in depth_stats:
        bases_within_threshs = depth_stats['bases_within_threshs']
        v_covered_bases_in_targ = bases_within_threshs.items()[0][1]
        v_percent_covered_bases_in_targ = 1.0 * (v_covered_bases_in_targ or 0) / target_info.bases_num if target_info.bases_num else None
        assert v_percent_covered_bases_in_targ <= 1.0 or v_percent_covered_bases_in_targ is None, str(v_percent_covered_bases_in_targ)

        report.add_record('Covered bases in ' + trg_type, v_covered_bases_in_targ)
        report.add_record('Percentage of ' + trg_type + ' covered by at least 1 read', v_percent_covered_bases_in_targ)

    if target_info.bed:
        info('Getting number of mapped reads on target...')
        # mapped_reads_on_target = number_mapped_reads_on_target(cnf, target_info.bed, bam_fpath)
        if 'mapped_on_target' in reads_stats:
            # report.add_record('Reads mapped on target', reads_stats['mapped_on_target'])
            percent_mapped_on_target = 1.0 * (reads_stats['mapped_on_target'] or 0) / reads_stats['mapped'] if reads_stats['mapped'] != 0 else None
            report.add_record('Percentage of reads mapped on target', percent_mapped_on_target)
            assert percent_mapped_on_target <= 1.0 or percent_mapped_on_target is None, str(percent_mapped_on_target)
            percent_mapped_off_target = 1.0 - percent_mapped_on_target
            report.add_record('Percentage of reads mapped off target ', percent_mapped_off_target)

        read_bases_on_targ = int(target_info.bases_num * depth_stats['ave_depth'])  # sum of all coverages
        report.add_record('Read bases mapped on target', read_bases_on_targ)

        if 'mapped_on_padded_target' in reads_stats:
            # report.add_record('Reads mapped on padded target', reads_stats['mapped_reads_on_padded_target'])
            percent_mapped_on_padded_target = 1.0 * (reads_stats['mapped_on_padded_target'] or 0) / reads_stats['mapped'] if reads_stats['mapped'] else None
            report.add_record('Percentage of reads mapped on padded target', percent_mapped_on_padded_target)
            assert percent_mapped_on_padded_target <= 1.0 or percent_mapped_on_padded_target is None, str(percent_mapped_on_padded_target)

    elif 'mapped_on_exome' in reads_stats:
        # report.add_record('Reads mapped on target', reads_stats['mapped_on_target'])
        percent_mapped_on_exome = 1.0 * (reads_stats['mapped_on_exome'] or 0) / reads_stats['mapped'] if reads_stats['mapped'] != 0 else None
        if percent_mapped_on_exome:
            report.add_record('Percentage of reads mapped on exome', percent_mapped_on_exome)
            assert percent_mapped_on_exome <= 1.0 or percent_mapped_on_exome is None, str(percent_mapped_on_exome)
            percent_mapped_off_exome = 1.0 - percent_mapped_on_exome
            report.add_record('Percentage of reads mapped off exome ', percent_mapped_off_exome)

    info('')
    report.add_record('Average ' + trg_type + ' coverage depth', depth_stats['ave_depth'])
    report.add_record('Std. dev. of ' + trg_type + ' coverage depth', depth_stats['stddev_depth'])
    # report.add_record('Minimal ' + trg_type + ' coverage depth', depth_stats['min_depth'])
    # report.add_record('Maximum ' + trg_type + ' coverage depth', depth_stats['max_depth'])
    if 'wn_20_percent' in depth_stats:
        report.add_record('Percentage of ' + trg_type + ' within 20% of mean depth', depth_stats['wn_20_percent'])
        assert depth_stats['wn_20_percent'] <= 1.0 or depth_stats['wn_20_percent'] is None, str( depth_stats['wn_20_percent'])

    if 'bases_within_threshs' in depth_stats:
        for depth, bases in depth_stats['bases_within_threshs'].items():
            percent_val = 1.0 * (bases or 0) / target_info.bases_num if target_info.bases_num else 0
            if percent_val > 0:
                report.add_record('Part of ' + trg_type + ' covered at least by ' + str(depth) + 'x', percent_val)
            assert percent_val <= 1.0 or percent_val is None, str(percent_val)
    info()

    report.add_record('Read mean length', reads_stats['ave_len'])
    report.add_record('Read min length', reads_stats['min_len'])
    report.add_record('Read max length', reads_stats['max_len'])
    report.add_record('Mean Mapping Quality', mm_indels_stats['mean_mq'])
    report.add_record('Mismatches', mm_indels_stats['mismatches'])
    report.add_record('Insertions', mm_indels_stats['insertions'])
    report.add_record('Deletions', mm_indels_stats['deletions'])
    report.add_record('Homopolymer indels', mm_indels_stats['homo_indels'])

    report.save_json(join(output_dir, sample.name + '.' + source.targetseq_name + '.json'))
    report.save_txt (join(output_dir, sample.name + '.' + source.targetseq_name + '.txt'))
    report.save_html(join(output_dir, sample.name + '.' + source.targetseq_name + '.html'), caption='Target coverage statistics for ' + sample.name)
    info()
    info('Saved to ')
    info('  ' + report.txt_fpath)
    return report


def make_per_gene_report(cnf, sample, bam_fpath, target_bed, exons_bed, exons_no_genes_bed, output_dir, gene_by_name):
    targets_or_exons = []

    if exons_no_genes_bed or target_bed:
        ready_target_bed = join(output_dir, 'target.bed')
        try:
            shutil.copy(target_bed or exons_bed, ready_target_bed)
        except OSError:
            err(traceback.format_exc())

        info()
        info('Calculation of coverage statistics for the regions in the input BED file...')
        targets_or_exons = bedcoverage_hist_stats(cnf, sample.name, bam_fpath, target_bed or exons_no_genes_bed)
        info()
        # info('Merging capture BED file to get total target cov statistics...')
        # total_merged_target_bed = total_merge_bed(cnf, target_bed or exons_no_genes_bed)
        # info()
        # info('Calculation of coverage statistics for total target...')
        # _, combined_region, _ = bedcoverage_hist_stats(cnf, sample.name, bam_fpath, total_merged_target_bed)

        # total_bed_size = get_total_bed_size(cnf, target_bed or exons_no_genes_bed)

    un_annotated_amplicons = []

    if exons_bed:
        if not target_bed:
            exons_with_optional_genes = targets_or_exons
        else:
            info()
            info('Calculating coverage statistics for exons...')
            exons_with_optional_genes = bedcoverage_hist_stats(cnf, sample.name, bam_fpath, exons_no_genes_bed)

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
            for ampl in targets_or_exons:
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
        per_gene_report = _generate_report_from_regions(cnf, sample, output_dir,
            gene_by_name.values(), un_annotated_amplicons)

    return per_gene_report


def _generate_report_from_regions(cnf, sample, output_dir, genes, un_annotated_amplicons):
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
    report.save_txt(sample.targetcov_detailed_txt)
    report.save_tsv(sample.targetcov_detailed_tsv)
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
        ])]
    )


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


