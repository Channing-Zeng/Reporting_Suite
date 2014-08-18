# coding=utf-8
from collections import OrderedDict
import copy
from itertools import izip, chain, repeat
from os.path import join, basename
from source.bcbio_structure import BCBioStructure, Sample
from source.calling_process import call, call_check_output, call_pipe
from source.file_utils import intermediate_fname, splitext_plus

from source.logger import step_greetings, critical, info, err
from source.reporting import Metric, Record, write_txt_reports, save_json, SampleReport
from source.targetcov.Region import Region
from source.tools_from_cnf import get_tool_cmdline
from source.utils import get_chr_len_fpath
from source.file_utils import file_transaction


def run_target_cov(cnf, bam, amplicons_bed):
    summary_report_fpath = None
    gene_report_fpath = None

    info()
    info('Calculation of coverage statistics for the regions in the input BED file...')
    amplicons, combined_region, max_depth, total_bed_size = _bedcoverage_hist_stats(cnf, bam, amplicons_bed)

    chr_len_fpath = get_chr_len_fpath(cnf)
    exons_bed = cnf['genome']['exons']

    if 'summary' in cnf['coverage_reports']['report_types']:
        step_greetings('Target coverage summary report')

        records = _run_header_report(
            cnf, amplicons_bed, bam, chr_len_fpath,
            cnf['coverage_reports']['depth_thresholds'], cnf['padding'],
            combined_region, max_depth, total_bed_size)

        save_json(records, join(cnf.output_dir, cnf.name + '.' + BCBioStructure.targetseq_name + '.json'))
        summary_report_fpath = write_txt_reports(
            cnf.output_dir, cnf.work_dir,
            [SampleReport(Sample(cnf.name, bam, amplicons_bed), '', records)],
            cnf.name + '.' + BCBioStructure.targetseq_name)
        info()
        info('Saved to ' + summary_report_fpath)


    if 'genes' in cnf['coverage_reports']['report_types']:
        if not exons_bed:
            if cnf['reports'] == 'genes':
                critical('Error: no genes or exons specified for the genome in system config, '
                         'cannot run per-exon report.')
            else:
                err('Warning: no genes or exons specified for the genome in system config, '
                    'cannot run per-exon report.')
        else:
            # log('Annotating amplicons.')
            # annotate_amplicons(amplicons, genes_bed)

            info('Sorting exons BED file.')
            exons_bed = sort_bed(cnf, exons_bed)
            info('Getting the exons that overlap amplicons.')
            overlapped_exons_bed = intersect_bed(cnf, exons_bed, amplicons_bed)
            info('Adding other exons for the genes of overlapped exons.')
            target_exons_bed = _add_other_exons(cnf, exons_bed, overlapped_exons_bed)

            info('Calculation of coverage statistics for exons of the genes ovelapping with the input regions...')
            exons, _, _, _ = _bedcoverage_hist_stats(cnf, bam, target_exons_bed)
            for exon in exons:
                exon.gene_name = exon.extra_fields[0]

            gene_report_fpath = join(cnf.output_dir,
                                     cnf.name + '.' + BCBioStructure.targetseq_dir +
                                     BCBioStructure.detail_gene_report_ending)
            info('Region cov report...')
            _run_region_cov_report(cnf, gene_report_fpath, cnf.name, cnf.coverage_reports.depth_thresholds,
                                   amplicons, exons)

    return summary_report_fpath, gene_report_fpath


def _add_other_exons(cnf, exons_bed, overlapped_exons_bed):
    gene_names = set()
    with open(overlapped_exons_bed) as overl_f:
        for line in overl_f:
            if not line or not line.strip() or line.startswith('#') or len(line.split()) < 4:
                continue
            gene_names.add(line.split()[3])

    new_overlp_exons_bed = intermediate_fname(cnf, overlapped_exons_bed, 'by_genes')
    with open(exons_bed) as exons_f, \
         open(new_overlp_exons_bed, 'w') as new_overl_f:
        for line in exons_f:
            if not line or not line.strip() or line.startswith('#') or len(line.split()) < 4:
                new_overl_f.write(line)
            elif line.split()[3] in gene_names:
                new_overl_f.write(line)

    return sort_bed(cnf, new_overlp_exons_bed)


_basic_metrics = Metric.to_dict([
    Metric('Reads'),
    Metric('Mapped reads', short_name='Mapped'),
    Metric('Percentage of mapped reads', short_name='%', unit='%'),

    Metric('Unmapped reads', short_name='Unmapped'),
    Metric('Percentage of unmapped reads', short_name='%'),

    Metric('Bases in target', short_name='Target bp'),

    Metric('Covered bases in target', short_name='Covered'),
    Metric('Percentage of target covered by at least 1 read', short_name='%', unit='%'),

    Metric('Reads mapped on target', short_name='Reads on trg'),
    Metric('Percentage of reads mapped on target', short_name='%', unit='%'),

    Metric('Reads mapped on padded target', 'On padded trg'),
    Metric('Percentage of reads mapped on padded target', short_name='%', unit='%'),

    Metric('Read bases mapped on target', short_name='Read bp on trg')])


_depth_metrics = Metric.to_dict([
    Metric('Average target coverage depth', short_name='Avg'),
    Metric('Std. dev. of target coverage depth', short_name='Std dev'),
    Metric('Maximum target coverage depth', short_name='Max'),
    Metric('Percentage of target within 20% of mean depth', short_name='&#177;20% avg', unit='%')])


def _add_depth_metrics(depth_thresholds):
    for depth in depth_thresholds:
        name = 'Part of target covered at least by ' + str(depth) + 'x'
        _depth_metrics[name] = Metric(name, short_name=str(depth) + 'x', description=name, unit='%')

def get_basic_metrics():
    return _basic_metrics

def get_depth_metrics(depth_thresholds):
    _add_depth_metrics(depth_thresholds)
    return _depth_metrics

def get_cov_metrics(depth_thresholds):
    _add_depth_metrics(depth_thresholds)
    return OrderedDict(_basic_metrics.items() + _depth_metrics.items())


def _run_header_report(cnf, bed, bam, chr_len_fpath,
                       depth_thresholds, padding,
                       combined_region, max_depth, total_bed_size):
    records = []
    cov_metrics = get_cov_metrics(depth_thresholds)

    def append_stat(metric_name, value):
        rec = Record(cov_metrics[metric_name.strip()], value)
        records.append(rec)
        info(metric_name + ': ' + rec.format())

    info('* General coverage statistics *')
    info('Getting number of reads...')
    v_number_of_reads = number_of_reads(cnf, bam)
    append_stat('Reads', v_number_of_reads)

    info('Getting number of mapped reads...')
    v_mapped_reads = number_of_mapped_reads(cnf, bam)
    v_percent_mapped = 100.0 * v_mapped_reads / v_number_of_reads if v_number_of_reads else None
    v_percent_unmapped = 100.0 * (v_number_of_reads - v_mapped_reads) / v_number_of_reads if v_number_of_reads else None
    append_stat('Mapped reads', v_mapped_reads)
    append_stat('Percentage of mapped reads', v_percent_mapped)
    append_stat('Unmapped reads', v_number_of_reads - v_mapped_reads)
    append_stat('Percentage of unmapped reads', v_percent_unmapped)
    info('')

    info('* Target coverage statistics *')
    append_stat('Bases in target', total_bed_size)

    bases_within_threshs, avg_depth, std_dev, percent_within_normal = combined_region.sum_up(depth_thresholds)

    v_covered_bases_in_targ = bases_within_threshs.items()[0][1]
    append_stat('Covered bases in target', v_covered_bases_in_targ)

    v_percent_covered_bases_in_targ = 100.0 * v_covered_bases_in_targ / total_bed_size if total_bed_size else None
    append_stat('Percentage of target covered by at least 1 read', v_percent_covered_bases_in_targ)
    info('Getting number of mapped reads on target...')

    v_mapped_reads_on_target = number_mapped_reads_on_target(cnf, bed, bam)
    append_stat('Reads mapped on target', v_mapped_reads_on_target)

    v_percent_mapped_on_target = 100.0 * v_mapped_reads_on_target / v_mapped_reads if v_mapped_reads else None
    append_stat('Percentage of reads mapped on target ', v_percent_mapped_on_target)

    info('Making bed file for padded regions...')
    padded_bed = get_padded_bed_file(cnf, bed, chr_len_fpath, padding)

    info('Getting number of mapped reads on padded target...')
    v_reads_on_padded_targ = number_mapped_reads_on_target(cnf, padded_bed, bam)
    append_stat('Reads mapped on padded target', v_reads_on_padded_targ)

    v_percent_mapped_on_padded_target = 100.0 * v_reads_on_padded_targ / v_mapped_reads if v_mapped_reads else None
    append_stat('Percentage of reads mapped on padded target', v_percent_mapped_on_padded_target)

    v_read_bases_on_targ = avg_depth * total_bed_size  # sum of all coverages
    append_stat('Read bases mapped on target', v_read_bases_on_targ)

    info('')
    append_stat('Average target coverage depth', avg_depth)
    append_stat('Std. dev. of target coverage depth', std_dev)
    append_stat('Maximum target coverage depth', max_depth)
    append_stat('Percentage of target within 20% of mean depth', percent_within_normal)

    for depth, bases in bases_within_threshs.items():
        percent = 100.0 * bases / total_bed_size if total_bed_size else 0
        append_stat('Part of target covered at least by ' + str(depth) + 'x', percent)

    return records


def _run_region_cov_report(cnf, report_fpath, sample_name, depth_threshs,
                          amplicons, exons):
    for ampl in amplicons:
        ampl.feature = 'Amplicon'
        ampl.sample = sample_name

    for exon in exons:
        exon.feature = 'Exon'
        exon.sample = sample_name

    info('Groupping exons per gene...')
    exon_genes = _get_exon_genes(cnf, exons)
    info()

    info('Groupping amplicons per gene...')
    amplicon_genes_by_name = _get_amplicon_genes(amplicons, exon_genes)
    info()

    result_regions = []
    info('Combining...')
    i = 0
    for exon_gene in exon_genes:
        if i and i % 10000 == 0:
            info('Processed {0:,} genes, current gene {1}'.format(i, exon_gene.gene_name))
        i += 1

        for exon in exon_gene.subregions:
            result_regions.append(exon)
        result_regions.append(exon_gene)

        amplicon_gene = amplicon_genes_by_name.get(exon_gene.gene_name)
        if amplicon_gene:
            for amplicon in amplicon_gene.subregions:
                result_regions.append(amplicon)
            result_regions.append(amplicon_gene)
    info('Processed {0:,} genes.'.format(i))

    info()
    info('Building final regions report...')
    return _build_regions_cov_report(
        cnf, report_fpath, depth_threshs, result_regions)


def _get_amplicon_genes(amplicons, exon_genes):
    amplicon_genes_by_name = dict()

    i = 0
    for exon_gene in exon_genes:
        if i and i % 1000 == 0:
            info('Processed {0:,} regions, current gene {1}'.format(i, exon_gene.gene_name))
        i += 1

        for amplicon in amplicons:
            if exon_gene.intersect(amplicon):
                amplicon_gene = amplicon_genes_by_name.get(exon_gene.gene_name)
                if amplicon_gene is None:
                    amplicon_gene = Region(
                        sample=amplicon.sample, chrom=amplicon.chrom,
                        start=amplicon.start, end=amplicon.end, size=0,
                        feature='Gene-' + amplicon.feature,
                        gene_name=exon_gene.gene_name)
                    amplicon_genes_by_name[amplicon_gene.gene_name] = amplicon_gene
                amplicon_copy = copy.copy(amplicon)
                amplicon_gene.add_subregion(amplicon_copy)
                amplicon_copy.gene_name = amplicon_gene.gene_name

    info('Processed {0:,} regions.'.format(i))

    return amplicon_genes_by_name


def _get_exon_genes(cnf, subregions):
    genes_by_name = dict()

    i = 0
    for exon in subregions:
        if not exon.gene_name:
            info()
            err('No gene name info in the record: ' + str(exon) + '. Skipping.')
            continue

        if i and i % 10000 == 0:
            info('Processed {0:,} exons, current gene {1}'.format(i, exon.gene_name))
        i += 1

        gene = genes_by_name.get(exon.gene_name)
        if gene is None:
            gene = Region(
                sample=exon.sample, chrom=exon.chrom,
                start=exon.start, end=exon.end, size=0,
                feature='Gene-' + exon.feature,
                gene_name=exon.gene_name)
            genes_by_name[exon.gene_name] = gene
        gene.add_subregion(exon)

    sorted_genes = sorted(genes_by_name.values(), key=lambda g: (g.chrom, g.start, g.end))

    info('Processed {0:,} exons.'.format(i))

    return sorted_genes


def _build_regions_cov_report(cnf, report_fpath, depth_threshs, regions,
                             extra_headers=list()):
    header_fields = ['SAMPLE', 'Chr', 'Start', 'End', 'Gene', 'Feature', 'Size',
                     'Mean Depth', 'Standard Dev.', 'Within 20% of Mean'] +\
                    ['{}x'.format(thres) for thres in depth_threshs]
    max_lengths = map(len, header_fields)

    all_values = []

    i = 0
    for region in regions:
        i += 1
        if i % 10000 == 0:
            info('Processed {0:,} regions.'.format(i))

        bases_within_threshs, avg_depth, std_dev, percent_within_normal = region.sum_up(depth_threshs)

        line_fields = map(str, [region.sample,
                                region.chrom,
                                region.start,
                                region.end,
                                region.gene_name])
        line_fields += [region.feature]
        line_fields += [str(region.get_size())]
        line_fields += ['{0:.2f}'.format(avg_depth)]
        line_fields += ['{0:.2f}'.format(std_dev)]
        line_fields += ['{0:.2f}%'.format(percent_within_normal)
                        if percent_within_normal is not None else '-']

        for depth_thres, bases in bases_within_threshs.items():
            if int(region.get_size()) == 0:
                percent_str = '-'
            else:
                percent = 100.0 * bases / region.get_size()
                percent_str = '{0:.2f}%'.format(percent)
            line_fields.append(percent_str)

        all_values.append(line_fields)
        max_lengths = map(max, izip(max_lengths, chain(map(len, line_fields), repeat(0))))
    info('Processed {0:,} regions.'.format(i))

    with file_transaction(cnf, report_fpath) as tx:
        with open(tx, 'w') as out, \
                open(join(cnf['work_dir'], basename(report_fpath)), 'w') as nice_out:
            out.write('\t'.join(header_fields) + '\n')

            for line_fields in all_values:
                out.write('\t'.join(line_fields) + '\n')

            for h, l in zip(header_fields, max_lengths):
                nice_out.write(h + ' ' * (l - len(h) + 2))
            nice_out.write('\n')

            for line_tokens in all_values:
                for v, l in zip(line_tokens, max_lengths):
                    nice_out.write(v + ' ' * (l - len(v) + 2))
                nice_out.write('\n')
    info('')
    info('Result: ' + report_fpath)
    return report_fpath


def _bedcoverage_hist_stats(cnf, bam, bed):
    regions, max_depth, total_bed_size = [], 0, 0

    bedtools = get_tool_cmdline(cnf, 'bedtools')
    cmdline = '{bedtools} coverage -abam {bam} -b {bed} -hist'.format(**locals())
    bedcov_output = join(cnf.work_dir, splitext_plus(basename(bed))[0] + '_' +
                                       splitext_plus(basename(bam))[0] + '_bedcov_output.txt')
    call(cnf, cmdline, bedcov_output)

    _total_regions_count = 0

    with open(bedcov_output) as f:
        for next_line in f:
            if not next_line.strip() or next_line.startswith('#'):
                continue

            line_tokens = next_line.strip().split()
            chrom = line_tokens[0]
            start, end = None, None
            try:
                depth, bases, region_size = map(int, line_tokens[-4:-1])
            except:
                critical('Undexpected error: incorrect line in the coverageBed output:\n' + next_line)

            if next_line.startswith('all'):
                max_depth = max(max_depth, depth)
                total_bed_size += bases
                extra_fields = []
            else:
                start, end = map(int, line_tokens[1:3])
                extra_fields = line_tokens[3:-4]

            line_region_key_tokens = (None, chrom, start, end)

            if regions == [] or hash(line_region_key_tokens) != regions[-1].key():
                region = Region(sample=None, chrom=chrom,
                                start=start, end=end, size=region_size,
                                extra_fields=extra_fields)
                regions.append(region)

                _total_regions_count += 1

                if _total_regions_count > 0 and _total_regions_count % 100000 == 0:
                    info('  Processed {0:,} regions'.format(_total_regions_count))

            regions[-1].add_bases_for_depth(depth, bases)

    if _total_regions_count % 100000 != 0:
        info('Processed {0:,} regions'.format(_total_regions_count))

    return regions[:-1], regions[-1], max_depth, total_bed_size


def sort_bed(cnf, bed_fpath):
    bedtools = get_tool_cmdline(cnf, 'bedtools')
    cmdline = '{bedtools} sort -i {bed_fpath}'.format(**locals())
    output_fpath = intermediate_fname(cnf, bed_fpath, 'sorted')
    call(cnf, cmdline, output_fpath)
    return output_fpath


def intersect_bed(cnf, bed1, bed2):
    bed1_fname, _ = splitext_plus(basename(bed1))
    bed2_fname, _ = splitext_plus(basename(bed2))
    output_fpath = join(cnf['work_dir'], bed1_fname + '__' + bed2_fname + '.bed')
    bedtools = get_tool_cmdline(cnf, 'bedtools')
    cmdline = '{bedtools} intersect -a {bed1} -b {bed2} -u'.format(**locals())
    call(cnf, cmdline, output_fpath)
    return output_fpath


# TODO very slow :(
def number_of_mapped_reads(cnf, bam):
    samtools = get_tool_cmdline(cnf, 'samtools')
    cmdline = '{samtools} view -c -F 4 {bam}'.format(**locals())
    res = call_check_output(cnf, cmdline)
    return int(res)


# TODO very slow :(
def number_of_unmapped_reads(cnf, bam):
    samtools = get_tool_cmdline(cnf, 'samtools')
    cmdline = '{samtools} view -c -f 4 {bam}'.format(**locals())
    res = call_check_output(cnf, cmdline)
    return int(res)


# TODO very slow :(
def number_of_reads(cnf, bam):
    samtools = get_tool_cmdline(cnf, 'samtools')
    cmdline = '{samtools} view -c {bam}'.format(**locals())
    res = call_check_output(cnf, cmdline)
    return int(res)


# TODO very slow :(
def number_mapped_reads_on_target(cnf, bed, bam):
    samtools = get_tool_cmdline(cnf, 'samtools')
    cmdline = '{samtools} view -c -F 4 -L {bed} {bam}'.format(**locals())
    res = call_check_output(cnf, cmdline)
    return int(res)


# TODO very slow :(
def number_bases_in_aligned_reads(cnf, bam):
    samtools = get_tool_cmdline(cnf, 'samtools')
    cmdline = '{samtools} depth {bam}'.format(**locals())
    proc = call_pipe(cnf, cmdline)
    count = 0
    while True:
        coverage_line = proc.stdout.readline()
        if coverage_line:
            values = coverage_line.strip().split('\t')
            count += int(values[2])
    return count


# TODO how to pass the data stream to samtools vs. creating file
def get_padded_bed_file(cnf, bed, genome, padding):
    bedtools = get_tool_cmdline(cnf, 'bedtools')
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


