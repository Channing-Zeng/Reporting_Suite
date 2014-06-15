import copy
from itertools import izip, chain, repeat
from os.path import join, basename

from source.logger import step_greetings, critical
from source.targetcov.Region import Region
from source.utils import intermediate_fname, get_tool_cmdline, info, err, \
    call_check_output, call_pipe, call, format_integer, format_decimal
from source.transaction import file_transaction
from source.utils_from_bcbio import splitext_plus


def run_target_cov(cnf, bam, bed):
    summary_report_fpath = None
    gene_report_fpath = None

    info('Calculation of coverage statistics for the regions in the input BED file...')
    amplicons, combined_region, max_depth, total_bed_size = _bedcoverage_hist_stats(cnf, bam, bed)

    chr_len_fpath = cnf['genome']['chr_lengths']
    exons_bed = cnf['genome']['exons']
    genes_bed = cnf['genome']['genes']

    if 'summary' in cnf['reports']:
        step_greetings('Target coverage summary report')
        summary_report_fpath = join(cnf['output_dir'], cnf['name'] + '.targetseq.summary.txt')
        _run_header_report(
            cnf, summary_report_fpath,
            bed, bam, chr_len_fpath,
            cnf['coverage_reports']['depth_thresholds'], cnf['padding'],
            combined_region, max_depth, total_bed_size)

    # if 'amplicons' in options['reports']:
    #     step_greetings('Coverage report for the input BED file regions')
    #     amplicons_report_fpath = join(output_dir, sample_name + '.targetseq.details.capture.txt')
    #     run_amplicons_cov_report(cnf, amplicons_report_fpath, sample_name, depth_threshs, amplicons)

    if 'genes' in cnf['reports']:
        if not genes_bed or not exons_bed:
            if cnf['reports'] == 'genes':
                critical('Error: no genes or exons specified for the genome in system config, '
                         'cannot run per-exon report.')
            else:
                err('Warning: no genes or exons specified for the genome in system config, '
                    'cannot run per-exon report.')
        else:
            # log('Annotating amplicons.')
            # annotate_amplicons(amplicons, genes_bed)

            info('Getting the gene regions that intersect with our capture panel.')
            bed = intersect_bed(cnf, genes_bed, bed)
            info('Getting the exons of the genes.')
            bed = intersect_bed(cnf, exons_bed, bed)
            info('Sorting final exon BED file.')
            bed = sort_bed(cnf, bed)

            info('Calculation of coverage statistics for exons of the genes ovelapping with the input regions...')
            exons, _, _, _ = _bedcoverage_hist_stats(cnf, bam, bed)
            for exon in exons:
                exon.gene_name = exon.extra_fields[0]

            gene_report_fpath = join(cnf['output_dir'], cnf['name'] + '.targetseq.details.gene.txt')
            _run_region_cov_report(cnf, gene_report_fpath, cnf['name'], cnf['coverage_reports']['depth_thresholds'],
                                  amplicons, exons)

    return summary_report_fpath, gene_report_fpath


def _run_header_report(cnf, result_fpath,
                       bed, bam, chr_len_fpath,
                       depth_thresholds, padding,
                       combined_region, max_depth, total_bed_size):
    stats = []

    def append_stat(stat):
        stats.append(stat)
        info(stat)

    info('* General coverage statistics *')
    info('Getting number of reads...')
    v_number_of_reads = number_of_reads(cnf, bam)
    append_stat(format_integer('Reads', v_number_of_reads))

    info('Getting number of mapped reads...')
    v_mapped_reads = number_of_mapped_reads(cnf, bam)
    append_stat(format_integer('Mapped reads', v_mapped_reads))
    append_stat(format_integer('Unmapped reads', v_number_of_reads - v_mapped_reads))

    v_percent_mapped = 100.0 * v_mapped_reads / v_number_of_reads if v_number_of_reads else None
    append_stat(format_decimal('Percentage of mapped reads', v_percent_mapped, '%'))
    info('')

    info('* Target coverage statistics *')
    append_stat(format_integer('Bases in target', total_bed_size))

    bases_within_threshs, avg_depth, std_dev, percent_within_normal = combined_region.sum_up(depth_thresholds)

    v_covered_bases_in_targ = bases_within_threshs.items()[0][1]
    append_stat(format_integer('Covered bases in target', v_covered_bases_in_targ))

    v_percent_covered_bases_in_targ = 100.0 * v_covered_bases_in_targ / total_bed_size if total_bed_size else None
    append_stat(format_decimal('Percentage of target covered by at least 1 read', v_percent_covered_bases_in_targ, '%'))
    info('Getting number of mapped reads on target...')

    v_mapped_reads_on_target = number_mapped_reads_on_target(cnf, bed, bam)
    append_stat(format_integer('Reads mapped on target', v_mapped_reads_on_target))

    v_percent_mapped_on_target = 100.0 * v_mapped_reads_on_target / v_mapped_reads if v_mapped_reads else None
    append_stat(format_decimal('Percentage of reads mapped on target ', v_percent_mapped_on_target, '%'))

    info('Making bed file for padded regions...')
    padded_bed = get_padded_bed_file(cnf, bed, chr_len_fpath, padding)

    info('Getting number of mapped reads on padded target...')
    v_reads_on_padded_targ = number_mapped_reads_on_target(cnf, padded_bed, bam)
    append_stat(format_integer('Reads mapped on padded target', v_reads_on_padded_targ))

    v_percent_mapped_on_padded_target = 100.0 * v_reads_on_padded_targ / v_mapped_reads if v_mapped_reads else None
    append_stat(format_decimal('Percentage of reads mapped on padded target', v_percent_mapped_on_padded_target, '%'))

    # v_aligned_read_bases = number_bases_in_aligned_reads(bam)
    # append_stat(format_integer('Total aligned bases in reads', v_aligned_read_bases))

    v_read_bases_on_targ = avg_depth * total_bed_size  # sum of all coverages
    append_stat(format_integer('Read bases mapped on target', v_read_bases_on_targ))

    info('')
    append_stat(format_decimal('Average target coverage depth', avg_depth))
    append_stat(format_decimal('Std. dev. of target coverage depth', std_dev))
    append_stat(format_integer('Maximum target coverage depth', max_depth))
    append_stat(format_decimal('Percentage of target within 20% of mean depth',
                               percent_within_normal, '%'))

    # v_percent_read_bases_on_targ = 100.0 * v_read_bases_on_targ / v_aligned_read_bases \
    #     if v_aligned_read_bases else None
    # format_decimal('Percent bases in reads on target', v_percent_read_bases_on_targ, '%'),

    # format_integer('Bases covered (at least 1x) in target', bases_per_depth[1]),

    # for depth, bases in bases_per_depth.items():
    #     append_stat(format_integer('Bases on target covered at least by ' + str(depth) +
    #                                ' read' + ('s' if depth != 1 else ''), bases))

    for depth, bases in bases_within_threshs.items():
        percent = 100.0 * bases / total_bed_size if total_bed_size else 0
        append_stat(format_decimal('Part of target covered at least by ' + str(depth) +
                                   'x', percent, '%'))

    max_len = max(len(l.rsplit(':', 1)[0]) for l in stats)
    with file_transaction(cnf.tmp_dir, result_fpath) as tx, open(tx, 'w') as out:
        for l in stats:
            text, val = l.rsplit(':', 1)
            # spaces = ' ' * (max_len - len(text) + 1)
            spaces = '\t'
            out.write(text + spaces + val + '\n')

    info('')
    info('Result: ' + result_fpath)
    return result_fpath


def _run_region_cov_report(cnf, report_fpath, sample_name, depth_threshs,
                          amplicons, exons):
    for ampl in amplicons:
        ampl.feature = 'Amplicon'
        ampl.sample = sample_name

    for exon in exons:
        exon.feature = 'Exon'
        exon.sample = sample_name

    exon_genes = _get_exon_genes(cnf, exons)
    amplicon_genes_by_name = _get_amplicon_genes(amplicons, exon_genes)

    result_regions = []
    for exon_gene in exon_genes:
        for exon in exon_gene.subregions:
            result_regions.append(exon)
        result_regions.append(exon_gene)

        amplicon_gene = amplicon_genes_by_name.get(exon_gene.gene_name)
        if amplicon_gene:
            for amplicon in amplicon_gene.subregions:
                result_regions.append(amplicon)
            result_regions.append(amplicon_gene)

    return _build_regions_cov_report(
        cnf, report_fpath, depth_threshs, result_regions)


def _get_amplicon_genes(amplicons, exon_genes):
    amplicon_genes_by_name = dict()

    for exon_gene in exon_genes:
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

    return amplicon_genes_by_name


# def run_amplicons_cov_report(cnf, report_fpath, sample_name, depth_threshs, regions):
#     for region in regions:
#         region.feature = 'Amplicon'
#         region.sample = sample_name
#
#     # exons_and_genes = add_genes_cov_analytics(regions, gene_pos=0)
#
#     return build_regions_cov_report(cnf, report_fpath, depth_threshs, regions)
#
#
# def run_exons_cov_report(cnf, report_fpath, sample_name, depth_threshs, regions):
#     extra_fields = ['Gene']
#
#     for region in regions:
#         region.feature = 'Exon'
#         region.sample = sample_name
#         region.extra_fields = region.extra_fields[:1]
#
#     exons_with_genes = add_genes_cov_analytics(regions, gene_pos=0)
#
#     return build_regions_cov_report(
#         cnf, report_fpath, depth_threshs,
#         exons_with_genes, extra_headers=extra_fields)


def _get_exon_genes(cnf, subregions):
    genes_by_name = dict()

    for exon in subregions:
        if not exon.gene_name:
            err('No gene name info in the record: ' +
                str(exon) + '. Skipping.')
            continue

        gene = genes_by_name.get(exon.gene_name)
        if gene is None:
            gene = Region(
                sample=exon.sample, chrom=exon.chrom,
                start=exon.start, end=exon.end, size=0,
                feature='Gene-' + exon.feature,
                gene_name=exon.gene_name)
            genes_by_name[exon.gene_name] = gene
        gene.add_subregion(exon)

    sorted_genes = sorted(genes_by_name.values(),
                          key=lambda g: (g.chrom, g.start, g.end))
    return sorted_genes


# def add_gene_names_to_amplicons(amplicons, genes):
#     for gene in genes:
#         for amplicon in amplicons:
#             if gene.intersect(amplicon):
#                 amplicon.gene_name = gene.gene_name


def _build_regions_cov_report(cnf, report_fpath, depth_threshs, regions,
                             extra_headers=list()):
    header_fields = ['SAMPLE', 'Chr', 'Start', 'End', 'Gene', 'Feature', 'Size',
                     'Mean Depth', 'Standard Dev.', 'Within 20% of Mean'] +\
                    ['{}x'.format(thres) for thres in depth_threshs]
    max_lengths = map(len, header_fields)

    all_values = []

    for i, region in enumerate(regions):
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

    with file_transaction(cnf['tmp_dir'], report_fpath) as tx:
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
    cmdline = '{bedtools} coverage -abam {bam} -b {bed} ' \
              '-hist'.format(**locals())

    _total_regions_count = 0

    for next_line in call_pipe(cnf, cmdline).stdout:
        if not next_line.strip() or next_line.startswith('#'):
            continue

        line_tokens = next_line.strip().split()
        chrom = line_tokens[0]
        start, end = None, None
        try:
            depth, bases, region_size = map(int, line_tokens[-4:-1])
        except:
            critical('Undexpected error: incorrect line in coverageBed output:\n' + next_line)

        if next_line.startswith('all'):
            max_depth = max(max_depth, depth)
            total_bed_size += bases
            extra_fields = []
        else:
            start, end = map(int, line_tokens[1:3])
            extra_fields = line_tokens[3:-4]

        line_region_key_tokens = (None, chrom, start, end)

        if regions == [] or hash(line_region_key_tokens) != regions[-1].key():
            _total_regions_count += 1

            region = Region(sample=None, chrom=chrom,
                            start=start, end=end, size=region_size,
                            extra_fields=extra_fields)
            regions.append(region)

        regions[-1].add_bases_for_depth(depth, bases)

        if _total_regions_count > 0 and _total_regions_count % 100000 == 0:
            info('processed %i regions' % _total_regions_count)

    if _total_regions_count % 100000 != 0:
        info('processed %i regions' % _total_regions_count)

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


