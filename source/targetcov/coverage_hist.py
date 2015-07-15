from os.path import join, basename

from source.calling_process import call
from source.file_utils import splitext_plus
from source.logger import critical, info, warn
from source.targetcov.Region import Region
from source.targetcov.bam_and_bed_utils import count_bed_cols
from source.tools_from_cnf import get_system_path
from source.utils import get_chr_len_fpath


def bedcoverage_hist_stats(cnf, sample_name, bam, bed, reuse=False):
    if not bam or not bed:
        info()
        msgs = []
        if not bam: msgs.append('BAM file is required.')
        if not bed: msgs.append('BED file is required.')
        if msgs:
            critical(msgs)

    bedcov_output = run_bedcoverage_hist_stats(cnf, bed, bam)
    bed_col_num = count_bed_cols(bed)
    return summarize_bedcoverage_hist_stats(bedcov_output, sample_name, bed_col_num)


def run_bedcoverage_hist_stats(cnf, bed, bam):
    bedtools = get_system_path(cnf, 'bedtools')
    chr_lengths = get_chr_len_fpath(cnf)

    cmdline = '{bedtools} coverage -sorted -g {chr_lengths} -a {bed} -b {bam} -hist'.format(**locals())
    bedcov_output = join(cnf.work_dir,
        splitext_plus(basename(bed))[0] + '_' +
        splitext_plus(basename(bam))[0] + '_bedcov_output.txt')
    # if reuse and file_exists(bedcov_output) and verify_file(bedcov_output):
    #     pass
    # else:
    res = call(cnf, cmdline, bedcov_output, exit_on_error=False)
    if not res:
        info()
        warn('Could not run bedtools, maybe old version, trying without -sorted -g [genome]')
        cmdline = '{bedtools} coverage -abam {bam} -b {bed} -hist'.format(**locals())
        info()
        res = call(cnf, cmdline, bedcov_output)
    return bedcov_output


def summarize_bedcoverage_hist_stats(bedcov_output, sample_name, bed_col_num):
    _total_regions_count = 0
    regions, max_depth, total_bed_size = [], 0, 0

    info()
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
                gene_name = line_tokens[3] if bed_col_num > 3 else None
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

    return regions, total_region, max_depth