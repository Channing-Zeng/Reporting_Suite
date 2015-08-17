from os.path import join, basename, splitext
from time import sleep

from source.calling_process import call
from source.file_utils import splitext_plus
from source.logger import critical, info, warn, send_email, err
from source.targetcov.Region import Region
from source.targetcov.bam_and_bed_utils import count_bed_cols, bedtools_version
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
    info()
    info('Anylising bedcoverage output...')
    return summarize_bedcoverage_hist_stats(bedcov_output, sample_name, bed_col_num)


# TODO:
# split bam and bed by chromosome
# run bedtools hist
# merge. for general stats, merge more sofisticated.


def run_bedcoverage_hist_stats(cnf, bed, bam):
    bedtools = get_system_path(cnf, 'bedtools')
    chr_lengths = get_chr_len_fpath(cnf)

    bedcov_output = join(cnf.work_dir,
        splitext_plus(basename(bed))[0] + '_' +
        splitext_plus(basename(bam))[0] + '_bedcov_output.txt')

    v = bedtools_version(bedtools)
    if v and v >= 24:
        cmdline = '{bedtools} coverage -sorted -g {chr_lengths} -a {bed} -b {bam} -hist'.format(**locals())
    else:
        cmdline = '{bedtools} coverage -abam {bam} -b {bed} -hist'.format(**locals())

    res = None
    tries = 0
    MAX_TRIES = 10
    err_fpath = join(cnf.work_dir, 'bedtools_cov_' + splitext(basename(bedcov_output))[0] + '.err')
    while True:
        stderr_dump = []
        res = call(cnf, cmdline, bedcov_output, stderr_dump=stderr_dump, exit_on_error=False)
        if res is not None:
            return res
        else:
            tries += 1
            msg = 'bedtools coverage crashed:\n' + cmdline + '\n' + \
                  (''.join(['\t' + l for l in stderr_dump]) if stderr_dump else '')
            if tries < MAX_TRIES:
                msg += '\n\nRerunning in 120 minutes (tries ' + str(tries) + '/10)'

            send_email(msg=msg,
                       subj='bedtools coverage crashed [' + str(cnf.project_name) + ']',
                       only_me=True)
            err(msg)
            if tries == MAX_TRIES:
                break
            sleep(120 * 60)
            info()

    return bedcov_output


def summarize_bedcoverage_hist_stats(bedcov_output_fpath, sample_name, bed_col_num):
    """
    :param bedcov_output_fpath: file path
    :param sample_name: sting
    :param bed_col_num: int
    :return: regions, total_region, max_depth
        regions: [Region(
            sample_name: string
            chrom: string
            start: int
            end: int
            size: int
            gene_name: string
            bases_by_depth: dict(depth:int->bases:int)
            extra_fields: [string])]
        total_region: Region(--|--)
        max_depth: int
    """

    _total_regions_count = 0
    regions, max_depth, total_bed_size = [], 0, 0

    with open(bedcov_output_fpath) as f:
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

                # if _total_regions_count > 0 and _total_regions_count % 100000 == 0:
                #     info('  Processed {0:,} regions'.format(_total_regions_count))

            regions[-1].add_bases_for_depth(depth, bases)

            if regions[-1].min_depth is None:
                regions[-1].min_depth = depth  # depth values go from lowest to highest, so if we are meeting the first record for this region, the depth would be the lowest

    # if _total_regions_count % 100000 != 0:
    #     info('  Processed {0:,} regions'.format(_total_regions_count))

    total_region = regions[-1]
    regions = regions[:-1]
    # info('Sorting genes...')
    # regions = sorted(regions[:-1], key=Region.get_order_key)

    return regions, total_region, max_depth