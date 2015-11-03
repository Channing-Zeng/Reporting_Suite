from genericpath import isfile
from os.path import join, basename, splitext
from random import random
from time import sleep

from source.calling_process import call
from source.file_utils import splitext_plus, add_suffix, safe_mkdir, file_transaction, verify_file
from source.logger import critical, info, warn, send_email, err
# from source.ngscat import coverageHisto
from source.targetcov.bam_and_bed_utils import verify_bam, verify_bed
from source.targetcov.Region import Region
from source.targetcov.bam_and_bed_utils import count_bed_cols, bedtools_version
from source.tools_from_cnf import get_system_path
from source.utils import get_chr_len_fpath


# def ngscat_bedcov_hist_stats(cnf, sample_name, bam_fpath, bed_fpath, reuse=False):
#     cov_out_fpath = join(cnf.work_dir, basename(bam_fpath)[:-len('bam')] + 'coverage')
#     bedgraph_fpath = join(cnf.work_dir, basename(bam_fpath)[:-len('bam')] + 'bedgraph')
#
#     info('Coveragefile = ' + cov_out_fpath)
#     from source.ngscat import bam_file
#     bam = bam_file.BamFile(bam_fpath, 'rb')
#
#     info('Launching coverageBed...')
#     bam.myCoverageBed(cnf, bed_fpath, writeToFile=cov_out_fpath, bedGraphFile=bedgraph_fpath)
#     info('\tDone.')
#
#     cov_hist_dirpath = join(cnf.work_dir, 'cov_hist')
#     safe_mkdir(cov_hist_dirpath)
#     info('Launching coverage distribution calculation...')
#     coverageHisto.histo_CV(cov_out_fpath, cov_hist_dirpath)


class BedCov:
    def __init__(self, chrom, chrom_len, bedcov_output_fpath):
        self.chrom = chrom
        self.chrom_len = chrom_len
        self.bedcov_output_fpath = bedcov_output_fpath


def bedcoverage_hist_stats(cnf, sample_name, bam_fpath, bed_fpath, reuse=False):
    if not bam_fpath or not bed_fpath:
        info()
        msgs = []
        if not bam_fpath: msgs.append('BAM file is required.')
        if not bed_fpath: msgs.append('BED file is required.')
        if msgs:
            critical(msgs)

    bed_col_num = count_bed_cols(bed_fpath)

    bedcov_output_fpath = join(cnf.work_dir,
        splitext_plus(basename(bed_fpath))[0] + '__' +
        splitext_plus(basename(bam_fpath))[0] + '_bedcov_output.txt')

    regions = []
    bamtools = get_system_path(cnf, 'bamtools')
    if not bamtools:
        info('Running bedcoverage -hist...')
        bedcov_output_fpath = launch_bedcoverage_hist(cnf, bed_fpath, bam_fpath, bedcov_output_fpath)
        info()
        info('Analysing bedcoverage -hist output...')
        regions = summarize_bedcoverage_hist_stats(bedcov_output_fpath, sample_name, bed_col_num)

    else:
        chroms = [l.split()[0] for l in open(get_chr_len_fpath(cnf)).readlines()]

        stub = join(cnf.work_dir, basename(splitext_plus(bam_fpath)[0]))
        if cnf.reuse_intermediate and all(verify_bam(stub + '.REF_' + chrom + '.bam', silent=True) for chrom in chroms):
            info('BAM ' + bam_fpath + ' is split, reusing...')
        else:
            info('Splitting the BAM file, writing as ' + stub + '.REF_#.bam')
            cmdline = '{bamtools} split -in {bam_fpath} -stub {stub} -reference'.format(**locals())
            call(cnf, cmdline)

        for chrom in chroms:
            # chrom_bed_fpath = add_suffix(bed_fpath, chrom)
            chrom_bed_fpath = '<(grep "^{chrom}\t" {bed_fpath})'.format(**locals())  # Note: "\t" will be converted to "	"
            # grep = get_system_path(cnf, 'grep')
            # cmdl = '{grep} "^{chrom}" {bed_fpath}'.format(**locals())
            # call(cnf, cmdl, output_fpath=chrom_bed_fpath, exit_on_error=False)
            chrom_bam_fpath = stub + '.REF_' + chrom + '.bam'
            if verify_bam(chrom_bam_fpath, silent=True):
                bedcov_output_fpath = launch_bedcoverage_hist(cnf, chrom_bed_fpath, chrom_bam_fpath)
                if not verify_file(bedcov_output_fpath):
                    info('No coverage for ' + chrom)
                else:
                    info('Anylising bedcoverage output for ' + str(chrom) + '...')
                    rs = summarize_bedcoverage_hist_stats(bedcov_output_fpath, sample_name, bed_col_num)
                    regions.extend(rs)

        # with open(bedcov_output_fpath, 'w') as f:
        #     for chrom, bedov_output in bedcov_by_chrom.items():

    # for r in regions:
    #TODO Sync real regions in {bed} and {regions}

    return regions


# TODO:
# split bam and bed by chromosome
# run bedtools hist
# merge. for general stats, merge more sofisticated.


# def grep(cnf, fpath, pattern):
#     res =


def launch_bedcoverage_hist(cnf, bed_fpath, bam_fpath, bedcov_output_fpath=None):
    # import pybedtools
    # bed = pybedtools.BedTool(bed_fpath)
    # bam = pybedtools.BedTool(bam_fpath)
    # return bed.coverage(bam)

    bedtools = get_system_path(cnf, 'bedtools')
    chr_lengths = get_chr_len_fpath(cnf)

    bedcov_output_fpath = bedcov_output_fpath or join(cnf.work_dir,
        splitext_plus(basename(bed_fpath))[0] + '__' +
        splitext_plus(basename(bam_fpath))[0] + '_bedcov_output.txt')

    v = bedtools_version(bedtools)
    if v and v >= 24:
        cmdline = '{bedtools} coverage -sorted -g {chr_lengths} -a {bed_fpath} -b {bam_fpath} -hist'.format(**locals())
    else:
        cmdline = '{bedtools} coverage -abam {bam_fpath} -b {bed_fpath} -hist'.format(**locals())

    res = None
    tries = 0
    MAX_TRIES = 2
    WAIT_MINUTES = int(random() * 60) + 30
    err_fpath = join(cnf.work_dir, 'bedtools_cov_' + splitext(basename(bedcov_output_fpath))[0] + '.err')
    while True:
        stderr_dump = []
        res = call(cnf, cmdline, bedcov_output_fpath, stderr_dump=stderr_dump, exit_on_error=False)
        if isfile(bedcov_output_fpath):
            return bedcov_output_fpath
        else:
            tries += 1
            msg = 'bedtools coverage crashed:\n' + cmdline + ' > ' + bedcov_output_fpath + '\n' + \
                  (''.join(['\t' + l for l in stderr_dump]) if stderr_dump else '')
            if tries < MAX_TRIES:
                msg += '\n\nRerunning in ' + str(WAIT_MINUTES) + ' minutes (tries ' + str(tries) + '/' + str(MAX_TRIES) + ' )'

            send_email(msg_other=msg,
                       subj='bedtools coverage crashed [' + str(cnf.project_name) + ']',
                       only_me=True)
            err(msg)
            if tries == MAX_TRIES:
                break
            sleep(WAIT_MINUTES * 60)
            info()

    return bedcov_output_fpath


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

    # for next_line in bedcov:
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

    regions = regions[:-1]
    # info('Sorting genes...')
    # regions = sorted(regions[:-1], key=Region.get_order_key)

    return regions  #, total_region, max_depth
