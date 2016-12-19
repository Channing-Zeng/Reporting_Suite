#!/usr/bin/env python
import bcbio_postproc
import sys
from collections import defaultdict
from os.path import basename, join
from optparse import OptionParser

from source import verify_file, BaseSample
from source.config import Config
from source.file_utils import file_transaction, safe_mkdir
from source.logger import info, critical, warn
from source.prepare_args_and_cnf import determine_sys_cnf
from source.targetcov.bam_and_bed_utils import get_bedgraph_coverage
from source.utils import mean, median
from source.variants.vcf_processing import bgzip_and_tabix


def write_coverage(cnf, output_dir, chrom, depths_by_pos, cov_thresholds, sample_index=None):
    coverage_data_fpath = join(output_dir, chrom + '.txt')
    if not cnf.reuse_intermediate or (not verify_file(coverage_data_fpath, silent=True) and
                                      not verify_file(coverage_data_fpath + '.gz', silent=True)):
        chrom_num = chrom.replace('chr', '')
        with file_transaction(cnf.work_dir, coverage_data_fpath) as tx:
            with open(tx, 'w') as f:
                fs = ['#chrom', 'pos', 'mean', 'median'] + [str(t) for t in cov_thresholds]
                f.write('\t'.join(fs) + '\n')
                sorted_positions = sorted(depths_by_pos.keys())
                for pos in sorted_positions:
                    depths = depths_by_pos[pos] if sample_index is None else [depths_by_pos[pos][sample_index]]
                    mean_coverage = mean(depths)
                    median_coverage = median(depths)
                    pcnt_samples_ge_threshold = [mean([1 if d >= t else 0 for d in depths]) for t in cov_thresholds]
                    res_line = chrom_num + '\t' + str(pos) + '\t' + str(mean_coverage) + '\t' + str(median_coverage)
                    for pcnt_samples in pcnt_samples_ge_threshold:
                        res_line += '\t' + str(pcnt_samples)
                    f.write(res_line + '\n')
    bgzip_and_tabix(cnf, coverage_data_fpath, tabix_parameters='-p bed')


def get_regions_coverage(cnf, samples):
    cov_thresholds = [1, 5, 10, 15, 20, 25, 30, 50, 100]
    depths_by_pos = defaultdict(lambda : [0] * len(samples))
    info()
    info('Coverage to bedgraph for ' + cnf.chrom)
    coverage_fpaths = []
    for index, sample in enumerate(samples):
        coverage_fpath = join(cnf.work_dir, sample.name + '_' + cnf.chrom + '.bedgraph')
        coverage_fpath = get_bedgraph_coverage(cnf, sample.bam, chr_len_fpath=cnf.chr_len_fpath, bed_fpath=cnf.bed, output_fpath=coverage_fpath, exit_on_error=False)
        if coverage_fpath and verify_file(coverage_fpath):
            coverage_fpaths.append(coverage_fpath)
            for line in open(coverage_fpath):
                if line.startswith('#'):
                    continue
                chrom, start, end, depth = line.split('\t')
                start, end, depth = map(int, (start, end, depth))
                for pos in xrange(start, end):
                    depths_by_pos[pos][index] = depth

    info()
    if not coverage_fpaths:
        warn(cnf.chrom + ' is not covered in all samples')
        return None

    info()
    info('Writing coverage for ' + cnf.chrom)
    write_coverage(cnf, cnf.output_dir, cnf.chrom, depths_by_pos, cov_thresholds)
    for index, sample in enumerate(samples):
        info('Writing coverage for ' + sample.name + ', ' + chrom)
        sample_output_dirpath = join(cnf.output_dir, sample.name)
        output_fpath = join(sample_output_dirpath, chrom + '.txt.gz')
        if cnf.reuse_intermediate and verify_file(output_fpath, silent=True):
            info(output_fpath + ' exists, reusing')
            continue
        write_coverage(cnf, sample_output_dirpath, cnf.chrom, depths_by_pos, cov_thresholds, sample_index=index)
        if not verify_file(output_fpath, silent=True):
            warn(sample.name + ' has no coverage at chromosome ' + chrom)
    return depths_by_pos


def main():
    info(' '.join(sys.argv))
    info()
    parser = OptionParser(usage='Usage: ' + basename(__file__) + ' --bed BED_file --bam BAM_file -g hg19 -o Output_BEDGRAPH_file '
                                                                 '--work-dir work_directory --chr chromosome')
    parser.add_option('-o', dest='output_dir')
    parser.add_option('--samples', dest='sample_names')
    parser.add_option('--bams', dest='bams')
    parser.add_option('--vcf', dest='vcf_fpath')
    parser.add_option('--chr', dest='chrom')
    parser.add_option('--bed', dest='bed', help='BED file.')
    parser.add_option('-g', '--genome', dest='chr_len_fpath', help='File with chromosomes lengths.')
    parser.add_option('--work-dir', dest='work_dir', help='Work directory.')
    (opts, args) = parser.parse_args(sys.argv[1:])

    cnf = Config(opts.__dict__, determine_sys_cnf(opts), {})
    samples = [BaseSample(sample_name, None, bam=bam) for (sample_name, bam) in zip(cnf.sample_names.split(','), cnf.bams.split(','))]

    if not cnf.output_dir or not cnf.bams:
        critical(parser.usage)

    safe_mkdir(cnf.output_dir)
    safe_mkdir(cnf.work_dir)
    get_regions_coverage(cnf, samples)
    info('Done.')


if __name__ == '__main__':
    main()
