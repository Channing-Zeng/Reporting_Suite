#!/usr/bin/env python
import bcbio_postproc
import os
import sys
from collections import defaultdict
from optparse import OptionParser

from os.path import join, basename

from source.bcbio.bcbio_structure import BCBioStructure, process_post_bcbio_args
from source.calling_process import call
from source.file_utils import safe_mkdir, adjust_path, file_transaction, verify_file, add_suffix, which
from source.logger import critical, info
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug, set_up_log
from source.qsub_utils import wait_for_jobs, submit_job
from source.targetcov.bam_and_bed_utils import call_sambamba
from source.tools_from_cnf import get_script_cmdline, get_system_path
from source.utils import mean, get_chr_len_fpath, is_us
from source.utils import median
from source.variants.filtering import combine_vcfs
from source.variants.vcf_processing import bgzip_and_tabix
from tools.txt2vcf import convert_txt_to_vcf


EXAC_FILES_DIRECTORY = '../exac_data/'


def run_bedtools_use_grid(cnf, bam_by_key, bed_fpath):
    output_by_key = dict()
    not_submitted_bams = bam_by_key.values()
    chr_len_fpath = get_chr_len_fpath(cnf)
    while not_submitted_bams:
        jobs_to_wait = []
        submitted_bams = []
        reused_bams = []

        for key, bam in bam_by_key.iteritems():
            if bam not in not_submitted_bams:
                continue
            uniq_name, dict_key = key
            output_fpath = join(cnf.work_dir, uniq_name + '_coverage.txt')
            output_by_key[dict_key] = output_fpath

            if cnf.reuse_intermediate and verify_file(output_fpath, silent=True):
                info(output_fpath + ' exists, reusing')
                reused_bams.append(bam)
                continue
            else:
                cmdline = get_script_cmdline(cnf, 'python', join('tools', 'get_region_coverage.py'), is_critical=True)
                cmdline += ' --bam {bam} -o {output_fpath} -g {chr_len_fpath} --work-dir {cnf.work_dir}'.format(**locals())
                if bed_fpath:
                     cmdline += ' --bed {bed_fpath}'.format(**locals())
                j = submit_job(cnf, cmdline,  basename(bam) + '_coverage')
                info()
                submitted_bams.append(bam)

                if not j.is_done:
                    jobs_to_wait.append(j)
                if len(jobs_to_wait) >= cnf.threads:
                    break
        if jobs_to_wait:
            info('Submitted ' + str(len(jobs_to_wait)) + ' jobs, waiting...')
            jobs_to_wait = wait_for_jobs(cnf, jobs_to_wait)
        else:
            info('No jobs to submit.')
        not_submitted_bams = [bam for bam in not_submitted_bams if
                                      bam not in submitted_bams and
                                      bam not in reused_bams]

    return output_by_key


def get_regions_depth(cnf, samples):
    depths_by_pos = defaultdict(lambda : defaultdict(list))
    bam_by_key = dict()
    for s in samples:
        bam_by_key[(s.name, s.name)] = s.bam
    output_by_sample = run_bedtools_use_grid(cnf, bam_by_key, cnf.bed)
    info()
    info('Parsing bedtools output...')
    for sample, coverage_fpath in output_by_sample.iteritems():
        for line in open(coverage_fpath):
            if line.startswith('#'):
                continue
            chrom, start, end, depth = line.split('\t')
            start, end, depth = map(int, (start, end, depth))
            for pos in xrange(start, end):
                depths_by_pos[chrom][pos].append(depth)
    return depths_by_pos


def dedup_and_sort_bams_use_grid(cnf, samples, do_sort=False):
    jobs_to_wait = []
    not_submitted_samples = [sample for sample in samples]
    done_samples = []
    while not_submitted_samples:
        jobs_to_wait = []
        submitted_samples = []
        reused_samples = []

        for sample in not_submitted_samples:
            if do_sort:
                output_bam_fpath = join(cnf.work_dir, sample.name + '.dedup.sorted.bam')
            else:
                output_bam_fpath = join(cnf.work_dir, sample.name + '.dedup.bam')

            if cnf.reuse_intermediate and verify_file(output_bam_fpath, silent=True):
                info(output_bam_fpath + ' exists, reusing')
                sample.bam = output_bam_fpath
                done_samples.append(sample)
                reused_samples.append(sample)
                continue
            else:
                if do_sort:
                    cmdline = 'sort {sample.bam} -o {output_bam_fpath}'.format(**locals())
                    j = call_sambamba(cnf, cmdline, output_fpath=output_bam_fpath, bam_fpath=sample.bam,
                                    use_grid=True, command_name='sort', sample_name=sample.name)
                else:
                    cmdline = 'view -f bam -F "not duplicate and not failed_quality_control" {sample.bam}'.format(**locals())
                    j = call_sambamba(cnf, cmdline, output_fpath=output_bam_fpath, bam_fpath=sample.bam,
                                      use_grid=True, command_name='dedup', sample_name=sample.name)
                info()
                sample.bam = output_bam_fpath
                done_samples.append(sample)
                submitted_samples.append(sample)

                if not j.is_done:
                    jobs_to_wait.append(j)
                if len(jobs_to_wait) >= cnf.threads:
                    break
        if jobs_to_wait:
            info('Submitted ' + str(len(jobs_to_wait)) + ' jobs, waiting...')
            jobs_to_wait = wait_for_jobs(cnf, jobs_to_wait)
        else:
            info('No jobs to submit.')
        not_submitted_samples = [sample for sample in not_submitted_samples if
                                          sample not in submitted_samples and
                                          sample not in reused_samples]
    return done_samples


def split_bam_files_use_grid(cnf, samples, combined_vcf_fpath, exac_features_fpath, exac_venv_pythonpath):
    samples = dedup_and_sort_bams_use_grid(cnf, samples, do_sort=False)
    samples = dedup_and_sort_bams_use_grid(cnf, samples, do_sort=True)

    chromosomes = ['chr%s' % x for x in range(1, 23)]
    chromosomes.extend(['chrX', 'chrY', 'chrM'])
    vcfs_by_chrom = dict()
    tabix = get_system_path(cnf, 'tabix')
    for chrom in chromosomes:
        vcf_fpath = join(cnf.work_dir, str(chrom) + '.vcf')
        cmdline = '{tabix} -h {combined_vcf_fpath} {chrom} > {vcf_fpath}'.format(**locals())
        call(cnf, cmdline)
        if verify_file(vcf_fpath):
            vcfs_by_chrom[chrom] = vcf_fpath

    output_dirpath = join(cnf.output_dir, 'combined_bams', cnf.project_name)
    safe_mkdir(output_dirpath)
    not_submitted_chroms = vcfs_by_chrom.keys()
    sample_names = ','.join(sample.name for sample in samples)
    sample_bams = ','.join(sample.bam for sample in samples)
    while not_submitted_chroms:
        jobs_to_wait = []
        submitted_chroms = []
        reused_chroms = []

        for chrom, vcf_fpath in vcfs_by_chrom.iteritems():
            if chrom not in not_submitted_chroms:
                continue
            output_fpaths = [join(output_dirpath, chrom.replace('chr', '') + '-' + sample.name.replace('-', '_') +
                                '.bam'.format(**locals())) for sample in samples]
            if cnf.reuse_intermediate and all(verify_file(output_fpath, silent=True) for output_fpath in output_fpaths):
                info('BAM files for ' + chrom + ' chromosome exists, reusing')
                reused_chroms.append(chrom)
                continue
            else:
                if exac_venv_pythonpath:  # to avoid compatibility problems with pysam and tabix
                    cmdline = 'PYTHONPATH= ' + exac_venv_pythonpath + ' ' + get_system_path(cnf,
                                                                            join('tools', 'split_bams_by_variants.py'))
                else:
                    cmdline = get_script_cmdline(cnf, 'python', join('tools', 'split_bams_by_variants.py'), is_critical=True)
                cmdline += ' --chr {chrom} --vcf {vcf_fpath} --samples {sample_names} --bams {sample_bams} ' \
                           '-o {output_dirpath} --work-dir {cnf.work_dir} -g {cnf.genome.name} '.format(**locals())
                if cnf.reuse_intermediate:
                    cmdline += ' --reuse'
                if exac_features_fpath and verify_file(exac_features_fpath):
                    cmdline += ' --features ' + exac_features_fpath
                j = submit_job(cnf, cmdline,  chrom + '_split')
                info()
                submitted_chroms.append(chrom)

                if not j.is_done:
                    jobs_to_wait.append(j)
                if len(jobs_to_wait) >= cnf.threads:
                    break
        if jobs_to_wait:
            info('Submitted ' + str(len(jobs_to_wait)) + ' jobs, waiting...')
            jobs_to_wait = wait_for_jobs(cnf, jobs_to_wait)
        else:
            info('No jobs to submit.')
        not_submitted_chroms = [chrom for chrom in not_submitted_chroms if
                                          chrom not in submitted_chroms and
                                          chrom not in reused_chroms]


def main():
    info(' '.join(sys.argv))
    info()

    description = 'This script prepare data for ExAC browser'
    parser = OptionParser(description=description)
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser)

    parser.add_option('--log-dir', dest='log_dir', default='-')
    parser.add_option('--bed', dest='bed', help='BED file.')
    parser.add_option('-o', dest='output_dir', help='Output directory with ExAC data.')

    cnf, bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths, tags, is_wgs_in_bcbio, is_rnaseq \
        = process_post_bcbio_args(parser)
    exac_dirpath = None
    exac_venv_pythonpath = None
    exac_features_fpath = None

    if is_us():
        if not cnf.genome:
            critical('Usage: ' + __file__ + ' -g hg19 project_bcbio_path [project_bcbio_path] [--bed bed_fpath] [-o output_dir]')
        exac_dirpath = '/ngs/usr/miheenko/git/exac_browser'
        exac_venv_pythonpath = join(exac_dirpath, 'venv_exac', 'bin', 'python')
        exac_features_fpath = os.path.join(exac_dirpath, EXAC_FILES_DIRECTORY, cnf.genome.name, 'all_features.bed.gz')
        cnf.output_dir = join('/ngs/usr/miheenko/git/exac_data', cnf.genome.name)  # temporary dir
    elif not cnf.output_dir:
        critical('Error! Please specify ExAC browser data directory')

    if len(bcbio_project_dirpaths) < 1:
        critical('Usage: ' + __file__ + ' -g hg19 project_bcbio_path [project_bcbio_path] [--bed bed_fpath] [-o output_dir]')

    info()
    info('*' * 70)
    bcbio_structures = []
    for bcbio_project_dirpath, bcbio_cnf, final_dirpath in zip(
            bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths):
        bs = BCBioStructure(cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath)
        bcbio_structures.append(bs)

    if not cnf.project_name:
        if len(bcbio_structures) == 1:
            cnf.project_name = bcbio_structures[0].project_name
        else:
            critical('If you combine multiple BCBIO projects you should specify new project name')
    cnf.caller_name = 'vardict'

    if cnf.output_dir is None:
        critical('Please specify path to ExAC data directory.')
    safe_mkdir(cnf.output_dir)

    cnf.log_dir = join(cnf.output_dir, cnf.project_name + '_log')
    info('log_dirpath: ' + cnf.log_dir)
    safe_mkdir(cnf.log_dir)
    set_up_log(cnf, 'prepare_for_exac', cnf.project_name, cnf.output_dir)

    cnf.work_dir = cnf.work_dir or adjust_path(join(cnf.output_dir, 'work'))
    safe_mkdir(cnf.work_dir)

    samples = [s for bs in bcbio_structures for s in bs.samples]

    vcf_fpath_by_sname = dict()
    for bs in bcbio_structures:
        for sample in bs.samples:
            vcf_fpath, pass_vcf_fpath = convert_txt_to_vcf(cnf, bs, sample, output_dir=cnf.work_dir)
            vcf_fpath_by_sname[sample.name] = vcf_fpath

    info()
    vcf_dirpath = join(cnf.output_dir, 'vardict')
    safe_mkdir(vcf_dirpath)
    combined_vcf_fpath = join(vcf_dirpath, cnf.project_name + '.vcf')
    combine_vcfs(cnf, vcf_fpath_by_sname, combined_vcf_fpath, additional_parameters='--genotypemergeoption UNSORTED')
    split_bam_files_use_grid(cnf, samples, combined_vcf_fpath + '.gz', exac_features_fpath, exac_venv_pythonpath)

    depths_by_pos = get_regions_depth(cnf, samples)
    cov_thresholds = [1, 5, 10, 15, 20, 25, 30, 50, 100]
    chromosomes = ['chr%s' % x for x in range(1, 23)]
    chromosomes.extend(['chrX', 'chrY', 'chrM'])

    info()
    info('Saving coverage')
    project_cov_dirpath = join(cnf.output_dir, 'coverage', cnf.project_name)
    safe_mkdir(project_cov_dirpath)
    for chrom in depths_by_pos.keys():
        if chrom not in chromosomes:
            continue
        coverage_data_fpath = join(project_cov_dirpath, chrom + '.txt')
        if cnf.reuse_intermediate and verify_file(coverage_data_fpath, silent=True):
            continue
        chrom_num = chrom.replace('chr', '')
        with file_transaction(cnf, coverage_data_fpath) as tx:
            with open(tx, 'w') as f:
                fs = ['#chrom', 'pos', 'mean', 'median'] + [str(t) for t in cov_thresholds]
                f.write('\t'.join(fs) + '\n')
                sorted_positions = sorted(depths_by_pos[chrom].keys())
                for pos in sorted_positions:
                    depths = depths_by_pos[chrom][pos]
                    if len(depths) < len(samples):
                        depths += [0] * (len(samples) - len(depths))
                    mean_coverage = mean(depths)
                    median_coverage = median(depths)
                    pcnt_samples_ge_threshold = [mean([1 if d >= t else 0 for d in depths]) for t in cov_thresholds]
                    res_line = chrom_num + '\t' + str(pos) + '\t' + str(mean_coverage) + '\t' + str(median_coverage)
                    for pcnt_samples in pcnt_samples_ge_threshold:
                        res_line += '\t' + str(pcnt_samples)
                    f.write(res_line + '\n')
        bgzip_and_tabix(cnf, coverage_data_fpath, tabix_parameters='-p bed')
    info()
    if exac_dirpath:
        info('Adding project to ExAC database')
        cmdline = 'PYTHONPATH= ' + exac_venv_pythonpath + ' ' + join(exac_dirpath, 'manage.py') + ' ' + 'add_project' + \
                  ' ' + cnf.project_name + ' ' + cnf.genome.name
        call(cnf, cmdline)

    info('Done.')


def log(msg=''):
    sys.stderr.write(msg + '\n')


if __name__ == '__main__':
    main()
