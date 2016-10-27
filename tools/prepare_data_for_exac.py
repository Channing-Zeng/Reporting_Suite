#!/usr/bin/env python
import bcbio_postproc
import os
import sys
from optparse import OptionParser

from os.path import join, basename

from source.bcbio.bcbio_structure import BCBioStructure, process_post_bcbio_args
from source.calling_process import call
from source.file_utils import safe_mkdir, adjust_path, verify_file
from source.logger import critical, info, is_local, warn
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug, set_up_log
from source.qsub_utils import wait_for_jobs, submit_job
from source.targetcov.bam_and_bed_utils import call_sambamba, verify_bam
from source.tools_from_cnf import get_script_cmdline, get_system_path
from source.utils import get_chr_len_fpath, is_us, get_ext_tools_dirname
from source.variants.filtering import combine_vcfs

from ngs_reporting.combine_reports import get_uniq_sample_key
from variant_filtering.txt2vcf_conversion import convert_vardict_txts_to_bcbio_vcfs


chromosomes = ['chr%s' % x for x in range(1, 23)]
chromosomes.extend(['chrX', 'chrY', 'chrM'])
# exac_us_url = 'http://172.18.72.170:5000/'
# exac_code_dir = '/ngs/usr/vlad/exac_browser'
# exac_data_dir = '/ngs/usr/vlad/exac_data'
exac_url = 'http://172.18.72.171:5000/'
exac_code_dir = '/ngs/oncology/exac/exac_browser'
exac_data_dir = '/ngs/oncology/exac/exac_data'
if is_local():
    exac_url = 'http://localhost:5000'
    exac_code_dir = '/Users/vlad/vagrant/exac_browser'
    exac_data_dir = '/Users/vlad/vagrant/exac_data'


def get_exac_us_url(genome, project_name):
    return exac_url + genome.split('-')[0] + '/' + project_name + '/'


def _submit_region_cov(cnf, work_dir, chrom, bam_fpaths, sample_names, output_dirpath, chr_len_fpath):
    if not bam_fpaths or not sample_names:
        return None

    cmdline = get_script_cmdline(cnf, 'python', join('tools', 'get_region_coverage.py'), is_critical=True)
    cmdline += (' --chr ' + chrom + ' --bams ' + bam_fpaths + ' --samples ' + sample_names +
                ' -o ' + output_dirpath + ' -g ' + chr_len_fpath + ' ' +
                ' --work-dir ' + work_dir)
    if cnf.bed:
         cmdline += ' --bed ' + cnf.bed
    return submit_job(cnf, cmdline, chrom + '_coverage_' + ('project' if (',' in sample_names) else sample_names))


def calculate_coverage_use_grid(cnf, samples, output_dirpath):
    assert len(samples) > 0

    sambamba = get_system_path(cnf, join(get_ext_tools_dirname(), 'sambamba'), is_critical=True)

    chr_len_fpath = get_chr_len_fpath(cnf)

    for chrom in chromosomes:
        info('Processing chromosome ' + chrom)
        avg_cov_output_fpath = join(output_dirpath, chrom + '.txt.gz')

        sliced_bam_by_sample = dict()

        sample_names = ','.join(sample.name for sample in samples)
        chrom_bams = []
        jobs_to_wait = []

        for sample in samples:
            assert verify_bam(sample.bam), 'BAM for ' + sample.name
            output_bam_fpath = join(cnf.work_dir, basename(sample.name) + '_' + str(chrom) + '.bam')
            cmdline = '{sambamba} slice {sample.bam} {chrom}'.format(**locals())
            call(cnf, cmdline, output_fpath=output_bam_fpath)
            if verify_file(output_bam_fpath):
                chrom_bams.append(output_bam_fpath)
                sliced_bam_by_sample[sample.name] = output_bam_fpath

        bam_fpaths = ','.join(chrom_bams)

        if cnf.reuse_intermediate and verify_file(avg_cov_output_fpath, silent=True):
            info(avg_cov_output_fpath + ' exists, reusing')
        else:
            j = _submit_region_cov(cnf, cnf.work_dir, chrom, bam_fpaths, sample_names, output_dirpath, chr_len_fpath)
            if j and not j.is_done:
                jobs_to_wait.append(j)
            info()

        for sample in samples:
            info('Coverage calculation for ' + sample.name + ', ' + chrom)
            sample_output_dirpath = join(output_dirpath, sample.name)
            safe_mkdir(sample_output_dirpath)
            output_fpath = join(sample_output_dirpath, chrom + '.txt.gz')
            if cnf.reuse_intermediate and verify_file(output_fpath, silent=True):
                info(output_fpath + ' exists, reusing')
                continue

            sample_bam_fpath = sliced_bam_by_sample[sample.name]

            work_dir = safe_mkdir(join(cnf.work_dir, 'get_region_coverage_' + sample.name + '_' + chrom))
            j = _submit_region_cov(cnf, work_dir, chrom, sample_bam_fpath, sample.name, sample_output_dirpath, chr_len_fpath)
            if j and not j.is_done:
                jobs_to_wait.append(j)
            info()

        if len(jobs_to_wait) >= cnf.threads:
            info('Submitted ' + str(len(jobs_to_wait)) + ' jobs, waiting...')
            jobs_to_wait = wait_for_jobs(cnf, jobs_to_wait)
        else:
            info('No jobs to submit.')

    for chrom in chromosomes:
        avg_cov_output_fpath = join(output_dirpath, chrom + '.txt.gz')
        if not verify_file(avg_cov_output_fpath, silent=True):
            warn('Project has no coverage at chromosome ' + chrom)
        for sample in samples:
            sample_output_dirpath = join(output_dirpath, sample.name)
            output_fpath = join(sample_output_dirpath, chrom + '.txt.gz')
            if not verify_file(output_fpath, silent=True):
                warn(sample.name + ' has no coverage at chromosome ' + chrom)


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
                                      use_grid=True, command_name='sort', sample_name=sample.name,
                                      stdout_to_outputfile=False)
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


def split_bam_files_use_grid(cnf, samples, combined_vcf_fpath, exac_features_fpath):
    samples = dedup_and_sort_bams_use_grid(cnf, samples, do_sort=False)
    samples = dedup_and_sort_bams_use_grid(cnf, samples, do_sort=True)

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
                # if exac_venv_pythonpath:  # to avoid compatibility problems with pysam and tabix
                #     cmdline = exac_venv_pythonpath + ' ' + get_system_path(cnf,
                #                                                             join('tools', 'split_bams_by_variants.py'))
                # else:
                cmdline = get_script_cmdline(cnf, 'python', join('tools', 'split_bams_by_variants.py'), is_critical=True)
                cmdline += (' --chr {chrom} --vcf {vcf_fpath} --samples {sample_names} ' +
                            '--bams {sample_bams} -o {output_dirpath} --work-dir {cnf.work_dir} ' +
                            '-g {cnf.genome.name} ').format(**locals())
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


def get_exac_dir(cnf):
    exac_dir = join(exac_data_dir, cnf.genome.name)  # temporary dir
    return exac_dir


def add_project_to_exac(cnf):
    info('Adding project to ExAC database')
    exac_venv_pythonpath = join(exac_code_dir, 'venv_exac', 'bin', 'python')
    if is_local():
        exac_venv_pythonpath = 'python'
    cmdline = exac_venv_pythonpath + ' ' + join(exac_code_dir, 'manage.py') + ' ' + 'add_project' + \
              ' ' + cnf.project_name + ' ' + cnf.genome.name
    call(cnf, cmdline)


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

    if not cnf.genome:
        critical('Usage: ' + __file__ + ' -g hg19 project_bcbio_path [project_bcbio_path] [--bed bed_fpath] [-o output_dir]')
    exac_venv_pythonpath = None
    if is_us():
        exac_venv_pythonpath = join(exac_code_dir, 'venv_exac', 'bin', 'python')
    cnf.output_dir = get_exac_dir(cnf)
    # if not cnf.output_dir:
    #     critical('Error! Please specify ExAC browser data directory')

    if len(bcbio_project_dirpaths) < 1:
        critical('Usage: ' + __file__ + ' -g hg19 project_bcbio_path [project_bcbio_path] [--bed bed_fpath] [-o output_dir]')

    info()
    info('*' * 70)
    bcbio_structures = []
    project_name = cnf.project_name
    cnf.project_name = None
    for bcbio_project_dirpath, bcbio_cnf, final_dirpath in zip(
            bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths):
        bs = BCBioStructure(cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath)
        bcbio_structures.append(bs)

    cnf.project_name = project_name
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

    cnf.work_dir = cnf.work_dir or adjust_path(join(cnf.output_dir, 'work', cnf.project_name))
    safe_mkdir(cnf.work_dir)

    sample_names = [s.name for bs in bcbio_structures for s in bs.samples]
    samples = []

    info()
    info('Preparing variants data')
    vcf_fpath_by_sname = dict()
    for bs in bcbio_structures:
        for sample in bs.samples:
            sample.name = get_uniq_sample_key(bs.project_name, sample)
            samples.append(sample)
            vcf_fpath, pass_vcf_fpath = convert_vardict_txts_to_bcbio_vcfs(
                    cnf.work_dir, cnf.genome.name, bs, sample, cnf.caller_name,
                    output_dirpath=cnf.work_dir, pass_only=False)
            if vcf_fpath:
                vcf_fpath_by_sname[sample.name] = vcf_fpath
            elif pass_vcf_fpath:
                vcf_fpath_by_sname[sample.name] = pass_vcf_fpath

    if not vcf_fpath_by_sname:
        info('No VCFs found, skipping preparing variants')
    else:
        info()
        variants_dirpath = join(cnf.output_dir, 'vardict')
        safe_mkdir(variants_dirpath)
        combined_vcf_fpath = join(variants_dirpath, cnf.project_name + '.vcf')
        combined_vcf_fpath = combine_vcfs(cnf, vcf_fpath_by_sname, combined_vcf_fpath, additional_parameters='--genotypemergeoption UNSORTED')

        info()
        info('Creating BAM files for IGV')
        exac_features_fpath = os.path.join(exac_data_dir, cnf.genome.name, 'all_features.bed.gz')
        split_bam_files_use_grid(cnf, samples, combined_vcf_fpath, exac_features_fpath)

    info()
    info('Saving coverage')
    project_cov_dirpath = join(cnf.output_dir, 'coverage', cnf.project_name)
    safe_mkdir(project_cov_dirpath)
    calculate_coverage_use_grid(cnf, samples, project_cov_dirpath)

    info()
    add_project_to_exac(cnf)
    info('Done.')


def log(msg=''):
    sys.stderr.write(msg + '\n')


if __name__ == '__main__':
    main()
