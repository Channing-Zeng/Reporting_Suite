from collections import defaultdict
import os
import sys
import shutil
import time
import datetime
from genericpath import exists
from os.path import join
from optparse import OptionParser, SUPPRESS_HELP
from yaml import dump as save_yaml


import source
from source.bcbio.bcbio_structure import BCBioStructure, process_post_bcbio_args
from source.bcbio.bcbio_runner import BCBioRunner
from source.config import defaults
from source.logger import info, err
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug, check_system_resources, set_up_log
from source.file_utils import safe_mkdir, adjust_path, safe_symlink_to, add_suffix, safe_remove
from source.targetcov import summarize_targetcov
from source.variants import summarize_qc
from source.variants.filtering import combine_results


def main():
    description = 'This script runs reporting suite on the bcbio final directory.'

    parser = OptionParser(description=description)
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser)

    parser.add_option('--load-mongo', '--mongo-loader', dest='load_mongo', action='store_true', default=defaults['load_mongo'], help='Load to Mongo DB')
    parser.add_option('--datahub-path', dest='datahub_path', help='DataHub directory path to upload final MAFs and CNV (can be remote).')
    parser.add_option('--email', dest='email', help='E-mail address to send notifications on errors and finished jobs.')
    parser.add_option('--reannotate', dest='reannotate', action='store_true', default=False, help='Re-annotate BED file with gene names')
    parser.add_option('--extended', dest='extended', action='store_true', default=False, help='Count flagged regions and missed variants')
    parser.add_option('--dedup', dest='dedup', action='store_true', default=False, help='Count duplicates in coverage metrics')
    parser.add_option('--seq2c-opts', dest='seq2c_opts', help='Options for the final lr2gene.pl script.')
    parser.add_option('--seq2c-controls', dest='seq2c_controls', help='Additional controls for Seq2C.')
    parser.add_option('--deep-seq', dest='deep_seq', action='store_true', default=False, help='Use run_info_DeepSeq.yaml')
    parser.add_option('--wgs', dest='is_wgs', action='store_true', default=None, help='Ignore sv_regions and run as WGS')
    parser.add_option('--only-summary', dest='only_summary', action='store_true', default=False, help='Only generate project-level report')
    parser.add_option('--jira', dest='jira', help='JIRA case path')
    parser.add_option('--bed', '--capture', '--amplicons', dest='bed', help='BED file to run targetSeq and Seq2C analysis on.')
    parser.add_option('--exons', '--exome', dest='exons', help='Exons BED file to make targetSeq exon/amplicon regions reports.')
    parser.add_option('--no-prep-bed', dest='prep_bed', help='do not fix input beds and exons', action='store_false', default=True)
    parser.add_option('--no-dedup', dest='no_dedup', action='store_true', help=SUPPRESS_HELP)
    parser.add_option('-f', '--freq', '--min-freq', dest='min_freq', type='float', help='Minimum allele frequency for the filtering.')
    parser.add_option('-o', dest='output_dir', help='Output directory for report combining.')
    parser.add_option('--no-bam2bigwig', dest='no_bam2bigwig', action='store_true', default=False, help=SUPPRESS_HELP)

    cnf, bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths, tags, is_wgs_in_bcbio, is_rnaseq \
        = process_post_bcbio_args(parser)
    is_wgs = cnf.is_wgs = cnf.is_wgs or is_wgs_in_bcbio

    cnf.run_date = time.localtime()
    cnf_project_name = cnf.project_name
    if len(bcbio_project_dirpaths) > 1:
        cnf.project_name = None

    info()
    info('*' * 70)
    bcbio_structures = []
    for bcbio_project_dirpath, bcbio_cnf, final_dirpath in zip(bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths):
        bs = BCBioStructure(cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath,
                            is_wgs=is_wgs, is_rnaseq=is_rnaseq)
        bcbio_structures.append(bs)

    # Post-processing one bcbio project as usually
    if len(bcbio_structures) == 1:
        if cnf.min_freq is not None:
            info('Min freq for filtering is %f' % cnf.min_freq)

        if cnf.steps and cnf.load_mongo and 'MongoLoader' not in cnf.steps:
            cnf.steps.append('MongoLoader')

        check_system_resources(cnf, required=['qsub'], optional='transcripts_fpath')

        bcbio_structure = bcbio_structures[0]

        bcbio_runner = BCBioRunner(cnf, bcbio_structure, cnf.bcbio_cnf)
        bcbio_runner.post_jobs()

    # Special case: multiple projects in input. No post-processing them, but rather combining summary reports together.
    elif len(bcbio_structures) > 1:
        if cnf_project_name:
            cnf.project_name = cnf_project_name
        else:
            # cnf.project_name = '_'.join([bs.project_name for bs in bcbio_structures])
            cnf.project_name = 'Combined_project'

        if cnf.output_dir is None:
            cnf.output_dir = join(os.getcwd(), cnf.project_name)

        safe_mkdir(cnf.output_dir)

        cnf.log_dir = join(cnf.output_dir, 'log')
        info('log_dirpath: ' + cnf.log_dir)
        safe_mkdir(cnf.log_dir)
        set_up_log(cnf, 'miltiple_projects', cnf.project_name, cnf.output_dir)

        cnf.work_dir = adjust_path(join(cnf.output_dir, 'work'))
        safe_mkdir(cnf.work_dir)
        safe_mkdir(adjust_path(join(cnf.output_dir, 'config')))

        combine_projects(cnf, bcbio_structures, tags)


def combine_varqc(cnf, bcbio_structures, tag_by_sample, varqc_dirname, varqc_name, caption):
    callers = []
    samples = []

    for bc in bcbio_structures:
        for vc in bc.variant_callers.values():
            if vc.name not in [c.name for c in callers]:
                callers.append(vc)

    jsons_by_sample_by_caller = defaultdict(dict)
    htmls_by_sample_by_caller = defaultdict(dict)
    for bc in bcbio_structures:
        for vc in bc.variant_callers.values():
            fpath_by_sample = vc.find_fpaths_by_sample(varqc_dirname, varqc_name, 'json', bc.final_dirpath)
            for sname, fpath in fpath_by_sample.items():
                jsons_by_sample_by_caller[vc.name][sname] = fpath
            fpath_by_sample = vc.find_fpaths_by_sample(varqc_dirname, varqc_name, 'html', bc.final_dirpath)
            for sname, fpath in fpath_by_sample.items():
                htmls_by_sample_by_caller[vc.name][sname] = fpath
            samples.extend(vc.samples)

    output_dir = join(cnf.output_dir, varqc_dirname)
    safe_mkdir(output_dir)

    if jsons_by_sample_by_caller and htmls_by_sample_by_caller:
        summarize_qc.make_summary_reports(cnf, 1, output_dir, callers, samples,
             jsons_by_sample_by_caller, htmls_by_sample_by_caller, tag_by_sample,
             varqc_name=varqc_name, caption=caption)
    else:
        err('Not JSON and HTML found, cannot generate summary reports.')


def combine_targqc(cnf, bcbio_structures, tag_by_sample):
    samples = [s for bs in bcbio_structures for s in bs.samples]
    output_dir = join(cnf.output_dir, BCBioStructure.targqc_summary_dir)
    safe_mkdir(output_dir)

    summarize_targetcov.summarize_targqc(cnf, cnf.threads or len(samples),
         output_dir, samples, tag_by_sample=tag_by_sample)


def combine_projects(cnf, bcbio_structures, tags=None):
    tag_by_sample = dict()
    if tags:
        for bs, tag in zip(bcbio_structures, tags):
            for s in bs.samples:
                tag_by_sample[s.name] = tag or bs.project_name
    # else:
    #     for bs in bcbio_structures:
    #         for s in bs.sampels:
    #             tag_by_sample[s.name] = bs.project_name

    final_dirpath = adjust_path(join(cnf.output_dir, 'final'))
    safe_mkdir(final_dirpath)

    merged_bcbio_cnf = merge_bcbio_yamls(cnf, bcbio_structures)

    samples = [s for bs in bcbio_structures for s in bs.samples]
    dirs_to_reprocess = [source.clinreport_dir, BCBioStructure.var_dir, source.varannotate_name, source.varfilter_name]
    for s in samples:
        sample_dir = join(final_dirpath, s.name)
        sample_var_dirpath = join(sample_dir, BCBioStructure.var_dir)
        safe_mkdir(sample_var_dirpath)
        for file_or_dir in os.listdir(s.dirpath):
            if file_or_dir not in dirs_to_reprocess:
                safe_symlink_to(join(s.dirpath, file_or_dir), sample_dir)
        for file in os.listdir(s.var_dirpath):
            safe_symlink_to(join(s.var_dirpath, file), sample_var_dirpath)

    merged_date_dir = join(final_dirpath, merged_bcbio_cnf['fc_date'] + '_' + merged_bcbio_cnf['fc_name'])
    merged_bs_var_dirpath = join(merged_date_dir, BCBioStructure.var_dir)
    merged_bs_raw_var_dirpath = join(merged_bs_var_dirpath, 'raw')
    safe_mkdir(merged_bs_raw_var_dirpath)
    for bs in bcbio_structures:
        for file in os.listdir(bs.raw_var_dirpath):
            safe_symlink_to(join(bs.raw_var_dirpath, file), merged_bs_raw_var_dirpath)

    variants_fpaths = []
    vardict_txt_fname = source.mut_fname_template.format(caller_name='vardict')
    variants_fpath = join(merged_bs_var_dirpath, vardict_txt_fname)
    pass_variants_fpath = add_suffix(variants_fpath, source.mut_pass_suffix)
    reject_variants_fpath = add_suffix(variants_fpath, source.mut_reject_suffix)

    cnf.reuse_intermediate = False
    cnf.steps = ['Variants']

    for bs_i, bs in enumerate(bcbio_structures):  # re-filtering, perform cohort-based filtering only within sub-projects
        correct_bs = BCBioStructure(cnf, cnf.output_dir, bs.bcbio_cnf, final_dirpath)
        bcbio_runner = BCBioRunner(cnf, correct_bs, bs.bcbio_cnf)
        bcbio_runner.post_jobs()
        bs_raw_variants_fpath = add_suffix(variants_fpath, str(bs_i))
        pass_bs_variants_fpath = add_suffix(bs_raw_variants_fpath, source.mut_pass_suffix)
        reject_bs_variants_fpath = add_suffix(bs_raw_variants_fpath, source.mut_reject_suffix)
        shutil.move(variants_fpath, bs_raw_variants_fpath)
        shutil.move(pass_variants_fpath, pass_bs_variants_fpath)
        shutil.move(reject_variants_fpath, reject_bs_variants_fpath)
        variants_fpaths.append(bs_raw_variants_fpath)

    merged_bs = BCBioStructure(cnf, cnf.output_dir, merged_bcbio_cnf, final_dirpath)
    merged_samples = [s for s in merged_bs.samples]

    cnf.variant_filtering.max_ratio = 1
    combine_results(cnf, merged_samples, variants_fpaths, variants_fpath, pass_variants_fpath=pass_variants_fpath)
    for variants_fpath in variants_fpaths:
        safe_remove(variants_fpath)
        pass_fpath = add_suffix(variants_fpath, source.mut_pass_suffix)
        safe_remove(pass_fpath)
        reject_fpath = add_suffix(variants_fpath, source.mut_reject_suffix)
        safe_remove(reject_fpath)

    cnf.reuse_intermediate = True
    cnf.steps = ['Seq2C', 'Summary']
    BCBioRunner(cnf, merged_bs, merged_bs.bcbio_cnf).post_jobs()


def merge_bcbio_yamls(cnf, bcbio_structures):
    today_date = datetime.datetime.now()
    today_bcbio_date = today_date.strftime("%Y-%m-%d")
    safe_mkdir(today_bcbio_date)

    bcbio_cnfs = [bs.bcbio_cnf for bs in bcbio_structures]
    merged_yaml_fpath = join(cnf.output_dir, 'config', 'bcbio.yaml')
    merged_bcbio_cnf = dict()
    merged_bcbio_cnf['fc_date'] = today_bcbio_date
    merged_bcbio_cnf['fc_name'] = 'bcbio'
    merged_bcbio_cnf['upload'] = bcbio_cnfs[0]['upload']
    merged_bcbio_cnf['details'] = []
    for bs_cnf in bcbio_cnfs:
        bs_cnf['fc_date'] = today_bcbio_date
        bs_cnf['fc_name'] = 'bcbio'
        merged_bcbio_cnf['details'].extend(bs_cnf['details'])
    with open(merged_yaml_fpath, 'w') as yaml_file:
        yaml_file.write(save_yaml(merged_bcbio_cnf))
    return merged_bcbio_cnf

