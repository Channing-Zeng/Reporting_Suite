from collections import defaultdict
import os
import sys
from os.path import join
from optparse import OptionParser

from source.bcbio.bcbio_structure import BCBioStructure, process_post_bcbio_args
from source.bcbio.bcbio_runner import BCBioRunner
from source.config import defaults
from source.logger import info
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug, check_system_resources, set_up_log
from source.file_utils import safe_mkdir, adjust_path
from source.targetcov import summarize_targetcov
from source.variants import summarize_qc


def main():
    info(' '.join(sys.argv))
    info()
    description = 'This script runs reporting suite on the bcbio final directory.'

    parser = OptionParser(description=description)
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser)

    parser.add_option('--load-mongo', '--mongo-loader', dest='load_mongo', action='store_true', default=defaults['load_mongo'], help='Load to Mongo DB')
    parser.add_option('--datahub-path', dest='datahub_path', help='DataHub directory path to upload final MAFs and CNV (can be remote).')
    parser.add_option('--email', dest='email', help='E-mail address to send notifications on errors and finished jobs.')
    parser.add_option('--reannotate', dest='reannotate', action='store_true', default=False, help='Re-annotate BED file with gene names')
    parser.add_option('--dedup', dest='dedup', action='store_true', default=False, help='Count duplicates in coverage metrics')
    parser.add_option('--seq2c-opts', dest='seq2c_opts', help='Options for the final lr2gene.pl script.')
    parser.add_option('--seq2c-controls', dest='seq2c_controls', help='Additional controls for Seq2C.')
    parser.add_option('--deep-seq', dest='deep_seq', action='store_true', default=False, help='Use run_info_DeepSeq.yaml')
    parser.add_option('--only-summary', dest='only_summary', action='store_true', default=False, help='Only generate project-level report')
    parser.add_option('--jira', dest='jira', help='JIRA case path')
    parser.add_option('--bed', '--capture', '--amplicons', dest='bed', help='BED file to run targetSeq and Seq2C analysis on.')
    parser.add_option('--exons', '--exome', dest='exons', help='Exons BED file to make targetSeq exon/amplicon regions reports.')
    parser.add_option('-f', '--freq', '--min-freq', dest='min_freq', type='float', help='Minimum allele frequency for the filtering. Default %f' % defaults['default_min_freq'])
    parser.add_option('-o', dest='output_dir', help='Output directory for report combining.')

    cnf, bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths, tags = process_post_bcbio_args(parser)

    cnf_project_name = cnf.project_name
    if len(bcbio_project_dirpaths) > 1:
        cnf.project_name = None

    info()
    info('*' * 70)
    bcbio_structures = []
    for bcbio_project_dirpath, bcbio_cnf, final_dirpath in zip(bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths):
        bs = BCBioStructure(cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath)
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
            cnf.project_name = '_'.join([bs.project_name for bs in bcbio_structures])

        if cnf.output_dir is None:
            cnf.output_dir = join(os.getcwd(), cnf.project_name)

        safe_mkdir(cnf.output_dir)

        cnf.log_dir = join(cnf.output_dir, 'log')
        info('log_dirpath: ' + cnf.log_dir)
        safe_mkdir(cnf.log_dir)
        set_up_log(cnf, 'miltiple_projects', cnf.project_name, cnf.output_dir)

        cnf.work_dir = cnf.work_dir or adjust_path(join(cnf.output_dir, 'work'))
        safe_mkdir(cnf.work_dir)

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

    summarize_qc.make_summary_reports(cnf, 1, output_dir, callers, samples,
         jsons_by_sample_by_caller, htmls_by_sample_by_caller, tag_by_sample,
         varqc_name=varqc_name, caption=caption)


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

    reuse = cnf.reuse_intermediate
    cnf.reuse_intermediate = True

    combine_varqc(cnf, bcbio_structures, tag_by_sample,
        varqc_dirname=BCBioStructure.varqc_dir,
        varqc_name=BCBioStructure.varqc_name,
        caption='Variant QC')
    combine_varqc(cnf, bcbio_structures, tag_by_sample,
        varqc_dirname=BCBioStructure.varqc_after_dir,
        varqc_name=BCBioStructure.varqc_after_name,
        caption='Variant QC after filtering')
    combine_targqc(cnf, bcbio_structures, tag_by_sample)

    cnf.reuse_intermediate = reuse

