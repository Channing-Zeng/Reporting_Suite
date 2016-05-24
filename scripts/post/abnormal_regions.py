#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc

import sys
import shutil
from source import BaseSample, verify_file
from source.bcbio.bcbio_structure import bcbio_summary_script_proc_params, BCBioStructure
from source.clinical_reporting.clinical_parser import parse_mutations, get_key_or_target_bed_genes
from source.file_utils import adjust_system_path
from source.prepare_args_and_cnf import check_system_resources
from source.targetcov.summarize_targetcov import _generate_summary_flagged_regions_report
from source.utils import info


def main(args):
    cnf, bcbio_structure = bcbio_summary_script_proc_params(
        BCBioStructure.targqc_name,
        BCBioStructure.targqc_summary_dir,
        extra_opts=[
            (['--mutations'], dict(
                dest='mutations_fpath',
            )),
            (['--bed', '--capture', '--amplicons'], dict(
                dest='bed',
                help='a BED file for capture panel or amplicons')
             ),
        ])

    check_system_resources(
        cnf,
        required=['samtools', 'bedtools'],
        optional=[])
    process_all(cnf, bcbio_structure)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


def process_all(cnf, bcbio_structure):
    samples = bcbio_structure.samples
    key_gene_by_name, use_custom_panel = get_key_or_target_bed_genes(cnf.bed, verify_file(adjust_system_path(cnf.key_genes), 'key genes'))
    key_or_target_genes = 'target' if use_custom_panel else 'key'
    mutations = {}
    for sample in samples:
        mutations[sample.name] = parse_mutations(cnf, sample, key_gene_by_name, cnf.mutations_fpath, key_or_target_genes,
                                                 for_flagged_report=True)
    _generate_summary_flagged_regions_report(cnf, bcbio_structure, samples, mutations, key_or_target_genes)
    pass
    # read all detail reports
    # normalize
    # extrac low-cov
    # report cov and missed vars


class Sample(BaseSample):
    def __init__(self, name, output_dir, **kwargs):
        BaseSample.__init__(self, name, output_dir, **kwargs)


# def process_one(cnf, output_dir):
#     sample = Sample(cnf.name, output_dir, bam=cnf.bam, bed=cnf.bed)
#     return make_flagged_regions_reports(cnf, output_dir, sample)


def finalize_one(cnf, *abnormal_regions_reports):
    msg = ['Regions with abnormal regions finished for ' + cnf.sample + ':']

    if abnormal_regions_reports:
        msg.append('Abnormal region reports: ')
        info('Abnormal region reports:')
        for rep in abnormal_regions_reports:
            msg.append('  ' + rep)
            info('  ' + rep)

    # send_email('\n'.join(msg))


if __name__ == '__main__':
    main(sys.argv)