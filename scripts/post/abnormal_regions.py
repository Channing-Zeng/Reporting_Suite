#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc

import sys
import shutil
from source import BaseSample
from source.bcbio.bcbio_structure import summary_script_proc_params, BCBioStructure
from source.clinical_reporting.clinical_parser import clinical_sample_info_from_cnf, get_ave_coverage, get_key_genes, \
    KeyGene, parse_mutations
from source.file_utils import adjust_path
from source.prepare_args_and_cnf import check_system_resources
from source.main import read_opts_and_cnfs
from source.targetcov.flag_regions import generate_flagged_regions_report
from source.targetcov.summarize_targetcov import _generate_summary_flagged_regions_report
from source.utils import info


def main(args):
    cnf, bcbio_structure = summary_script_proc_params(
        BCBioStructure.targqc_name,
        BCBioStructure.targqc_summary_dir,
        extra_opts=[
            (['--mutations'], dict(
                dest='mutations_fpath',
            )),
        ])

    check_system_resources(
        cnf,
        required=['samtools', 'bedtools'],
        optional=[])
    process_all(cnf, bcbio_structure.samples)

    if not cnf['keep_intermediate']:
        shutil.rmtree(cnf['work_dir'])


def process_all(cnf, samples):
    key_gene_by_name = dict()
    for gene_name in get_key_genes(cnf.key_genes):
        key_gene_by_name[gene_name] = KeyGene(gene_name)
    mutations = {}
    for sample in samples:
        mutations[sample.name] = parse_mutations(cnf, sample, key_gene_by_name, cnf.mutations_fpath, for_flagged_report=True)
    _generate_summary_flagged_regions_report(cnf.output_dir, samples, cnf, mutations)
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