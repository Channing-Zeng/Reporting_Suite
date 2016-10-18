#!/usr/bin/env python
# noinspection PyUnresolvedReferences

import os
import re
import sys
from collections import OrderedDict
from optparse import OptionParser
from os.path import join, basename

from ngs_reporting.utils import get_target_genes

import bcbio_postproc
from source.bcbio.bcbio_structure import BCBioStructure, process_post_bcbio_args
from source.file_utils import safe_mkdir, adjust_path, add_suffix, file_transaction, verify_file, open_gzipsafe
from source.logger import info, warn
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug
from source.variants import vcf_parser as vcf
from source.variants.vcf_processing import bgzip_and_tabix, verify_vcf

from ngs_reporting.clinical_parser import get_record_from_vcf, parse_mutations, get_mutations_fpath_from_bs
from ngs_reporting.combine_reports import get_rejected_mutations_fpaths


filter_descriptions_dict = {
    'not canonical transcript': 'Transcript',
    'PASS=False': 'REJECT',
    'PROTEIN_PROTEIN_CONTACT': 'PROT_PROT',
    'MSI fail': 'MSI',
    'snp in snpeffect export polymorphic': 'PolymorphicSNP',
    'in snpeffect snp': 'SNP_snpeff',
    'common SNP': 'SNP',
    'not act': 'Not_act',
    'not known': 'Not_known',
    'unknown': 'Unknown',
    'act germline': 'Act_germ',
    'act somatic': 'Act_som',
    'SYNONYMOUS': 'Synonymous',
    'not ClnSNP_known': 'Not_ClnSNP_known',
    'not ClnSNP known': 'Not_ClnSNP_known',
    'in filter common snp': 'filterSNP',
    'in filter artifacts': 'Artifact',
    'AF < 0.35': 'f0.35',
    'AF >= 0.5': 'F0.5',
    'clnsig dbSNP': 'clnsig',
    'in snpeff_snp': 'SnpEff_snp',
    'dbSNP': 'dbSNP',
    'in INTRON': 'Intron',
    'no aa ch\g': 'no_aa_chg',
    'SPLICE': 'Splice',
    'UPSTREAM': 'Upstream',
    'DOWNSTREAM': 'Downstream',
    'INTERGENIC': 'Intergenic',
    'INTRAGENIC': 'Intragenic',
    'not UTR /CODON': 'not_UTR_/Codon',
    'NON CODING': 'Non_coding',
    'fclass=NON CODING': 'Non_coding',
    'variants occurs after last known critical amino acid': 'CritAA',
    'blacklist gene': 'Blacklist',
}

filter_patterns_dict = {
    re.compile('depth < (\d+)'): 'd',
    re.compile('VD < (\d+)'): 'v',
    re.compile('AF < ([0-9.]+)'): 'f',
    re.compile('all GMAF > ([0-9.]+)'): 'GMAF',
    re.compile('cohort freq > ([0-9.]+)'): 'CohortFreq',
}


filt_vcf_ending = '.filt.vcf'


def main():
    info(' '.join(sys.argv))
    info()
    description = 'This script converts Vardict TXT file to VCF.'

    parser = OptionParser(description=description, usage='Usage: ' + basename(__file__) + ' [-o Output_directory -c Var_caller_name] Project_directory')
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser)
    parser.add_option('--log-dir', dest='log_dir', default='-')
    parser.add_option('-c', '--caller', dest='caller_name', default='vardict')
    parser.add_option('-o', dest='output_dir', help='Output directory.')

    cnf, bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths, tags, is_wgs_in_bcbio, is_rnaseq \
        = process_post_bcbio_args(parser)

    if not bcbio_project_dirpaths:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    bcbio_structures = []
    for bcbio_project_dirpath, bcbio_cnf, final_dirpath in zip(
            bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths):
        bs = BCBioStructure(cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath)
        bcbio_structures.append(bs)

    cnf.work_dir = cnf.work_dir or adjust_path(join(cnf.output_dir, 'work'))
    safe_mkdir(cnf.work_dir)

    info('')
    info('*' * 70)
    for bs in bcbio_structures:
        for sample in bs.samples:
            if sample.phenotype != 'normal':
                convert_vardict_txts_to_bcbio_vcfs(cnf, bs, sample)


def get_mutation_dicts(cnf, bs, sample, pass_only=False):
    pass_mut_dict = dict()
    reject_mut_dict = dict()
    filter_values = set()

    target_gene_names_chroms = None
    if bs.sv_bed:
        target_gene_names_chroms = get_target_genes(bs.sv_bed)
    pass_mutations_fpath, _ = get_mutations_fpath_from_bs(bs)
    pass_mutations = parse_mutations(cnf.genome.name, sample, dict(), pass_mutations_fpath, '',
                                     target_gene_names_chroms=target_gene_names_chroms)
    for mut in pass_mutations:
        pass_mut_dict[(mut.chrom, mut.pos, mut.transcript)] = mut

    if not pass_only:
        for reject_mutations_fpath in get_rejected_mutations_fpaths(pass_mutations_fpath):
            if verify_file(reject_mutations_fpath, silent=True):
                reject_mutations = parse_mutations(cnf, sample, dict(), reject_mutations_fpath, '')
                for mut in reject_mutations:
                    reject_mut_dict[(mut.chrom, mut.pos, mut.transcript)] = mut
                    filt_values = mut.reason.split(' and ')
                    for val in filt_values:
                        filter_values.add(val)
        for filt_val in filter_values:
            for pattern, description in filter_patterns_dict.iteritems():
                val = pattern.findall(filt_val)
                if val:
                    filter_descriptions_dict[filt_val] = description + val[0]
    return pass_mut_dict, reject_mut_dict, list(filter_values)


def combine_mutations(pass_mut_dict, reject_mut_dict):
    combined_dict = pass_mut_dict.copy()
    combined_dict.update(reject_mut_dict)
    sorted_combined_dict = OrderedDict(sorted(combined_dict.items(), key=lambda x: ([x[0][j] for j in range(len(x[0]))])))
    return sorted_combined_dict


def add_keys_to_header(vcf_reader, filter_values):
    vcf_reader.formats['ADJAF'] = vcf._Format('ADJAF', 1, 'Float', 'Adjusted AF for indels due to local realignment')
    vcf_reader.formats['BIAS'] = vcf._Format('BIAS', 1, 'String', 'Strand Bias Info')
    vcf_reader.formats['HIAF'] = vcf._Format('HIAF', 1, 'Float', 'Allele frequency using only high quality bases')
    vcf_reader.formats['MQ'] = vcf._Format('MQ', 1, 'Float', 'Mean Mapping Quality')
    vcf_reader.formats['NM'] = vcf._Format('NM', 1, 'Float', 'Mean mismatches in reads')
    vcf_reader.formats['ODDRATIO'] = vcf._Format('ODDRATIO', 1, 'Float', 'Strand Bias Oddratio')
    vcf_reader.formats['PMEAN'] = vcf._Format('PMEAN', 1, 'Float', 'Mean position in reads')
    vcf_reader.formats['PSTD'] = vcf._Format('PSTD', 1, 'Float', 'Position STD in reads')
    vcf_reader.formats['QUAL'] = vcf._Format('QUAL', 1, 'Float', 'Mean quality score in reads')
    vcf_reader.formats['QSTD'] = vcf._Format('QSTD', 1, 'Float', 'Quality score STD in reads')
    vcf_reader.formats['SBF'] = vcf._Format('SBF', 1, 'Float', 'Strand Bias Fisher p-value')
    vcf_reader.formats['SN'] = vcf._Format('SN', 1, 'Float', 'Signal to noise')
    vcf_reader.infos['SOR'] = vcf._Info('SOR', 1, 'Float', 'Odds ratio')
    vcf_reader.infos['SSF'] = vcf._Info('SSF', 1, 'Float', 'P-value')
    vcf_reader.infos['StrongSomatic'] = vcf._Info('StrongSomatic', 0, 'Flag', 'Variant is somatic')
    vcf_reader.infos['Signif'] = vcf._Info('Signif', '.', 'String', 'Significance')
    vcf_reader.infos['Status'] = vcf._Info('Status', '.', 'String', 'Status')
    vcf_reader.infos['Reason'] = vcf._Info('Reason', '.', 'String', 'Reason')
    for filt_val in filter_values:
        filter_id = filter_descriptions_dict[filt_val.replace('_', ' ')] if filt_val in filter_descriptions_dict else filt_val
        filt = vcf._Filter(filter_id, filt_val)
        vcf_reader.filters[filter_id] = filt
    return vcf_reader


def txt_to_vcf():
    pass


def convert_vardict_txts_to_bcbio_vcfs(cnf, bs, sample, output_dir=None, pass_only=False):
    info('')
    info('Preparing data for ' + sample.name)
    anno_filt_vcf_fpath = sample.find_filt_vcf_by_callername(cnf.caller_name)
    if not anno_filt_vcf_fpath:
        return None, None

    if not output_dir:
        output_dir = cnf.output_dir or os.path.dirname(anno_filt_vcf_fpath)
    output_vcf_fpath = join(output_dir, sample.name + '-' + cnf.caller_name + filt_vcf_ending)
    pass_output_vcf_fpath = add_suffix(output_vcf_fpath, 'pass')
    if cnf.reuse_intermediate and verify_vcf(output_vcf_fpath + '.gz') and verify_vcf(pass_output_vcf_fpath + '.gz'):
        info(output_vcf_fpath + '.gz and ' + pass_output_vcf_fpath + '.gz exists, reusing')
        return output_vcf_fpath + '.gz', pass_output_vcf_fpath + '.gz'

    info('Parsing PASS and REJECT mutations...')
    pass_mut_dict, reject_mut_dict, filter_values = get_mutation_dicts(cnf, bs, sample, pass_only=pass_only)
    sorted_mut_dict = combine_mutations(pass_mut_dict, reject_mut_dict)

    info('')
    info('Writing VCFs')
    vcf_reader = vcf.Reader(open_gzipsafe(anno_filt_vcf_fpath, 'r'))
    vcf_reader = add_keys_to_header(vcf_reader, filter_values)
    with file_transaction(cnf.work_dir, output_vcf_fpath) as filt_tx, \
        file_transaction(cnf.work_dir, pass_output_vcf_fpath) as pass_tx:
        vcf_writer = None
        if not pass_only:
            vcf_writer = vcf.Writer(open(filt_tx, 'w'), template=vcf_reader)
        vcf_pass_writer = vcf.Writer(open(pass_tx, 'w'), template=vcf_reader)
        for key, mut in sorted_mut_dict.items():
            record = get_record_from_vcf(vcf_reader, mut)
            if record:
                if key in pass_mut_dict:
                    record.FILTER = ['PASS']
                    if mut.reason:
                        record.INFO['Reason'] = mut.reason.replace(' ', '_')
                elif pass_only:
                    continue
                elif key in reject_mut_dict:
                    if not mut.reason:
                        continue
                    reject_reason_ids = [filter_descriptions_dict[reason] if reason in filter_descriptions_dict else reason
                                         for reason in mut.reason.split(' and ')]
                    record.FILTER = [';'.join(reject_reason_ids)]
                if mut.signif:
                    record.INFO['Signif'] = mut.signif
                if mut.status:
                    record.INFO['Status'] = mut.status
                if vcf_writer:
                    vcf_writer.write_record(record)
                if key in pass_mut_dict:
                    vcf_pass_writer.write_record(record)
            else:
                warn('No record was found in ' + anno_filt_vcf_fpath + ' for mutation ' + str(mut))

    output_gzipped_vcf_fpath = None
    if vcf_writer:
        vcf_writer.close()
        output_gzipped_vcf_fpath = bgzip_and_tabix(cnf, output_vcf_fpath)
        info('VCF file for vardict.txt is saved to ' + output_gzipped_vcf_fpath)
    vcf_pass_writer.close()
    output_gzipped_pass_vcf_fpath = bgzip_and_tabix(cnf, pass_output_vcf_fpath)
    info('VCF file for vardict.PASS.txt is saved to ' + output_gzipped_pass_vcf_fpath)
    return output_gzipped_vcf_fpath, output_gzipped_pass_vcf_fpath


if __name__ == '__main__':
    main()
