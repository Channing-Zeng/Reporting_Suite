#!/usr/bin/env python
# noinspection PyUnresolvedReferences
from os.path import isfile

import bcbio_postproc

import sys
from traceback import format_exc

import source
from source.bcbio.bcbio_structure import BCBioStructure, bcbio_summary_script_proc_params
from source.logger import info, step_greetings, err
from source.rnaseq.gene_expression import make_gene_expression_heatmaps
from tools.bed_processing.make_exons import _rm_quotes


def get_gene_transcripts_id(cnf):
    genes_dict = dict()
    transcripts_dict = dict()
    if not cnf.genome.all_transcripts or not isfile(cnf.genome.all_transcripts):
        err('File with transcripts and genes ID was not found! Heatmaps cannot be created.')
    info('Getting transcripts ID and genes ID from ' + cnf.genome.all_transcripts)

    with open(cnf.genome.all_transcripts) as f:
        for i, l in enumerate(f):
            chrom, _, feature, start, end, _, strand, _, props_line = l.replace('\n', '').split('\t')
            if feature != 'transcript':
                continue
            try:
                _prop_dict = dict((t.strip().split(' ')[0], ' '.join(t.strip().split(' ')[1:]))
                                  for t in props_line.split(';') if t.strip())
            except ValueError:
                sys.stderr.write(format_exc())
                sys.stderr.write(l)

            gene_symbol = _rm_quotes(_prop_dict['gene_name'])
            gene_id = _rm_quotes(_prop_dict['gene_id'])
            transcript_id = _rm_quotes(_prop_dict['transcript_id'])
            #gene = Gene(gene_symbol, chrom=chrom, gene_id=gene_id, transcript_id=transcript_id)
            genes_dict[gene_id] = gene_symbol
            transcripts_dict[transcript_id] = gene_symbol
    return genes_dict, transcripts_dict


def main():
    info(' '.join(sys.argv))
    info()
    cnf, bcbio_structure = bcbio_summary_script_proc_params(BCBioStructure.gene_counts_name, BCBioStructure.gene_counts_dir)

    step_greetings('Gene expression heatmaps summary for all samples')
    counts_fpaths = [bcbio_structure.gene_counts_fpath, bcbio_structure.exon_counts_fpath,
                              bcbio_structure.gene_tpm_fpath, bcbio_structure.isoform_tpm_fpath]
    report_fpaths = [bcbio_structure.gene_counts_report_fpath, bcbio_structure.exon_counts_report_fpath,
                              bcbio_structure.gene_tpm_report_fpath, bcbio_structure.isoform_tpm_report_fpath]
    report_caption_names = ['Gene counts', 'Exon counts', 'Gene TPM', 'Isoform TPM']
    genes_dict, transcripts_dict = get_gene_transcripts_id(cnf)
    for counts_fpath, report_fpath, report_caption_name in zip(counts_fpaths, report_fpaths, report_caption_names):
        is_isoforms = counts_fpath == bcbio_structure.isoform_tpm_fpath
        using_dict = transcripts_dict if is_isoforms else genes_dict
        make_gene_expression_heatmaps(cnf, bcbio_structure, counts_fpath, using_dict, report_fpath, report_caption_name,
                                      not_rename_genes=is_isoforms)

    info('Done')


if __name__ == '__main__':
    main()
