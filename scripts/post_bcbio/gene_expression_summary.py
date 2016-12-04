#!/usr/bin/env python
# noinspection PyUnresolvedReferences
from os.path import isfile, join

import bcbio_postproc

import sys
from traceback import format_exc

import source
from source.bcbio.bcbio_structure import BCBioStructure, bcbio_summary_script_proc_params
from source.logger import info, step_greetings, err, critical
from source.rnaseq.gene_expression import make_gene_expression_heatmaps, annotate_gene_counts
from source.file_utils import verify_file, open_gzipsafe, safe_mkdir


def _rm_quotes(l):
    return l[1:-1]


def _get_gene_transcripts_id(cnf):
    genes_dict = dict()
    transcripts_dict = dict()

    if not cnf.genome.all_transcripts:
        critical('File with transcripts and genes ID ' + cnf.genome.name + ' was not found! Heatmaps cannot be created.')
    if not verify_file(cnf.genome.all_transcripts):
        critical('File with transcripts and genes ID ' + cnf.genome.name + ' at ' + cnf.genome.all_transcripts + ' was not found! Heatmaps cannot be created.')

    info('Getting transcripts ID and genes ID from ' + cnf.genome.all_transcripts)

    with open_gzipsafe(cnf.genome.all_transcripts) as f:
        for i, l in enumerate(f):
            if l.startswith('#'):
                continue
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
    cnf, bcbio_structure = bcbio_summary_script_proc_params(
            'expression', BCBioStructure.expression_dir)

    step_greetings('Gene expression heatmaps summary for all samples')
    report_caption_names = ['Gene counts', 'Exon counts', 'Gene TPM', 'Isoform TPM']
    genes_dict, transcripts_dict = _get_gene_transcripts_id(cnf)
    for counts_fname, report_caption_name in zip(bcbio_structure.counts_names, report_caption_names):
        counts_fpath = join(bcbio_structure.expression_dirpath, counts_fname)
        if not verify_file(counts_fpath, silent=True):
            raw_counts_fpath = join(bcbio_structure.expression_dirpath, 'raw', 'combined.' + counts_fname)
            info('Annotating ' + report_caption_name + ' from ' + raw_counts_fpath)
            annotate_gene_counts(cnf, raw_counts_fpath, counts_fpath, genes_dict)
        verify_file(counts_fpath, is_critical=True, description=counts_fname)

        isoforms_found = counts_fname == 'isoform.sf.tpm' and counts_fpath
        used_dict = transcripts_dict if isoforms_found else genes_dict
        report_fpath = join(safe_mkdir(join(bcbio_structure.expression_dirpath, 'html')), counts_fname + '.html')

        make_gene_expression_heatmaps(cnf, bcbio_structure, counts_fpath, used_dict, report_fpath,
                                      report_caption_name, keep_gene_names=isoforms_found)
    info('Done')


if __name__ == '__main__':
    main()
