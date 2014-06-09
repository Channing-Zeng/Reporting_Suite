#!/usr/bin/env python

from __future__ import print_function
import sys
import os
from source.runner import run_all
from source.targetcov.cov import bedcoverage_hist_stats, run_header_report, intersect_bed, sort_bed, \
    run_region_cov_report

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from os.path import join, expanduser, splitext, basename, isdir, abspath


#downlad hg19.genome
#https://github.com/arq5x/bedtools/tree/master/genomes

#TODO
# log file
# take folder name as a sample name (first column on the report)
from shutil import rmtree
from source.main import read_opts_and_cnfs, check_system_resources, load_genome_resources
from source.utils import verify_file, critical, step_greetings, rmtx, info, make_tmpdir, err


REPORT_TYPES = 'summary,genes'


def main(args):
    required = ['bam', 'bed']
    optional = []
    genome_resources = ['chr_lengths', 'genes', 'exons']

    cnf, options = read_opts_and_cnfs(
        extra_opts=[
            (['--bam'], 'align.bam', {
                'dest': 'bam',
                'help': 'used to generate some annotations by GATK'}),

            (['--bed', '--capture', '--amplicons'], 'amplicons.bed', {
                'dest': 'bed',
                'help': 'capture panel/amplicons'}),

            (['--padding'], '250', {
                'dest': 'padding',
                'help': 'input regions will be extended by this value in both directions',
                'default': 250}),

            (['--reports'], '', {
                'dest': 'reports',
                'help': 'default: --reports ' + REPORT_TYPES,
                'default': REPORT_TYPES}),
        ],
        required_inputs=required)

    check_system_resources(cnf, ['samtools', 'bedtools'])
    load_genome_resources(cnf, genome_resources)

    read_samples_info(cnf, options)

    run_all(cnf, None, required + genome_resources, [], process_one, None, None)


def read_samples_info(cnf, options):
    info('')
    info('Processing input details...')

    cnf['bed'] = options.get('bed') or cnf.get('bed')
    cnf['bam'] = options.get('bam') or cnf.get('bam')
    cnf['genes'] = options.get('genes') or cnf.get('genes') or cnf['genome'].get('genes')
    cnf['exons'] = options.get('exons') or cnf.get('exons') or cnf['genome'].get('exons')
    cnf['chr_lengths'] = cnf.get('chr_lengths') or cnf['genome'].get('chr_lengths')

    if not cnf['genes']:
        critical(cnf.get('log'), 'Specify sorted genes bed file in system info or in run info.')
    if not cnf['exons']:
        critical(cnf.get('log'), 'Specify sorted exons bed file in system info or in run info.')
    if not cnf['chr_lengths']:
        critical(cnf.get('log'), 'Specify chromosome lengths for the genome in system info or in run info.')
    if not cnf['bam']:
        critical(cnf.get('log'), 'Specify bam file by --bam option or in run_config.')
    if not cnf['bed']:
        critical(cnf.get('log'), 'Specify capture file by --capture option or in run_config.')

    info(cnf.get('log'), 'using genes ' + cnf['genes'])
    info(cnf.get('log'), 'using exons ' + cnf['exons'])
    info(cnf.get('log'), 'using chr lengths ' + cnf['chr_lengths'])
    info(cnf.get('log'), 'using bam ' + cnf['bam'])
    info(cnf.get('log'), 'using capture panel/amplicons ' + cnf['bed'])

    genes_bed = expanduser(cnf['genes'])
    exons_bed = expanduser(cnf['exons'])
    chr_len_fpath = expanduser(cnf['chr_lengths'])
    bam = expanduser(cnf['bam'])
    capture_bed = expanduser(cnf['bed'])

    if not verify_file(genes_bed): exit(1)
    if not verify_file(exons_bed): exit(1)
    if not verify_file(chr_len_fpath): exit(1)
    if not verify_file(bam): exit(1)
    if not verify_file(capture_bed): exit(1)

    cnf['reports'] = options['reports']
    cnf['padding'] = options.get('padding') or cnf.get('padding') or 250
    cnf['depth_thresholds'] = cnf.get('depth_thresholds') or [5, 10, 25, 50, 100, 500, 1000, 5000,
                                                              10000, 50000, 100000, 500000]


def process_one(cnf, bam, capture_bed, chr_len_fpath, genes_bed, exons_bed):
    # sample_name, _ = splitext(basename(bam))
    sample_name = os.path.basename(os.path.dirname(bam))

    summary_report_fpath = None
    gene_report_fpath = None

    with make_tmpdir(cnf):
        info(cnf.get('log'), 'Calculation of coverage statistics for the regions in the input BED file...')
        amplicons, combined_region, max_depth, total_bed_size = \
            bedcoverage_hist_stats(cnf, capture_bed, bam)

        if 'summary' in cnf['reports']:
            step_greetings('Target coverage summary report')
            summary_report_fpath = join(cnf['output_dir'], sample_name + '.targetseq.summary.txt')
            run_header_report(
                cnf, summary_report_fpath,
                capture_bed, bam, chr_len_fpath,
                cnf['depth_thresholds'], cnf['padding'],
                combined_region, max_depth, total_bed_size)

        # if 'amplicons' in options['reports']:
        #     step_greetings('Coverage report for the input BED file regions')
        #     amplicons_report_fpath = join(output_dir, sample_name + '.targetseq.details.capture.txt')
        #     run_amplicons_cov_report(cnf, amplicons_report_fpath, sample_name, depth_threshs, amplicons)

        if 'genes' in cnf['reports']:
            if not genes_bed or not exons_bed:
                if cnf['reports'] == 'genes':
                    critical('Error: no genes or exons specified for the genome in system config, '
                             'cannot run per-exon report.')
                else:
                    err(cnf.get('log'), 'Warning: no genes or exons specified for the genome in system config, '
                        'cannot run per-exon report.')
            else:
                # log('Annotating amplicons.')
                # annotate_amplicons(amplicons, genes_bed)

                info(cnf.get('log'), 'Getting the gene regions that intersect with our capture panel.')
                bed = intersect_bed(cnf, genes_bed, capture_bed)
                info(cnf.get('log'), 'Getting the exons of the genes.')
                bed = intersect_bed(cnf, exons_bed, bed)
                info(cnf.get('log'), 'Sorting final exon BED file.')
                bed = sort_bed(cnf, bed)

                info(cnf.get('log'), 'Calculation of coverage statistics for exons of the genes ovelapping with the input regions...')
                exons, _, _, _ = bedcoverage_hist_stats(cnf, bed, bam)
                for exon in exons:
                    exon.gene_name = exon.extra_fields[0]

                gene_report_fpath = join(cnf['output_dir'], sample_name + '.targetseq.details.gene.txt')
                run_region_cov_report(cnf, gene_report_fpath, sample_name, cnf['depth_thresholds'],
                                      amplicons, exons)

        rmtx(cnf['work_dir'])
        rmtx(cnf['output_dir'])

        info(cnf.get('log'), '')
        info(cnf.get('log'), '*' * 70)
        if summary_report_fpath:
            info(cnf.get('log'), 'Summary report: ' + summary_report_fpath)
        if gene_report_fpath:
            info(cnf.get('log'), 'Exons coverage report: ' + gene_report_fpath)


if __name__ == '__main__':
    main(sys.argv)