#!/usr/bin/env python
import bcbio_postproc

import os
import re
import sys
import pysam
from optparse import OptionParser
from os.path import basename, join
from collections import defaultdict

from source import BaseSample
from source.calling_process import call
from source.config import Config
from source.file_utils import which, open_gzipsafe, verify_file, safe_mkdir
from source.logger import info, critical, err
from source.prepare_args_and_cnf import determine_sys_cnf, add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug
from source.utils import get_chr_lengths_from_seq


def parse_variants(vcf_fpath):
    variants_by_chrom = defaultdict(list)
    with open_gzipsafe(vcf_fpath) as vcf:
        for line in vcf:
            line = line.strip('\n')
            if line.startswith('##INFO=<ID=ANN'):
                ann_field_names = line.split('Format: ')[-1].strip('">').split('|')
                ann_field_names = [f.strip() for f in ann_field_names]
                ann_field_names[0] = ann_field_names[0].split('\'')[1]
            if line.startswith('#'):
                continue

            fields = line.split('\t')
            info_field = dict([(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', fields[7])])
            annotation_array = info_field['ANN'].split(',') if 'ANN' in info_field else []
            all_annotations = [dict(zip(ann_field_names, x.split('|'))) for x in annotation_array if len(ann_field_names) == len(x.split('|'))]
            coding_annotations = [ann for ann in all_annotations if ann['Feature_ID'].startswith('NM')]

            variant = dict()
            variant['chrom'] = fields[0]
            alt_alleles = fields[4].split(',')

            # different variant for each alt allele
            for i, alt_allele in enumerate(alt_alleles):
                annotations = [ann for ann in coding_annotations if (ann['Allele']) == alt_allele]

                variant['pos'], variant['ref'], variant['alt'] = get_minimal_representation(fields[1], fields[3], alt_allele)
                variant['transcripts'] = set()
                for annotation in annotations:
                    transcript = annotation['Feature_ID'].split('.')[0]
                    variant['transcripts'].add(transcript)

                variant['transcripts'] = list(variant['transcripts'])
                variants_by_chrom[variant['chrom']].append(variant)
    return variants_by_chrom


def split_bams(cnf, samples, vcf_fpath):
    variants_by_chrom = parse_variants(vcf_fpath)
    temp_output_dirpath = join(cnf.work_dir, 'temp')
    safe_mkdir(temp_output_dirpath)
    info('Splitting BAM files...')
    for chrom, variants in variants_by_chrom.iteritems():
        chr_lengths = get_chr_lengths_from_seq(cnf.genome.seq)
        chr_lengths_dict = dict((c, l) for (c, l) in chr_lengths)
        chr_length = chr_lengths_dict[chrom]
        transcripts = get_transcipts_with_exons_from_features(cnf.features, cur_chrom=chrom)
        bams_created_before = []
        bams_by_sample = defaultdict(list)
        for variant in variants:
            variant_bams_by_sample = extract_variant_from_bams(cnf, temp_output_dirpath, transcripts, chr_length,
                                                               samples, chrom, variant, bams_created_before)
            bams_created_before.extend(variant_bams_by_sample.values())
            for sample_name, bam_fpath in variant_bams_by_sample.iteritems():
                bams_by_sample[sample_name].append(bam_fpath)
        chrom = chrom.replace('chr', '')
        for sample_name, bam_fpaths in bams_by_sample.iteritems():
            bam_fname = '{chrom}-{sample_name}.bam'.format(**locals())
            temp_combined_bam_fpath = join(temp_output_dirpath, bam_fname)
            combined_bam_fpath = join(cnf.output_dir, bam_fname)
            generate_combined_bam(cnf, bam_fpaths, temp_combined_bam_fpath, combined_bam_fpath)


def get_minimal_representation(pos, ref, alt):
    """
    Get the minimal representation of a variant, based on the ref + alt alleles in a VCF
    This is used to make sure that multiallelic variants in different datasets,
    with different combinations of alternate alleles, can always be matched directly.
    """
    pos = int(pos)
    # If it's a simple SNV, don't remap anything
    if len(ref) == 1 and len(alt) == 1:
        return pos, ref, alt
    else:
        # strip off identical suffixes
        while alt[-1] == ref[-1] and min(len(alt), len(ref)) > 1:
            alt = alt[:-1]
            ref = ref[:-1]
        # strip off identical prefixes and increment position
        while alt[0] == ref[0] and min(len(alt), len(ref)) > 1:
            alt = alt[1:]
            ref = ref[1:]
            pos += 1
        return pos, ref, alt


def extract_variant_from_bams(cnf, out_dirpath, transcripts, chr_length, samples, chrom, variant, bams_created_before):
    padding = 500
    sambamba = which('sambamba')
    transcript_names = None
    pos, ref, alt, variant_transcripts = variant['pos'], variant['ref'], variant['alt'], variant['transcripts']
    for transcript in variant_transcripts:
        transcript_exons = transcripts[(transcript, chrom)]
        for idx, exon in enumerate(transcript_exons):
            if exon['start'] <= pos <= exon['stop']:
                start, end = exon['start'], exon['stop']
                transcript_names = ','.join(variant_transcripts)
                bam_prefix = '{chrom}-{transcript_names}-{idx}-'.format(**locals())
        if transcript_names:
            break
    if not transcript_names:
        start, end = max(1, pos - padding), min(chr_length, pos + padding)
        bam_prefix = '{chrom}-{pos}-{ref}-{alt}-'.format(**locals())
    bams_by_sample = dict()
    for sample in samples:
        sample_name = sample.name.replace('-', '_')
        output_bam_fpath = join(out_dirpath, bam_prefix + '{sample_name}.bam'.format(**locals()))
        if output_bam_fpath in bams_created_before:
            continue
        if not cnf.reuse_intermediate or not verify_file(output_bam_fpath, silent=True):
            cmdline = '{sambamba} slice -h {sample.bam} {chrom}:{start}-{end} > {output_bam_fpath}'.format(**locals())
            call(cnf, cmdline)
            cmdline = '{sambamba} index {output_bam_fpath}'.format(**locals())
            call(cnf, cmdline)
        bams_by_sample[sample.name] = output_bam_fpath
    return bams_by_sample


def bam_path_to_fields(bam_path):
    return os.path.basename(bam_path).replace('.bam', '').replace('_', '-').split('-')


def bam_path_to_read_group_id(bam_path):
    return os.path.basename(bam_path.replace('chr', '').replace('.bam', ''))


def bam_path_to_dict(bam_path):
    # for example: /read_viz/22/5822/chr22-46615822-A-G.bam
    return dict(zip(['chrom', 'pos', 'ref', 'alt'], bam_path_to_fields(bam_path)[:4]))


def generate_combined_bam(cnf, bam_fpaths, temp_combined_bam_fpath, combined_bam_fpath):
    info('Combining %s bams into %s' % (len(bam_fpaths), combined_bam_fpath))

    # sorted_bam_paths = sorted(bam_fpaths, key=lambda bam_path: int(bam_path_to_dict(bam_path)['pos']))

    read_group_ids = map(bam_path_to_read_group_id, bam_fpaths)
    read_groups = [{'ID': read_group_id, 'SM': 0} for read_group_id in read_group_ids]

    out_bam = None
    for bam_fpath in bam_fpaths:
        try:
            ibam = pysam.AlignmentFile(bam_fpath, 'rb')

            if out_bam is None:
                header = {
                    'HD': {'VN': '1.4'},
                    'SQ': ibam.header['SQ'],
                    'RG': read_groups,
                }

                out_bam = pysam.AlignmentFile(temp_combined_bam_fpath, 'wb', header=header)

            # iterate over the reads
            rg_tag = (('RG', bam_path_to_read_group_id(bam_fpath)), )
            for r in ibam:
                r.tags = rg_tag
                out_bam.write(r)

            ibam.close()

        except (IOError, ValueError) as e:
            err('ERROR on file %s: %s', bam_fpath, e)
    if out_bam is not None:
        out_bam.close()

    sambamba = which('sambamba')
    cmdline = '{sambamba} sort -t {cnf.threads} {temp_combined_bam_fpath} -o {combined_bam_fpath}'.format(**locals())
    call(cnf, cmdline)
    cmdline = '{sambamba} index -o {combined_bam_fpath}'.format(**locals())
    call(cnf, cmdline)
    print combined_bam_fpath + ' saved!'


def get_transcipts_with_exons_from_features(features_file, cur_chrom=None):
    transcripts = defaultdict(list)
    with open_gzipsafe(features_file) as in_f:
        for line in in_f:
            if line.startswith('#'):
                continue
            fields = line.strip('\n').split('\t')

            chrom = fields[0]
            if cur_chrom and chrom != cur_chrom:
                continue

            feature_type = fields[6]
            if feature_type not in ['Exon', 'CDS', 'UTR']:
                continue

            start = int(fields[1])
            stop = int(fields[2])
            transcript_id = fields[8]

            exon = {
                'transcript_id': transcript_id,
                'chrom': chrom,
                'start': start,
                'stop': stop
            }
            transcripts[(transcript_id, chrom)].append(exon)
    return transcripts


def main():
    info(' '.join(sys.argv))
    info()

    parser = OptionParser(usage='Usage: ' + basename(__file__) + ' --chr chr --vcf VCF_file --samples Sample1,Sample2 '
                                                                 '--bams BAM_file1,BAM_file2 -o Output_directory '
                                                                 '--features BED_file')
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser)
    parser.add_option('-o', dest='output_dir')
    parser.add_option('--samples', dest='sample_names')
    parser.add_option('--bams', dest='bams')
    parser.add_option('--vcf', dest='vcf_fpath')
    parser.add_option('--chr', dest='chrom')
    parser.add_option('--features', dest='features', help='BED file with real CDS/Exon/Gene/Transcript regions with '
                                                          'annotations (default "features" is in system_config)')
    (opts, args) = parser.parse_args(sys.argv[1:])

    cnf = Config(opts.__dict__, determine_sys_cnf(opts), {})
    cnf.verbose = True

    if not cnf.output_dir or not cnf.vcf or not cnf.features or not cnf.chrom:
        critical(parser.usage)

    samples = [BaseSample(sample_name, None, bam=bam) for (sample_name, bam) in zip(cnf.sample_names.split(','), cnf.bams.split(','))]
    split_bams(cnf, samples, cnf.vcf_fpath)
    info('Done.')


if __name__ == '__main__':
    main()


