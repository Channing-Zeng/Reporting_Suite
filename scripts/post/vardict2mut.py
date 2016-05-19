#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc

from optparse import OptionParser
from os.path import join, abspath, dirname
from os.path import exists
import time
import sys

from source import verify_file
from source.config import Config
from source import logger
from source.file_utils import adjust_path, verify_dir
from source.logger import info, critical, err, warn, debug
from source.prepare_args_and_cnf import determine_run_cnf, check_genome_resources, \
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug
from source.prepare_args_and_cnf import determine_sys_cnf
from source.variants.vardict2mut_src import Filtration


def get_args():
    info(' '.join(sys.argv))
    info()
    description = (
        'The program will filter the VarDict output after vcf2txt.pl to '
        'candidate interpretable mutations, somatic or germline.')
    parser = OptionParser(description=description)
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser, threads=1)

    parser.add_option('-o', dest='output_file')
    parser.add_option('--o-all-transcripts', dest='all_transcripts_output_file')
    parser.add_option('--o-fm', dest='fm_output_file')

    parser.add_option('--cohort-freqs', dest='cohort_freqs_fpath')

    parser.add_option('-D', '--min-depth', dest='filt_depth', type='int', help='The minimum total depth')
    parser.add_option('-V', '--min-vd', dest='min_vd', type='int', help='The minimum reads supporting variant')
    parser.add_option('--gmaf', dest='min_gmaf', type='float',
                      help='When the GMAF is greater than specified, it\'s considered common SNP and filtered out.')
    parser.add_option('-f', '--min-freq', dest='min_freq', type='float',
                      help='The minimum allele frequency for regular variants. Default: 0.05')
    parser.add_option('-F', '--min-freq-hs', dest='min_hotspot_freq', type='float',
                      help='The minimum allele frequency hotspot somatic mutations, typically lower then -f. '
                           'Default: 0.01 or half -f, whichever is less')
    parser.add_option('-N', '--keep-utr-intronic', dest='keep_utr_intronic', action='store_true',
                      help='Keep all intronic and UTR in the output, but will be set as "unknown".')

    parser.add_option('-p', '--platform', dest='platform',
                      help='The platform, such as WXS, WGS, RNA-Seq, VALIDATION, etc. No Default. '
                           'Used for output in FM\'s format')

    parser.set_usage('Usage: ' + __file__ + ' vcf2txt_res_fpath [opts] -o output_fpath')

    (opts, args) = parser.parse_args()
    if len(args) < 1:
        critical('Provide the first argument - output from vcf2txt.pl')
    logger.is_debug = opts.debug

    vcf2txt_res_fpath = verify_file(args[0], is_critical=True)

    run_cnf = determine_run_cnf(opts)
    cnf = Config(opts.__dict__, determine_sys_cnf(opts), run_cnf)
    if not cnf.genome:
        critical('Please, specify the --genome option (e.g. --genome hg19)')

    check_genome_resources(cnf)

    if not cnf.output_file:
        critical('Please, specify the output fpath with -o')

    info()

    return cnf, vcf2txt_res_fpath


def main():
    cnf, vcf2txt_res_fpath = get_args()

    info('-' * 70)
    info('Writing to ' + cnf.output_file)
    if cnf.all_transcripts_output_file:
        info('Writing info for all transcripts to ' + cnf.all_transcripts_output_file)
    if cnf.fm_output_file:
        info('Writing in FM format to ' + cnf.fm_output_file)

    f = Filtration(cnf)

    input_f = open(verify_file(vcf2txt_res_fpath))
    output_f = open(adjust_path(cnf.output_file), 'w')
    fm_output_f = open(adjust_path(cnf.fm_output_file), 'w') if cnf.fm_output_file else None
    all_transcripts_output_f = open(adjust_path(cnf.all_transcripts_output_file), 'w') if cnf.all_transcripts_output_file else None

    info()
    info('-' * 70)
    info('Running filtering...')
    f.do_filtering(input_f, output_f, fm_output_f, all_transcripts_output_f)

    input_f.close()
    output_f.close()
    if fm_output_f:
        fm_output_f.close()
    if all_transcripts_output_f:
        all_transcripts_output_f.close()

    info()
    info('Saved to ' + cnf.output_file)


if __name__ == '__main__':
    main()
