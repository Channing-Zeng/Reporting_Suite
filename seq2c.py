#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc

from collections import OrderedDict
from shutil import copyfile
from traceback import format_exc

import sys
import os
from os.path import join, basename, splitext, dirname
from optparse import OptionParser, SUPPRESS_HELP

import source
from source import logger
from source.calling_process import call
from source.config import Config
from source.copy_number import run_seq2c
from source.file_utils import verify_dir, adjust_system_path, verify_file, adjust_path, safe_mkdir
from source.logger import info, critical, err
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug, determine_run_cnf, \
    determine_sys_cnf, set_up_dirs, check_genome_resources
from source.targetcov.bam_and_bed_utils import verify_bam, verify_bed, count_bed_cols, prepare_beds
from source.tools_from_cnf import get_system_path

from ngs_utils.proc_args import find_bams

'''
cd /Users/vlad/Dropbox/az/analysis/Seq2C
sample2bam.tsv --bed /Users/vlad/Dropbox/reference_data/bed/NexteraRapidCapture-CRUK_SMP2-1-2-hg38.bed -g hg19-chr21 -c Patient-E7801004-normal -o .
'''


def proc_args(argv):
    info(' '.join(sys.argv))
    info()

    description = 'This script generates target QC reports for each BAM provided as an input. ' \
                  'Usage: ' + basename(__file__) + ' sample2bam.tsv --bed target.bed --contols sample1:sample2 -o results_dir'
    parser = OptionParser(description=description, usage=description)
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser)
    parser.add_option('-o', dest='output_dir', metavar='DIR', default=join(os.getcwd(), 'seq2c'))
    parser.add_option('--bed', dest='bed', help='BED file to run Seq2C analysis')
    parser.add_option('-c', '--controls', dest='controls', help='Optional control sample names for Seq2C. For multiple controls, separate them using :')
    parser.add_option('--seq2c-opts', dest='seq2c_opts', help='Options for the final lr2gene.pl script.')
    parser.add_option('--no-prep-bed', dest='prep_bed', help=SUPPRESS_HELP, action='store_false', default=True)

    (opts, args) = parser.parse_args()
    logger.is_debug = opts.debug

    if len(args) == 0:
        parser.print_usage()
        sys.exit(1)
    if len(args) == 1 and not args[0].endswith('.bam'):
        sample_names, bam_fpaths = read_samples(verify_file(args[0], is_critical=True, description='Input sample2bam.tsv'))
        bam_by_sample = OrderedDict()
        for s, b in zip(sample_names, bam_fpaths):
            bam_by_sample[s] = b
    else:
        bam_by_sample = find_bams(args)

    run_cnf = determine_run_cnf(opts, is_wgs=not opts.__dict__.get('bed'))
    cnf = Config(opts.__dict__, determine_sys_cnf(opts), run_cnf)
    check_genome_resources(cnf)

    cnf.output_dir = adjust_path(cnf.output_dir)
    verify_dir(dirname(cnf.output_dir), is_critical=True)
    safe_mkdir(cnf.output_dir)

    if not cnf.project_name:
        cnf.project_name = basename(cnf.output_dir)
    info('Project name: ' + cnf.project_name)

    cnf.proc_name = 'Seq2C'
    set_up_dirs(cnf)

    samples = [
        source.TargQC_Sample(name=s_name, dirpath=join(cnf.output_dir, s_name), bam=bam_fpath)
            for s_name, bam_fpath in bam_by_sample.items()]
    info('Samples: ')
    for s in samples:
        info('  ' + s.name)
    samples.sort(key=lambda _s: _s.key_to_sort())

    target_bed = verify_bed(cnf.bed, is_critical=True) if cnf.bed else None

    if not cnf.only_summary:
        cnf.qsub_runner = adjust_system_path(cnf.qsub_runner)
        if not cnf.qsub_runner: critical('Error: qsub-runner is not provided is sys-config.')
        verify_file(cnf.qsub_runner, is_critical=True)

    return cnf, samples, target_bed, cnf.output_dir


def read_samples(sample2bam_fpath):
    bam_fpaths = []
    sample_names = []
    bad_bam_fpaths = []

    info('Reading sample info from ' + sample2bam_fpath)
    with open(sample2bam_fpath) as f:
        for l in f:
            if l.startswith('#'):
                continue
            l = l.replace('\n', '')
            if not l:
                continue
            sample_name = None
            if len(l.split('\t')) == 2:
                sample_name, bam_fpath = l.split('\t')
            else:
                sample_name, bam_fpath = None, l
            if not verify_bam(bam_fpath):
                bad_bam_fpaths.append(bam_fpath)
            bam_fpath = verify_bam(bam_fpath, is_critical=True)
            bam_fpaths.append(bam_fpath)

            if sample_name is None:
                sample_name = basename(splitext(bam_fpath)[0])
                if sample_name.endswith('-ready'):
                    sample_name = sample_name.split('-ready')[0]
            sample_names.append(sample_name)
            info(sample_name + ': ' + bam_fpath)

    if bad_bam_fpaths:
        critical('BAM files cannot be found, empty or not BAMs:' + ', '.join(bad_bam_fpaths))

    return sample_names, bam_fpaths


def main():
    cnf, samples, bed_fpath, output_dir = proc_args(sys.argv)
    info('Processing ' + str(len(samples)) + ' samples')

    if cnf.prep_bed is not False:
        if not bed_fpath:
            info('No input BED is specified, using CDS instead from ' + str(cnf.genome.cds))
            bed_fpath = verify_bed(cnf.genome.cds, 'CDS bed file for ' + cnf.genome.name)

        seq2c_bed_fname = basename(bed_fpath)

        bed_cols = count_bed_cols(bed_fpath)
        if bed_cols < 4:
            check_genome_resources(cnf)
            _, _, _, bed_fpath = prepare_beds(cnf, None, None, bed_fpath)

        try:
            copyfile(bed_fpath, join(output_dir, seq2c_bed_fname))
        except OSError:
            err(format_exc())
            info()
        else:
            info('Seq2C bed file is saved in ' + join(output_dir, seq2c_bed_fname))

    bed_fpath = verify_bed(bed_fpath, is_critical=True, description='Input BED file')
    info('Using target ' + bed_fpath)

    run_seq2c(cnf, output_dir, samples, bed_fpath, cnf.is_wgs)


# def standalone_cnv(cnf, samples, target_bed):
#     sample2bam_tsv_fpath = join(cnf.work_dir, 'seq2c_sample2bam.tsv')
#     with open(sample2bam_tsv_fpath, 'w') as f:
#         for s in samples:
#             f.write(s.name + '\t' + s.bam + '\n')
#
#     seq2c_sh = get_system_path(cnf, join('Seq2C', 'seq2c.sh'))
#     samtools = get_system_path(cnf, 'samtools')
#     cmdl = '{seq2c_sh} {sample2bam_tsv_fpath} {target_bed} {cnf.controls} {cnf.seq2c_opts} "-q {cnf.queue}" {samtools}'.format(**locals())
#     seq2c_report_fpath = join(cnf.output_dir, source.seq2c_name + '.tsv')
#     return call(cnf, cmdl, output_fpath=seq2c_report_fpath)


if __name__ == '__main__':
    main()

