#!/usr/bin/env python
import math
import os
from collections import defaultdict, OrderedDict
from os.path import join, splitext, basename, dirname, abspath, isfile
from shutil import copyfile
from time import sleep
from traceback import format_exc

from ext_modules.simplejson import load
from joblib import Parallel, delayed

import source
from source.bcbio.bcbio_structure import BCBioStructure
from source.calling_process import call
from source.file_utils import verify_file, adjust_path, safe_mkdir, expanduser, file_transaction
from source.logger import info, err, step_greetings, critical
from source.targetcov.bam_and_bed_utils import verify_bed, number_of_mapped_reads, sambamba_depth
from source.qsub_utils import submit_job, wait_for_jobs
from source.reporting.reporting import SampleReport
from source.targetcov.Region import Region
from source.targetcov.bam_and_bed_utils import count_bed_cols, prepare_beds
from source.tools_from_cnf import get_script_cmdline, get_system_path
from source.utils import OrderedDefaultDict, get_chr_len_fpath
from source.utils import median, mean
from tools.bed_processing.find_ave_cov_for_regions import save_regions_to_seq2cov_output__nocnf


def run_seq2c_bcbio_structure(cnf, bcbio_structure):
    step_greetings('Coverage statistics for each gene for all samples')

    if cnf.prep_bed is not False:
        info('Preparing BED files')
        features_bed_fpath = cnf.features or cnf.genome.features  # only for annotation
        if cnf.bed or bcbio_structure.bed:
            _, _, _, seq2c_bed = \
                prepare_beds(cnf, features_bed=features_bed_fpath,
                    target_bed=bcbio_structure.bed, seq2c_bed=bcbio_structure.sv_bed)
        else:
            seq2c_bed = verify_bed(cnf.genome.cds)
    else:
        seq2c_bed = verify_bed(cnf.bed)

    info('Calculating normalized coverages for CNV...')
    cnv_report_fpath = run_seq2c(
        cnf, join(bcbio_structure.date_dirpath, BCBioStructure.cnv_dir),
        bcbio_structure.samples, seq2c_bed, is_wgs=cnf.is_wgs)

    # if not verify_module('matplotlib'):
    #     warn('No matplotlib, skipping plotting Seq2C')
    # else:
    #     Parallel(n_jobs=cnf.threads) \
    #         (delayed(draw_seq2c_plot)(CallCnf(cnf.__dict__), cnv_report_fpath, s.name,
    #                 cnf.output_dir, chr_lens=get_chr_lengths(cnf))
    #             for s in bcbio_structure.samples)
    #
    #     for s in bcbio_structure.samples:
    #         plot_fpath = draw_seq2c_plot(cnf, cnv_report_fpath, s.name, cnf.output_dir)
    info()
    info('*' * 70)
    if cnv_report_fpath:
        info('Seq2C:')
        if cnv_report_fpath:
            info('   ' + cnv_report_fpath)

    return [cnv_report_fpath]


# class Region:
#     def __init__(self, symbol, chrom, start, end, size, ave_depth):
#         self.symbol = symbol
#         self.chrom = chrom
#         self.start = start
#         self.end = end
#         self.size = size
#         self.ave_depth = ave_depth


# def _read_seq2c_regions_from_targetcov_report(detailed_gene_report_fpath, is_wgs=False):
#     amplicons = []
#
#     info('Parsing amplicons from from ' + detailed_gene_report_fpath)
#
#     with open(detailed_gene_report_fpath, 'r') as f:
#         for i, line in enumerate(f):
#             if (not is_wgs and 'Capture' in line) or (is_wgs and ('CDS' in line)):
#                 ts = line.split('\t')
#                 # Chr  Start  End  Size  Gene  Strand  Feature  Biotype  Min depth  Ave depth  Std dev.  W/n 20% of ave  ...
#                 chrom, s, e, size, symbol, strand, feature, biotype, min_depth, ave_depth = ts[:10]
#                 ampl = Region(
#                     gene_name=symbol, chrom=chrom, strand=strand, feature=feature,
#                     start=int(s), end=int(e), size=int(size), avg_depth=float(ave_depth))
#                 amplicons.append(ampl)
#
#     if not amplicons:
#         critical('No ' + ('"Capture"' if not is_wgs else '"CDS"') + ' record was found in ' + detailed_gene_report_fpath)
#
#     return amplicons


def seq2c_seq2cov(cnf, seq2cov, samtools, sample, bam_fpath, amplicons_bed, seq2c_output):
    sample_name = sample.name

    cmdline = '{seq2cov} -m {samtools} -z -b {bam_fpath} -N {sample_name} {amplicons_bed}'.format(**locals())
    res = call(cnf, cmdline, seq2c_output)
    if not res:
        err('Could not run seq2cov.pl for ' + sample.name)
        return None


def run_seq2c(cnf, output_dirpath, samples, seq2c_bed, is_wgs):
    step_greetings('Running Seq2C')

    bams_by_sample = dict()
    for s in samples:
        if not s.bam:
            err('No BAM file for ' + s.name)
            continue
        bams_by_sample[s.name] = s.bam
        # cnf.work_dir = join(ori_work_dir, source.targqc_name + '_' + s.name)
        # safe_mkdir(cnf.work_dir)
        # s.dedup_bam = intermediate_fname(cnf, s.bam, source.dedup_bam)
        # dedupped_bam_by_sample[s.name] = s.dedup_bam
        # if verify_bam(s.dedup_bam, silent=True):
        #     info(s.dedup_bam + ' exists')
        # else:
        #     info('Deduplicating bam file ' + s.dedup_bam)
        #     dedup_jobs.append(remove_dups(cnf, s.bam, s.dedup_bam, use_grid=True))

    # cnf.work_dir = ori_work_dir
    # wait_for_jobs(cnf, dedup_jobs)
    #
    # ok = True
    # for s in samples:
    #     if not dedupped_bam_by_sample.get(s.name) or not verify_bam(dedupped_bam_by_sample[s.name]):
    #         err('No BAM file for ' + s.name)
    #         ok = False
    # if not ok:
    #     err('No BAM files found for any sample, cannot run Seq2C.')
    #     return None

    info('Getting reads and cov stats')
    mapped_read_fpath = join(output_dirpath, 'mapped_reads_by_sample.tsv')
    mapped_read_fpath, samples = __get_mapped_reads(cnf, samples, bams_by_sample, mapped_read_fpath)
    info()
    if not mapped_read_fpath:
        return None

    combined_gene_depths_fpath = join(output_dirpath, 'cov.tsv')
    combined_gene_depths_fpath = __seq2c_coverage(cnf, samples, bams_by_sample, seq2c_bed, is_wgs, combined_gene_depths_fpath)
    info()
    if not combined_gene_depths_fpath:
        return None

    seq2c_report_fpath = join(output_dirpath, source.seq2c_name + '.tsv')
    seq2c_report_fpath = __final_seq2c_scripts(cnf, mapped_read_fpath, combined_gene_depths_fpath, seq2c_report_fpath)
    if not seq2c_report_fpath:
        return None

    info('Done. The results is ' + seq2c_report_fpath)
    return seq2c_report_fpath


def __final_seq2c_scripts(cnf, read_stats_fpath, combined_gene_depths_fpath, output_fpath):
    cov2lr = get_script_cmdline(cnf, 'perl', join('Seq2C', 'cov2lr.pl'), is_critical=True)
    cov2lr_output = join(cnf.work_dir, splitext(basename(output_fpath))[0] + '.cov2lr.tsv')

    controls = ''
    lr2gene_opt = ''
    if cnf.controls:
        controls = '-c ' + cnf.controls  # ':'.join([adjust_path(fpath) for fpath in cnf.controls.split(':')])
        lr2gene_opt = '-c'

    cmdline = '{cov2lr} -a {controls} {read_stats_fpath} {combined_gene_depths_fpath}'.format(**locals())
    call(cnf, cmdline, cov2lr_output, exit_on_error=False)
    info()

    if not verify_file(cov2lr_output):
        return None

    seq2c_opts = cnf.seq2c_opts or ''

    lr2gene = get_script_cmdline(cnf, 'perl', join('Seq2C', 'lr2gene.pl'), is_critical=True)
    cmdline = '{lr2gene} {lr2gene_opt} {seq2c_opts} {cov2lr_output}'.format(**locals())
    res = call(cnf, cmdline, output_fpath, exit_on_error=False)
    info()

    if not verify_file(output_fpath):
        return None

    return res


# def _run_cov2cnv(cnf, seq2cov, samtools, sample, bed_fpath, dedupped_bam_by_sample):
#     if not verify_file(sample.seq2cov_output_fpath):
#         safe_mkdir(dirname(sample.seq2cov_output_fpath))
#         seq2c_seq2cov(cnf, seq2cov, samtools, sample, sample.bam, bed_fpath, sample.seq2cov_output_fpath)
#
#     if not verify_file(sample.seq2cov_output_dup_fpath):
#         safe_mkdir(dirname(sample.seq2cov_output_dup_fpath))
#         seq2c_seq2cov(cnf, seq2cov, samtools, sample, dedupped_bam_by_sample[sample.name], bed_fpath, sample.seq2cov_output_dup_fpath)


def __seq2c_coverage(cnf, samples, bams_by_sample, bed_fpath, is_wgs, output_fpath):
    if cnf.reuse_intermediate and verify_file(output_fpath, silent=True):
        info(output_fpath + ' exists, reusing')
        return output_fpath

    jobs_by_sample = dict()
    depth_output_by_sample = dict()
    seq2cov_output_by_sample = dict()
    seq2c_work_dirpath = join(cnf.work_dir, source.seq2c_name)
    safe_mkdir(seq2c_work_dirpath)
    info()
    for s in samples:
        info('*' * 50)
        info(s.name + ':')
        seq2cov_output_by_sample[s.name] = join(seq2c_work_dirpath, s.name + '.seq2cov.txt')

        if not cnf.reuse_intermediate and isfile(seq2cov_output_by_sample[s.name]):
            os.remove(seq2cov_output_by_sample[s.name])

        if cnf.reuse_intermediate and verify_file(seq2cov_output_by_sample[s.name], silent=True):
            info(seq2cov_output_by_sample[s.name] + ' exists, reusing')

        elif verify_file(s.targetcov_detailed_tsv, silent=True):
            info('Using targetcov detailed output for Seq2C coverage.')
            info(s.name + ': using targetseq output')
            targetcov_details_to_seq2cov(cnf, s.targetcov_detailed_tsv, seq2cov_output_by_sample[s.name], s.name, is_wgs=is_wgs)

        else:
            info(s.name + ': ' + s.targetcov_detailed_tsv + ' does not exist: submitting sambamba depth')
            bam_fpath = bams_by_sample[s.name]
            depth_output = join(seq2c_work_dirpath, s.name + '_depth' + '.txt')
            depth_output_by_sample[s.name] = depth_output
            if cnf.reuse_intermediate and verify_file(depth_output, silent=True):
                info(depth_output + ' exists, reusing')
            else:
                j = sambamba_depth(cnf, bed_fpath, bam_fpath, depth_output, use_grid=True)
                jobs_by_sample[s.name] = j
        info()
    info('*' * 50)

    wait_for_jobs(cnf, jobs_by_sample.values())
    for s_name, j in jobs_by_sample.items():
        if j.is_done and not j.is_failed:
            info(s_name + ': summarizing bedcoverage output ' + depth_output_by_sample[s_name])
            bed_col_num = count_bed_cols(bed_fpath)
            sambamba_depth_to_seq2cov(cnf, j.output_fpath, seq2cov_output_by_sample[s_name], s_name, bed_col_num)
        else:
            err('ERROR: ' + s_name + ' could not get coverage stats, log saved to ' + j.log_fpath)

            # script = get_script_cmdline(cnf, 'python', join('tools', 'bed_processing', 'find_ave_cov_for_regions.py'),
            #                             is_critical=True)
            # bedcov_hist_fpath = depth_output_by_sample[s_name]
            # cmdline = '{script} {bedcov_hist_fpath} {s_name} {bed_col_num}'.format(**locals())
            # j = submit_job(cnf, cmdline, s_name + '_bedcov_2_seq2cov', output_fpath=seq2cov_output_by_sample[s_name])
            # sum_jobs_by_sample[s_name] = j

    # sum_jobs_by_sample = dict()
    # info('* Submitting seq2cov output *')
    # for s_name, j in jobs_by_sample.items():
    #     if not verify_file(seq2cov_output_by_sample[s_name], silent=True):
    #         info(s_name + ': summarizing bedcoverage output ' + depth_output_by_sample[s_name])
    #
    #         script = get_script_cmdline(cnf, 'python', join('tools', 'bed_processing', 'find_ave_cov_for_regions.py'),
    #                                     is_critical=True)
    #         bedcov_hist_fpath = depth_output_by_sample[s_name]
    #         bed_col_num = count_bed_cols(seq2c_bed)
    #         cmdline = '{script} {bedcov_hist_fpath} {s_name} {bed_col_num}'.format(**locals())
    #         j = submit_job(cnf, cmdline, s_name + '_bedcov_2_seq2cov', output_fpath=seq2cov_output_by_sample[s_name])
    #         sum_jobs_by_sample[s_name] = j
    #
    # wait_for_jobs(cnf, sum_jobs_by_sample.values())

    info()
    info('Done')
    info('*' * 50)
    info()
    info('Combining seq2cov output')
    with open(output_fpath, 'w') as out:
        for i, s in enumerate(samples):
            verify_file(seq2cov_output_by_sample[s.name], description='seq2cov_output for ' + s.name, is_critical=True)
            with open(seq2cov_output_by_sample[s.name]) as inp:
                for l in inp:
                    out.write(l)

    verify_file(output_fpath, description='__simulate_cov2cnv_w_bedtools output_fpath', is_critical=True)
    info('Saved combined seq2cov output to ' + output_fpath)
    info()
    return output_fpath


def targetcov_details_to_seq2cov(cnf, targetcov_detials_fpath, output_fpath, sample_name, is_wgs=False):
    info('Parsing coverage from targetcov per-regions: ' + targetcov_detials_fpath + ' -> ' + output_fpath)
    if not is_wgs:
        info('Filtering by "Capture" tag')
        fn_keep_only = lambda l: 'Capture' in l
    else:
        info('Filtering by "CDS" tag')
        fn_keep_only = lambda l: 'CDS' in l

    convert_to_seq2cov(cnf, targetcov_detials_fpath, output_fpath, sample_name,
                       chrom_col=0, start_col=1, end_col=2, gene_col=4, ave_depth_col=10,
                       keep_only=fn_keep_only)

''' chr20_normal.targetSeq.details.gene.tsv:
#Chr    Start   End     Size   Gene    Strand  Feature    Biotype         Transcript   Min depth  Ave depth      Std dev W/n 20% of ave depth    1x      5x      10x     25x     50x     100x    500x    1000x   5000x   10000x  50000x
chr20   68345   68413   68     DEFB125 .       Capture    .               .            28         32.5           1.66716 .       1.0     1.0     1.0     1.0     0       0       0       0       0       0       0
chr20   76640   77301   661    DEFB125 .       Capture    .               .            24         36.9213        5.74231 .       1.0     1.0     1.0     0.995461        0.00302572      0       0       0       0       0       0
chr20   68350   68408   58     DEFB125 +       CDS        protein_coding  NM_153325    28         32.5           1.75431 .       1.0     1.0     1.0     1.0     0       0       0       0       0       0       0
chr20   76645   77214   569    DEFB125 +       CDS        protein_coding  NM_153325    24         36.819         6.13761 .       1.0     1.0     1.0     0.994728        0.00351494      0       0       0       0       0       0
chr20   68312   77214   627    DEFB125 .       Gene-Exon  protein_coding  NM_153325    24         36.4194752791   6.00301802873   .       1.0     1.0     1.0     0.995215681021  0.00318979403509        0       0       0       0       0       0
'''

def sambamba_depth_to_seq2cov(cnf, sambamba_depth_output_fpath, output_fpath, sample_name, bed_col_num):
    info('Converting sambamba depth output to seq2cov output: ' + sambamba_depth_output_fpath + ' -> ' + output_fpath)
    assert bed_col_num >= 4, bed_col_num
    convert_to_seq2cov(cnf, sambamba_depth_output_fpath, output_fpath, sample_name,
                       chrom_col=0, start_col=1, end_col=2, gene_col=3, ave_depth_col=bed_col_num + 2)

''' sambamba_depth_output_fpath:
# chrom chromStart  chromEnd  F3       readCount  minDepth  meanCoverage  stdDev   percentage1  percentage5  percentage10  ...  sampleName
chr20   68345       68413     DEFB125  56         28        32.5          1.66716  100          100          100           ...  chr20_tumor
chr20   76640       77301     DEFB125  279        24        36.9213       5.74231  100          100          100           ...  chr20_tumor
'''
''' seq2cov:
chr20_tumor_1   DEFB125   chr20   68346   68413   Amplicon    68   28.0
chr20_tumor_1   DEFB125   chr20   76641   77301   Amplicon    661  24.0
chr20_tumor_1   DEFB125   chr20   68346   77301   Whole-Gene  729  24.3731138546
chr20_tumor_1   DEFB126   chr20   123247  123332  Amplicon    86   40.0
'''
def convert_to_seq2cov(cnf, input_fpath, output_fpath, sample_name,
                       chrom_col, start_col, end_col, gene_col, ave_depth_col,
                       keep_only=None):
    if cnf.reuse_intermediate and isfile(output_fpath) and verify_file(output_fpath):
        info(output_fpath  + ' exists, reusing')
        return output_fpath

    info('First round: collecting gene ends')
    gene_end_by_gene = defaultdict(lambda: -1)
    with open(input_fpath) as f:
        for l in f:
            if l.startswith('#'): continue
            if keep_only and not keep_only(l): continue
            fs = l.replace('\n', '').split('\t')
            if any(fs[i] == '.' for i in [chrom_col, start_col, end_col, gene_col, ave_depth_col]): continue
            end = int(fs[end_col])
            gene_name = fs[gene_col]
            gene_end_by_gene[gene_name] = max(gene_end_by_gene[gene_name], end)

    info('Second round: calculating coverage')
    total_cov_by_gene = dict()
    gene_start_by_gene = dict()
    total_size_by_gene = dict()
    with file_transaction(None, output_fpath) as tx:
        with open(input_fpath) as f, open(tx, 'w') as out:
            for l in f:
                if l.startswith('#'): continue
                if keep_only and not keep_only(l): continue
                fs = l.replace('\n', '').split('\t')
                if any(fs[i] == '.' for i in [chrom_col, start_col, end_col, gene_col, ave_depth_col]): continue

                chrom = fs[chrom_col]
                start = int(fs[start_col])
                end = int(fs[end_col])
                gene_name = fs[gene_col]
                ave_depth = float(fs[ave_depth_col])

                if gene_name not in gene_start_by_gene:
                    gene_start_by_gene[gene_name] = start
                    total_cov_by_gene[gene_name] = 0
                    total_size_by_gene[gene_name] = 0
                else:
                    gene_start_by_gene[gene_name] = min(start, gene_start_by_gene[gene_name])
                total_cov_by_gene[gene_name] += ave_depth * (end - start)
                total_size_by_gene[gene_name] += end - start

                fs = [sample_name, gene_name, chrom, str(start + 1), str(end), 'Amplicon', str(end - start), str(ave_depth)]
                out.write('\t'.join(fs) + '\n')

                if end >= gene_end_by_gene[gene_name]:
                    assert end == gene_end_by_gene[gene_name], (end, gene_end_by_gene[gene_name])
                    start = gene_start_by_gene[gene_name]
                    ave_depth = total_cov_by_gene[gene_name] / total_size_by_gene[gene_name]
                    size = total_size_by_gene[gene_name]
                    fs = [sample_name, gene_name, chrom, str(start + 1), str(end), 'Whole-Gene', str(size), str(ave_depth)]
                    out.write('\t'.join(fs) + '\n')
    info('Done, saved to ' + output_fpath)
    return output_fpath


def save_regions_to_seq2cov_output(cnf, sample_name, regions, output_fpath):
    info('*' * 50)
    info(sample_name + ':')
    if cnf.reuse_intermediate and verify_file(output_fpath, silent=True):
        info('reusing ' + output_fpath)
        return output_fpath
    else:
        info(sample_name + ': summing up whole-genes')
        with file_transaction(cnf.work_dir, output_fpath) as tx_fpath:
            save_regions_to_seq2cov_output__nocnf(sample_name, regions, tx_fpath)

        info(sample_name + ': saved seq2cov to ' + output_fpath)
        return output_fpath


def __cov2cnv(cnf, target_bed, samples, dedupped_bam_by_sample, combined_gene_depths_fpath):
    info()
    # info('Combining gene depths...')
    seq2c_work_dir = join(cnf.work_dir, source.seq2c_name)
    safe_mkdir(seq2c_work_dir)
    seq2cov_fpath_by_sample = {s.name: join(seq2c_work_dir, s.name) for s in samples}

    result = []
    features_bed_fpath = adjust_path(cnf.features) if cnf.features else adjust_path(cnf.genome.features)
    # print any(not verify_file(seq2cov_fpath_by_sample[s.name], silent=True) for s in samples)
    seq2c_bed = None
    if any(not verify_file(seq2cov_fpath_by_sample[s.name], description='seq2cov_fpath for ' + s.name,
            silent=True) for s in samples) or not cnf.reuse_intermediate:
        if cnf._prep_bed is not False:
            _, _, _, seq2c_bed = \
                prepare_beds(cnf, features_bed=features_bed_fpath, target_bed=target_bed, seq2c_bed=target_bed)
        else:
            seq2c_bed = target_bed or verify_bed(cnf.genome.cds)

    # info('Running first for the de-dupped version, then for the original version.')
    # Parallel(n_jobs=cnf.threads) \
    #     (delayed(_run_cov2cnv)(CallCnf(cnf.__dict__), seq2cov, samtools, s, bed_fpath, dedupped_bam_by_sample)
    #         for s in samples)

    qsub = get_system_path(cnf, 'qsub')
    runner_script = abspath(expanduser(cnf.qsub_runner))
    queue = cnf.queue
    bash = get_system_path(cnf, 'bash')

    seq2cov_wrap = get_system_path(cnf, 'bash', join('Seq2C', 'seq2cov_wrap.sh'), is_critical=True)
    seq2cov      = get_system_path(cnf, join('Seq2C', 'seq2cov.pl'), is_critical=True)
    wait_vardict = get_system_path(cnf, 'perl', join('Seq2C', 'waitVardict.pl'), is_critical=True)
    samtools = get_system_path(cnf, 'samtools')

    to_redo_samples = []
    for s in samples:
        if verify_file(seq2cov_fpath_by_sample[s.name], silent=True) and cnf.reuse_intermediate:
            info(seq2cov_fpath_by_sample[s.name] + ' already exist, reusing')
        else:
            if isfile(seq2cov_fpath_by_sample[s.name]):
                os.remove(seq2cov_fpath_by_sample[s.name])
            to_redo_samples.append(s)

    tx_output_fpath_by_sn = dict()
    for i, s in enumerate(to_redo_samples):
        safe_mkdir(dirname(seq2cov_fpath_by_sample[s.name]))
        bam_fpath = dedupped_bam_by_sample[s.name]
        seq2cov_output_log = join(cnf.log_dir, s.name + '.seq2cov.err')
        tx_output_fpath = join(cnf.work_dir, s.name + '.seq2cov.tx')
        tx_output_fpath_by_sn[s.name] = tx_output_fpath
        done_marker = join(cnf.work_dir, 'seq2c.done.' + s.name)
        # with file_transaction(cnf.work_dir, seq2cov_fpath_by_sample[s.name]) as tx_fpath:
        cmdline = (
            '{seq2cov_wrap} {bam_fpath} {s.name} {seq2c_bed} {s.name} {seq2cov} '
            '{samtools} {tx_output_fpath} {done_marker}').format(**locals())
        # print str(cnf.project_name)
        qsub_cmdline = (
            '{qsub} -pe smp 1 -S {bash} -q {queue} '
            '-j n -o {seq2cov_output_log} -e {seq2cov_output_log} -hold_jid \'_\' '
            '-N SEQ2C_seq2cov_{cnf.project_name}_{s.name} {runner_script} {done_marker} '
            '"{cmdline}"').format(**locals())

        info('Sumbitting seq2cov.pl for ' + s.name)
        info(qsub_cmdline)
        if isfile(done_marker):
            os.remove(done_marker)
        call(cnf, qsub_cmdline, silent=True)
        info()

    # cnt = str(len(to_redo_samples))
    # cmdline = '{wait_vardict} seq2c {cnt}'.format(**locals())
    # info('Sumbitting waitVardict.pl')
    # info(cmdline)
    # call(cnf, cmdline, silent=True)
    for s in to_redo_samples:
        tx_output_fpath = tx_output_fpath_by_sn[s.name]
        info('Waiting for ' + tx_output_fpath + ' to be written...')
        done_marker = join(cnf.work_dir, 'seq2c.done.' + s.name)
        while not isfile(done_marker):
            sleep(30)
        verify_file(tx_output_fpath, description='tx_output_fpath for ' + s.name,  is_critical=True)
        # info('tx_output_fpath = ' + tx_output_fpath)
        # info('seq2cov_fpath_by_sample[s.name] = ' + seq2cov_fpath_by_sample[s.name])
        os.rename(tx_output_fpath, seq2cov_fpath_by_sample[s.name])
        verify_file(seq2cov_fpath_by_sample[s.name], description='seq2cov_fpath for ' + s.name, is_critical=True)
        info('Done and saved to ' + seq2cov_fpath_by_sample[s.name])
        os.remove(done_marker)

    # for s in samples:
    #     if not verify_file(s.seq2cov_output_dup_fpath):
    #         safe_mkdir(dirname(s.seq2cov_output_dup_fpath))
    #         qsub_cmdline = (
    #             '{qsub} -pe smp 1 -q {queue} -V -N {cnf.project_name}_seq2cov_dup_{s.name} '
    #             '-S /bin/bash {seq2cov_wrap} {s.bam} {s.bam} {bed_fpath} {s.name} {seq2cov} {samtools} {s.seq2cov_output_dup_fpath}'
    #         ).format(**locals())
    #
    #         info('Sumbitting seq2cov.pl on duplicated bam file for ' + s.name)
    #         info(qsub_cmdline)
    #         call(cnf, qsub_cmdline, silent=True)
    #         info()
    #
    # cmdline = '{wait_vardict} seq2c $COUNTER'
    # info('Sumbitting seq2cov.pl for ' + s.name)
    # info(cmdline)
    # call(cnf, cmdline, silent=True)

    # while read i; do
    #   a=(${i//\\t/})
    #   qsub $SGE_OPT -pe smp 1 -cwd -V -N seq2c_sam${COUNTER} -S /bin/bash ${DIR}/seq2cov_wrap.sh ${a[1]} ${a[0]} $BED $COUNTER ${DIR}/seq2cov.pl
    # done < $SAM2BAM
    # perl ${DIR}/waitVardict.pl seq2c $COUNTER
    # cat cov.txt.* > cov.txt

    # for suf in '', '_dup':
    with open(combined_gene_depths_fpath, 'w') as out:
        for i, sample in enumerate(samples):
            seq2cov_fpath = seq2cov_fpath_by_sample[sample.name]

            if not verify_file(seq2cov_fpath):
                return None
            with open(seq2cov_fpath) as inp:
                for l in inp:
                    if l.startswith('Sample\tGene\tChr\t'):
                        if i == 0:
                            out.write(l)
                    else:
                        out.write(l)

    # result.append(combined_gene_depths_fpath)
    # if combined_gene_depths_fpath:
    info('Saved to ' + combined_gene_depths_fpath)
    info()

    return combined_gene_depths_fpath

    # # READ STATS
    # info('Getting read counts...')
    # bam2reads = get_script_cmdline(cnf, 'perl', join('Seq2C', 'bam2reads.pl'))
    # if not bam2reads: sys.exit(1)
    #
    # bam2reads_list_of_bams_fpath = join(cnf.work_dir, 'seq2c_list_of_bams.txt')
    # with open(bam2reads_list_of_bams_fpath, 'w') as f:
    #     for sample in samples:
    #         dedup_bam_fpath = remove_dups(cnf, sample.bam)
    #         f.write(sample.name + '\t' + dedup_bam_fpath + '\n')
    #
    # samtools = get_system_path(cnf, 'samtools')
    # read_stats_fpath = join(cnf.work_dir, 'seq2c_read_stats.txt')
    # cmdline = '{bam2reads} -m {samtools} {bam2reads_list_of_bams_fpath}'.format(**locals())
    # if not call(cnf, cmdline, read_stats_fpath): return None, None
    #
    # return read_stats_fpath, combined_gene_depths_fpath


def __get_mapped_reads(cnf, samples, bam_by_sample, output_fpath):
    if cnf.reuse_intermediate and verify_file(output_fpath, silent=True):
        info(output_fpath + ' exists, reusing')
        return output_fpath, samples

    mapped_reads_by_sample = OrderedDict()

    job_by_sample = dict()

    for s in samples:
        if verify_file(s.targetcov_json_fpath, silent=True):
            info('Parsing targetSeq output ' + s.targetcov_json_fpath)
            with open(s.targetcov_json_fpath) as f:
                data = load(f, object_pairs_hook=OrderedDict)
            cov_report = SampleReport.load(data, s)
            mapped_reads = next(rec.value for rec in cov_report.records if rec.metric.name == 'Mapped reads')
            info(s.name + ': ')
            info('  Mapped reads: ' + str(mapped_reads))
            mapped_reads_by_sample[s.name] = mapped_reads

        else:
            if s.name not in bam_by_sample:
                err('No BAM for ' + s.name + ', not running Seq2C')
                return None, None

            info('Submitting a sambamba job to get mapped read numbers')
            bam_fpath = bam_by_sample[s.name]
            j = number_of_mapped_reads(cnf, bam_fpath, dedup=True, use_grid=True)
            job_by_sample[s.name] = j

    wait_for_jobs(cnf, job_by_sample.values())
    for s_name, j in job_by_sample.items():
        if j and j.is_done and not j.is_failed:
            with open(j.output_fpath) as f:
                mapped_reads = int(f.read().strip())
                info(s_name + ': ')
                info('  Mapped reads: ' + str(mapped_reads))
                mapped_reads_by_sample[s_name] = mapped_reads
        else:
            err('ERROR: ' + s_name + ' could not get mapped reads, log saved to ' + j.log_fpath)

    with open(output_fpath, 'w') as f:
        for sample_name, mapped_reads in mapped_reads_by_sample.items():
            f.write(sample_name + '\t' + str(mapped_reads) + '\n')

    verify_file(output_fpath, is_critical=True)
    successful_samples = [s for s in samples if s.name in mapped_reads_by_sample]
    info('Samples processed: ' + str(len(samples)) + ', successfully: ' + str(len(successful_samples)))
    return output_fpath, successful_samples
