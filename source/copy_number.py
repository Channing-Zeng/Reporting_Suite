#!/usr/bin/env python
import math
from collections import defaultdict, OrderedDict
import os
from os.path import join, splitext, basename, dirname, abspath, isfile
from time import sleep
from joblib import Parallel, delayed
import sys
from ext_modules.simplejson import load
from subprocess import check_output
from source.bcbio_structure import BCBioStructure
from source.calling_process import call_subprocess, call_pipe, call
from source.config import CallCnf
from source.file_utils import verify_file, adjust_path, iterate_file, safe_mkdir, expanduser, file_transaction, \
    add_suffix, splitext_plus
from source.logger import info, err, step_greetings, critical, send_email, warn
from source.ngscat.bed_file import verify_bed, verify_bam
from source.qsub_utils import submit_job, wait_for_jobs
from source.reporting import write_tsv_rows, Record, SampleReport
from source.targetcov.Region import Region, GeneInfo
from source.targetcov.bam_and_bed_utils import sort_bed, count_bed_cols, annotate_amplicons, cut, \
    group_and_merge_regions_by_gene, bedtools_version
from source.targetcov.cov import make_and_save_general_report, make_targetseq_reports, remove_dups, \
    summarize_bedcoverage_hist_stats
from source.tools_from_cnf import get_script_cmdline, get_system_path
from source.utils import OrderedDefaultDict, get_chr_len_fpath
from source.utils import median, mean
import source


def cnv_reports(cnf, bcbio_structure):
    step_greetings('Coverage statistics for each gene for all samples')

    info('Calculating normalized coverages for CNV...')
    cnv_report_fpath = _seq2c(cnf, bcbio_structure)

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


def _read_amplicons_from_targetcov_report(detailed_gene_report_fpath):
    amplicons = []

    info('Reading from ' + detailed_gene_report_fpath)

    with open(detailed_gene_report_fpath, 'r') as f:
        for i, line in enumerate(f):
            if 'Capture' in line:
                ts = line.split('\t')
                # Chr  Start  End  Size  Gene  Strand  Feature  Biotype  Min depth  Ave depth  Std dev.  W/n 20% of ave  ...
                chrom, s, e, size, symbol, strand, feature, _, ave_depth = ts[:9]
                ampl = Region(
                    gene_name=symbol, chrom=chrom, strand=strand, feature=feature,
                    start=int(s), end=int(e), size=int(size), avg_depth=float(ave_depth))
                amplicons.append(ampl)

    if not amplicons:
        critical('No Capture is found in ' + detailed_gene_report_fpath)

    return amplicons


def seq2c_seq2cov(cnf, seq2cov, samtools, sample, bam_fpath, amplicons_bed, seq2c_output):
    sample_name = sample.name

    cmdline = '{seq2cov} -m {samtools} -z -b {bam_fpath} -N {sample_name} {amplicons_bed}'.format(**locals())
    res = call(cnf, cmdline, seq2c_output)
    if not res:
        err('Could not run seq2cov.pl for ' + sample.name)
        return None


def _seq2c(cnf, bcbio_structure):
    """
    Normalize the coverage from targeted sequencing to CNV log2 ratio. The algorithm assumes the medium
    is diploid, thus not suitable for homogeneous samples (e.g. parent-child).
    """
    info('Deduplicating BAMs if needed')
    dedup_bam_dirpath = join(cnf.work_dir, source.dedup_bam)
    safe_mkdir(dedup_bam_dirpath)
    dedupped_bam_by_sample = dict()
    dedup_jobs = []
    for s in bcbio_structure.samples:
        dedup_bam_fpath = join(dedup_bam_dirpath, add_suffix(basename(s.bam), source.dedup_bam))
        dedupped_bam_by_sample[s.name] = dedup_bam_fpath
        if verify_bam(dedup_bam_fpath, silent=True):
            info(dedup_bam_fpath + ' exists')
        else:
            info('Deduplicating bam file ' + dedup_bam_fpath)
            dedup_jobs.append(remove_dups(cnf, s.bam, dedup_bam_fpath))
    dedup_jobs = wait_for_jobs(dedup_jobs)

    info('Getting reads and cov stats')
    mapped_read_fpath = join(cnf.work_dir, 'mapped_reads_by_sample.txt')
    __get_mapped_reads(cnf, bcbio_structure, dedupped_bam_by_sample, mapped_read_fpath)
    info()
    #  __cov2cnv(cnf, bcbio_structure.sv_bed or bcbio_structure.bed, bcbio_structure.samples, dedupped_bam_by_sample)
    combined_gene_depths_fpath = join(cnf.work_dir, 'gene_depths.seq2c.txt')
    __simulate_cov2cnv_w_bedtools(cnf, bcbio_structure, bcbio_structure.samples, dedupped_bam_by_sample, combined_gene_depths_fpath)
    info()

    seq2c_report_fpath = join(cnf.output_dir, BCBioStructure.seq2c_name + '.tsv')
    __new_seq2c(cnf, mapped_read_fpath, combined_gene_depths_fpath, seq2c_report_fpath)

    # info('Getting old way reads and cov stats, but with amplicons')
    # info()
    # cnv_gene_ampl_report_fpath = __cov2cnv2(cnf, read_stats_fpath, combined_gene_depths_fpath, cnv_gene_ampl_report_fpath)
    # cnv_gene_ampl_report_fpath__mine = __new_seq2c(cnf, read_stats_fpath, combined_gene_depths_fpath__mine, cnv_gene_ampl_report_fpath__mine)

    # run_copy_number__cov2cnv2(cnf, mapped_reads_by_sample, amplicon_summary_lines, cnv_ampl_report_fpath)

    # save_results_separate_for_samples(results)

    return seq2c_report_fpath


def __new_seq2c(cnf, read_stats_fpath, combined_gene_depths_fpath, output_fpath):
    cov2lr = get_script_cmdline(cnf, 'perl', join('Seq2C', 'cov2lr.pl'), is_critical=True)
    cov2lr_output = join(cnf.work_dir, splitext(basename(output_fpath))[0] + '.cov2lr.tsv')

    controls = ''
    lr2gene_opt = ''
    if cnf.controls:
        controls = '-c ' + cnf.controls  # ':'.join([adjust_path(fpath) for fpath in cnf.controls.split(':')])
        lr2gene_opt = '-c'
    #TODO: Sakina, what is usually passed to Seq2C as controls? Samples that are part of the project, or something outside of the project?

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


def __prep_bed(cnf, bed_fpath, exons_bed):
    info()
    info('Preparing BED file for seq2c...')

    info()
    info('Sorting regions by (chrom, gene name, start)')
    bed_fpath = sort_bed(cnf, bed_fpath)

    cols = count_bed_cols(bed_fpath)
    if cols < 4:
        info('Sorting exons by (chrom, gene name, start)')
        exons_bed = sort_bed(cnf, exons_bed)
        info(str(cols) + ' columns (Less than 4). Annotating amplicons with gene names from Ensembl...')
        bed_fpath = annotate_amplicons(cnf, bed_fpath, exons_bed)

    elif 8 > cols > 4:
        bed_fpath = cut(cnf, bed_fpath, 4)

    elif cols > 8:
        bed_fpath = cut(cnf, bed_fpath, 4)  # TODO: address 8-column bed-files for Seq2C

    # removing regions with no gene annotation
    def f(l, i):
        if l.split('\t')[3].strip() == '.': return None
        else: return l
    bed_fpath = iterate_file(cnf, bed_fpath, f, 'filt')

    info('Done: ' + bed_fpath)
    return bed_fpath


# def _run_cov2cnv(cnf, seq2cov, samtools, sample, bed_fpath, dedupped_bam_by_sample):
#     if not verify_file(sample.seq2cov_output_fpath):
#         safe_mkdir(dirname(sample.seq2cov_output_fpath))
#         seq2c_seq2cov(cnf, seq2cov, samtools, sample, sample.bam, bed_fpath, sample.seq2cov_output_fpath)
#
#     if not verify_file(sample.seq2cov_output_dup_fpath):
#         safe_mkdir(dirname(sample.seq2cov_output_dup_fpath))
#         seq2c_seq2cov(cnf, seq2cov, samtools, sample, dedupped_bam_by_sample[sample.name], bed_fpath, sample.seq2cov_output_dup_fpath)


def __simulate_cov2cnv_w_bedtools(cnf, bcbio_structure, samples, dedupped_bam_by_sample, output_fpath):
    if cnf.reuse_intermediate and verify_file(output_fpath, silent=True):
        info(output_fpath + ' exists, reusing')
        return output_fpath

    info('Preparing BED files')
    bed_fpath = bcbio_structure.sv_bed or bcbio_structure.bed
    exons_bed_fpath = cnf.exons if cnf.exons else cnf.genome.exons
    bed_fpath = adjust_path(bed_fpath or exons_bed_fpath)
    verify_bed(bed_fpath, is_critical=True)
    bed_fpath = __prep_bed(cnf, bed_fpath, exons_bed_fpath)

    regions_by_sample = defaultdict(dict)
    jobs_to_wait = []
    bedcov_output_by_sample = dict()
    seq2cov_output_by_sample = dict()
    chr_lengths = get_chr_len_fpath(cnf)
    seq2c_work_dirpath = join(cnf.work_dir, source.seq2c_name)
    safe_mkdir(seq2c_work_dirpath)
    info()
    for s in samples:
        info('*' * 50)
        info(s.name + ':')
        seq2cov_output_by_sample[s.name] = join(seq2c_work_dirpath, s.name + '.seq2cov.txt')
        if cnf.reuse_intermediate and verify_file(seq2cov_output_by_sample[s.name], silent=True):
            info(seq2cov_output_by_sample[s.name] + ' exists, reusing')

        elif verify_file(s.targetcov_detailed_tsv, silent=True):
            info(s.name + ': parsing targetseq output')
            amplicons = _read_amplicons_from_targetcov_report(s.targetcov_detailed_tsv)
            regions_by_sample[s.name] = amplicons

        else:
            info(s.name + ': submitting bedcoverage hist')
            bam_fpath = dedupped_bam_by_sample[s.name]
            job_name = s.name + '_bedcov'
            bedcov_output = join(seq2c_work_dirpath, job_name + '.txt')
            bedcov_output_by_sample[s.name] = bedcov_output
            if cnf.reuse_intermediate and verify_file(bedcov_output, silent=True):
                info(bedcov_output + ' exists, reusing')
            else:
                bedtools = get_system_path(cnf, 'bedtools')
                v = bedtools_version(bedtools)
                if v and v >= 24:
                    cmdline = '{bedtools} coverage -sorted -g {chr_lengths} -a {bed_fpath} -b {bam_fpath} -hist'.format(**locals())
                else:
                    cmdline = '{bedtools} coverage -abam {bam_fpath} -b {bed_fpath} -hist'.format(**locals())
                j = submit_job(cnf, cmdline, job_name, sample=s, output_fpath=bedcov_output)
                jobs_to_wait.append(j)
        info()
    info('*' * 50)

    info('* Making seq2cov output *')
    jobs_to_wait = wait_for_jobs(jobs_to_wait)
    for s in samples:
        if not regions_by_sample[s.name] and not verify_file(seq2cov_output_by_sample[s.name], silent=True):
            info(s.name + ': summarizing bedcoverage output ' + bedcov_output_by_sample[s.name])
            amplicons, _, _ = summarize_bedcoverage_hist_stats(bedcov_output_by_sample[s.name], s.name, count_bed_cols(bed_fpath))
            amplicons = sorted(amplicons, key=lambda a: (a.chrom, a.gene_name, a.start))
            for r in amplicons:
                r.calc_avg_depth()
            regions_by_sample[s.name] = amplicons

    def __sum_up_gene(g):
        g.start = g.amplicons[0].start
        g.end = g.amplicons[-1].end
        g.size = sum(a.end - a.start for a in g.amplicons)
        g.avg_depth = sum(float(a.size) * a.avg_depth for a in g.amplicons) / g.get_size() if g.get_size() else 0

    for s in samples:
        info('*' * 50)
        info(s.name + ':')
        if cnf.reuse_intermediate and verify_file(seq2cov_output_by_sample[s.name], silent=True):
            info('reusing ' + seq2cov_output_by_sample[s.name])
        else:
            info(s.name + ': summing up whole-genes')
            final_regions = []
            gene = None
            for a in regions_by_sample[s.name]:
                a.sample_name = s.name
                a.feature = 'Amplicon'
                if gene and gene.gene_name != a.gene_name:
                    __sum_up_gene(gene)
                    final_regions.append(gene)
                    gene = None
                if not gene:
                    gene = GeneInfo(sample_name=s.name, gene_name=a.gene_name, chrom=a.chrom,
                                    strand=a.strand, feature='Whole-Gene')
                gene.add_amplicon(a)
                final_regions.append(a)
            if gene:
                __sum_up_gene(gene)
                final_regions.append(gene)

            coverage_info = []
            for r in final_regions:
                if r.avg_depth is not None and r.avg_depth != 0:
                    coverage_info.append([s.name, r.gene_name, r.chrom, r.start, r.end, r.feature, r.size, r.avg_depth])

            with open(seq2cov_output_by_sample[s.name], 'w') as f:
                for fs in coverage_info:
                    f.write('\t'.join(map(str, fs)) + '\n')
            info(s.name + ': saved seq2cov to ' + seq2cov_output_by_sample[s.name])
    info()
    info('*' * 50)

    with open(output_fpath, 'w') as out:
        for i, s in enumerate(samples):
            verify_file(seq2cov_output_by_sample[s.name], is_critical=True)
            with open(seq2cov_output_by_sample[s.name]) as inp:
                for l in inp:
                    out.write(l)

    verify_file(output_fpath, is_critical=True)
    info('Saved to ' + output_fpath)
    info()
    return output_fpath


def __cov2cnv(cnf, bed_fpath, samples, dedupped_bam_by_sample, combined_gene_depths_fpath):
    info()
    # info('Combining gene depths...')

    result = []
    exons_bed_fpath = adjust_path(cnf.exons) if cnf.exons else adjust_path(cnf.genome.exons)
    bed_fpath = bed_fpath or exons_bed_fpath
    # print any(not verify_file(s.seq2cov_output_fpath, silent=True) for s in samples)
    if any(not verify_file(s.seq2cov_output_fpath, silent=True) for s in samples) or not cnf.reuse_intermediate:
        verify_bed(bed_fpath, is_critical=True)
        bed_fpath, exons_bed = __prep_bed(cnf, bed_fpath, exons_bed_fpath)

    # info('Running first for the de-dupped version, then for the original version.')
    # Parallel(n_jobs=cnf.threads) \
    #     (delayed(_run_cov2cnv)(CallCnf(cnf.__dict__), seq2cov, samtools, s, bed_fpath, dedupped_bam_by_sample)
    #         for s in samples)

    qsub = get_system_path(cnf, 'qsub')
    runner_script = abspath(expanduser(cnf.qsub_runner))
    queue = cnf.queue
    bash = get_system_path(cnf, 'bash')

    # seq2cov_wrap = get_system_path(cnf, 'bash', join('Seq2C', 'seq2cov_wrap.sh'), is_critical=True)
    # seq2cov      = get_system_path(cnf, join('Seq2C', 'seq2cov.pl'), is_critical=True)
    # wait_vardict = get_system_path(cnf, 'perl', join('Seq2C', 'waitVardict.pl'), is_critical=True)
    samtools = get_system_path(cnf, 'samtools')

    to_redo_samples = []
    for s in samples:
        if verify_file(s.seq2cov_output_fpath, silent=True) and cnf.reuse_intermediate:
            info(s.seq2cov_output_fpath + ' already exist, reusing')
        else:
            if isfile(s.seq2cov_output_fpath):
                os.remove(s.seq2cov_output_fpath)
            to_redo_samples.append(s)

    tx_output_fpath_by_sn = dict()
    for i, s in enumerate(to_redo_samples):
        safe_mkdir(dirname(s.seq2cov_output_fpath))
        bam_fpath = dedupped_bam_by_sample[s.name]
        seq2cov_output_log = join(cnf.log_dir, s.name + '.seq2cov.err')
        tx_output_fpath = join(cnf.work_dir, s.name + '.seq2cov.tx')
        tx_output_fpath_by_sn[s.name] = tx_output_fpath
        done_marker = join(cnf.work_dir, 'seq2c.done.' + s.name)
        # with file_transaction(cnf.work_dir, s.seq2cov_output_fpath) as tx_fpath:
        cmdline = (
            '{seq2cov_wrap} {bam_fpath} {s.name} {bed_fpath} {s.name} {seq2cov} '
            '{samtools} {tx_output_fpath} {done_marker}').format(**locals())
        print str(cnf.project_name)
        qsub_cmdline = (
            '{qsub} -pe smp 1 -S {bash} -q {queue} '
            '-j n -o {seq2cov_output_log} -e {seq2cov_output_err} -hold_jid \'_\' '
            '-N SEQ2C_seq2cov_{cnf.project_name}_{s.name} {runner_script} {done_marker} "{cmdline}"'
        ).format(**locals())

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
        verify_file(tx_output_fpath, is_critical=True)
        os.rename(tx_output_fpath, s.seq2cov_output_fpath)
        verify_file(s.seq2cov_output_fpath, is_critical=True)
        info('Done and saved to ' + s.seq2cov_output_fpath)
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
            seq2cov_fpath = sample.seq2cov_output_fpath

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


def __get_mapped_reads(cnf, bcbio_structure, dedupped_bam_by_sample, output_fpath):
    if cnf.reuse_intermediate and verify_file(output_fpath, silent=True):
        info(output_fpath + ' exists, reusing')
        return output_fpath

    mapped_reads_by_sample = OrderedDict()

    jobs_to_wait = []
    for s in bcbio_structure.samples:
        if verify_file(s.targetcov_json_fpath, silent=True):
            info('Parsing targetSeq output ' + s.targetcov_json_fpath)
            with open(s.targetcov_json_fpath) as f:
                data = load(f, object_pairs_hook=OrderedDict)
            cov_report = SampleReport.load(data, s, bcbio_structure)
            mapped_reads = next(rec.value for rec in cov_report.records if rec.metric.name == 'Mapped reads')
            info(s.name + ': ')
            info('  Mapped reads: ' + str(mapped_reads))
            mapped_reads_by_sample[s.name] = mapped_reads

        else:
            info('targetSeq output for ' + s.name + ' was not found; submitting a flagstat job')
            samtools = get_system_path(cnf, 'samtools')
            flagstat_fpath = join(cnf.work_dir, basename(dedupped_bam_by_sample[s.name]) + '_flag_stats')
            bam_fpath = dedupped_bam_by_sample[s.name]
            cmdline = '{samtools} flagstat {bam_fpath}'.format(**locals())
            j = submit_job(cnf, cmdline, 'flagstat_' + s.name, sample=s, output_fpath=flagstat_fpath)
            jobs_to_wait.append(j)

    # if running falgstat ourselves, finally parse its output
    jobs_to_wait = wait_for_jobs(jobs_to_wait)
    for j in jobs_to_wait:
        with open(j.output_fpath) as f:
            lines = f.readlines()
            mapped_reads = int(next(l.split()[0] for l in lines if 'mapped' in l))
            info(j.sample.name + ': ')
            info('  Mapped reads: ' + str(mapped_reads))
            mapped_reads_by_sample[j.sample.name] = mapped_reads

    with open(output_fpath, 'w') as f:
        for sample_name, mapped_reads in mapped_reads_by_sample.items():
            f.write(sample_name + '\t' + str(mapped_reads) + '\n')

    verify_file(output_fpath, is_critical=True)
    return output_fpath


# def save_results_separate_for_samples(results):
#     header = results[0]
#
#     results_per_sample = OrderedDict()
#
#     for fields in results[1:]:
#         sample_name = fields[0]
#         if sample_name not in results_per_sample:
#             results_per_sample[sample_name] = [header]
#
#         results_per_sample[sample_name].append(fields)
#
#     for sample_name, fields in results_per_sample.items():


def __proc_sample(sample, norm_depths_by_sample, factors_by_gene, med_depth, median_depth_by_sample):
    info('  Processing ' + sample + '...')
    norm_depths_by_gene = OrderedDefaultDict(dict)
    norm2 = OrderedDefaultDict(dict)
    norm3 = OrderedDefaultDict(dict)

    i = 0
    for gene, norm_depth_by_sample in norm_depths_by_sample.items():
        if i % 1000 == 0:
            info('    ' + sample + ': processing gene #{0:,}: {1}'.format(i, str(gene)))
        i += 1

        if sample not in norm_depth_by_sample:
            continue

        gene_norm_depth = norm_depth_by_sample[sample] * factors_by_gene[gene] + 0.1  # norm1b

        norm_depths_by_gene[gene][sample] = gene_norm_depth

        norm2[gene][sample] = math.log(gene_norm_depth / med_depth, 2) if med_depth else 0

        norm3[gene][sample] = math.log(gene_norm_depth / median_depth_by_sample[sample], 2) if \
            median_depth_by_sample[sample] else 0

    info('  ' + sample + ': Done. Processed {0:,} genes.'.format(i))

    return norm_depths_by_gene, norm2, norm3


def my_copy_number(mapped_reads_by_sample, gene_depth):
    mapped_reads_by_sample = {k: v for k, v in mapped_reads_by_sample.items() if 'Undetermined' not in k}

    info('Parsing rows...')
    records = _report_row_to_objects(gene_depth)

    info('Calculating depths normalized by samples...')
    norm_depths_by_sample = _get_norm_depths_by_sample(mapped_reads_by_sample, records)

    med_depth = median([depth for gn, vs in norm_depths_by_sample.items() for sn, depth in vs.items()])

    info('Getting factors by genes...')
    factors_by_gene = _get_factors_by_gene(norm_depths_by_sample, med_depth)

    info('Getting median norm...')
    median_depth_by_sample = _get_med_norm_depths(mapped_reads_by_sample, norm_depths_by_sample)

    info()
    info('Norm depth by gene, norm2, norm3...')

    samples = set(rec.sample_name for rec in records)

    norm_depths_by_gene = OrderedDefaultDict(dict)
    norm2 = OrderedDefaultDict(dict)
    norm3 = OrderedDefaultDict(dict)
    for sample in samples:
        info('  Processing ' + sample + '...')
        i = 0
        for gene, norm_depth_by_sample in norm_depths_by_sample.items():
            if i % 1000 == 0:
                info('    ' + sample + ': processing gene #{0:,}: {1}'.format(i, str(gene)))
            i += 1

            if sample not in norm_depth_by_sample:
                continue

            gene_norm_depth = norm_depth_by_sample[sample] * factors_by_gene[gene] + 0.1  # norm1b
            norm_depths_by_gene[gene][sample] = gene_norm_depth

            norm2[gene][sample] = math.log(gene_norm_depth / med_depth, 2) if med_depth else 0
            norm3[gene][sample] = math.log(gene_norm_depth / median_depth_by_sample[sample], 2) if median_depth_by_sample[sample] else 0

        info('  ' + sample + ': Done. Processed {0:,} genes.'.format(i))

    return _make_report_data(records, norm2, norm3, norm_depths_by_gene, norm_depths_by_sample)


def _get_med_norm_depths(mapped_reads_by_sample, norm_depths_by_sample):
    mean_depths_by_sample = defaultdict(list)

    for sample_name in mapped_reads_by_sample:
        for gn, v_by_sample in norm_depths_by_sample.items():
            if sample_name not in v_by_sample:
                info(sample_name + ' - no ' + gn)
            mean_depths_by_sample[sample_name].append(v_by_sample.get(sample_name) or 0)

    median_depth_by_sample = OrderedDict()
    for sample_name in mapped_reads_by_sample:
        median_depth_by_sample[sample_name] = median(mean_depths_by_sample[sample_name])

    return median_depth_by_sample


# mean_reads =  mean for all the samples mapped reads
# return factor_by_sample = sample mapped read/mean for all the samples mapped reads
def _get_factors_by_sample(mapped_reads_by_sample):
    factor_by_sample = dict()

    mean_reads = mean(mapped_reads_by_sample.values())
    for sample, mapped_reads in mapped_reads_by_sample.items():
        factor_by_sample[sample] = mean_reads / mapped_reads if mapped_reads else 0

    return factor_by_sample


#med_depth = median for all gene
#return factor_by_gene = median gene for all samples/ med_depth
def _get_factors_by_gene(norm_depths_1, med_depth):
    # mean_depth_by_genes = defaultdict(list)
    # for gene_name, values in norm_depths_by_gene:
    #     mean_depth_by_genes[gene_name].append(rec.mean_depth)

    factors_by_gene = dict()
    for gene_name, values in norm_depths_1.items():
        med = median(values.values())
        factors_by_gene[gene_name] = med_depth / med if med else 0

    return factors_by_gene


#MeanDepth_Norm1
# gene -> { sample -> [] }
def _get_norm_depths_by_sample(mapped_reads_by_sample, records):
    norm_depths = defaultdict(dict)

    factor_by_sample = _get_factors_by_sample(mapped_reads_by_sample)

    # Initialization
    for sample, factor in factor_by_sample.items():
        for rec in records:
            norm_depths[rec.gene_name][sample] = 0

    # Filling the dict with available balues
    for sample, factor in factor_by_sample.items():
        info('  ' + sample)

        for rec in [rec for rec in records if rec.sample_name == sample]:
            norm_depths[rec.gene_name][sample] = rec.mean_depth * factor

    return norm_depths


def _make_report_data(records, norm2, norm3, norm_depths_by_gene, norm_depths_by_seq_distr):
    header = ["Sample", "Gene", "Chr", "Start", "Stop", "Length", "MeanDepth", "MeanDepth_Norm1",
              "MeanDepth_Norm2", "log2Ratio_norm1", "log2Ratio_norm2"]
    report_data = []

    for rec in records:
        gene_name = rec.gene_name
        sample = rec.sample_name
        report_data.append(map(str,
           (sample, gene_name, rec.chrom, rec.start_position,
            rec.end_position, rec.size,
            '{0:.3f}'.format(rec.mean_depth, norm_depths_by_seq_distr[gene_name][sample]),
            '{0:.3f}'.format(norm_depths_by_gene[gene_name][sample]),
            '{0:.3f}'.format(norm2[gene_name][sample]),
            '{0:.3f}'.format(norm3[gene_name][sample]))))

    # report_data = sorted(report_data, key=itemgetter(2, 1))
    report_data = [header] + report_data

    return report_data


# sample_order = '''PC9_IRLR PC9_9291-L0B_1 PC9_IR-GM PC9_9291_6 PC9_9291-5 PC9_IR-4
# PC9_IR-6 PC9_9291-2 PC9 PC9_9291-3 gDNA'''.split()
#
# gene_order = 'PIK3CA KRAS NRAS HRAS BRAF EGFR'.split()

def _report_row_to_objects(gene_depth):
    inputs = OrderedDefaultDict(OrderedDict)
    # for gene_name in gene_order:
    #     for sample_name in sample_order:
    #         inputs[gene_name][sample_name] = None

    for read in gene_depth:
        rec = CovRec(*read)
        if 'Undetermined' not in rec.sample_name:
            inputs[rec.gene_name][rec.sample_name] = rec

    details_list = []
    for gene_dict in inputs.values():
        for rec in gene_dict.values():
            details_list.append(rec)

    return details_list


class CovRec:
    def __init__(self, sample_name=None, chrom=None, start_position=None, end_position=None, gene_name=None,
                 type="Gene-Capture", size=None, mean_depth=None):
        self.sample_name = sample_name
        self.chrom = chrom
        self.start_position = int(start_position)
        self.end_position = int(end_position)
        self.gene_name = gene_name
        self.type = type
        self.size = int(size)
        self.mean_depth = float(mean_depth)

    def __str__(self):
        values = [self.sample_name, self.chrom, self.start_position, self.end_position, self.gene_name, self.type, self.size,
                  self.mean_depth]
        return '"' + '\t'.join(map(str, values)) + '"'

    def __repr__(self):
        return repr(
            (self.gene_name, self.sample_name))

