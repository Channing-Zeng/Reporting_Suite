#!/usr/bin/env python
import math
from collections import defaultdict, OrderedDict
import os
from os.path import join, splitext, basename, dirname, abspath, isfile
from shutil import copyfile
from time import sleep
from traceback import format_exc

from ext_modules.simplejson import load

from joblib import Parallel, delayed
from source.bcbio.bcbio_structure import BCBioStructure
from source.calling_process import call
from source.config import CallCnf
from source.file_utils import verify_file, adjust_path, safe_mkdir, expanduser, file_transaction, \
    verify_module, intermediate_fname
from source.logger import info, err, step_greetings, critical, warn
from source.targetcov.bam_and_bed_utils import verify_bam, bam_to_bed
from source.qsub_utils import submit_job, wait_for_jobs
from source.reporting.reporting import SampleReport
from source.targetcov.Region import Region
from source.targetcov.bam_and_bed_utils import count_bed_cols, bedtools_version, prepare_beds, remove_dups
from source.targetcov.coverage_hist import launch_bedcoverage_hist
from source.tools_from_cnf import get_script_cmdline, get_system_path
from source.utils import OrderedDefaultDict, get_chr_len_fpath, get_chr_lengths
from source.utils import median, mean
import source
from tools.bed_processing.find_ave_cov_for_regions import save_regions_to_seq2cov_output__nocnf


def cnv_reports(cnf, bcbio_structure):
    step_greetings('Coverage statistics for each gene for all samples')

    info('Calculating normalized coverages for CNV...')
    cnv_report_fpath = _seq2c(cnf, bcbio_structure)

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


def _read_amplicons_from_targetcov_report(detailed_gene_report_fpath, is_wgs=False):
    amplicons = []

    info('Parsing amplicons from from ' + detailed_gene_report_fpath)

    with open(detailed_gene_report_fpath, 'r') as f:
        for i, line in enumerate(f):
            if (not is_wgs and 'Capture' in line) or (is_wgs and ('\tExon' in line or 'CDS' in line)):
                ts = line.split('\t')
                # Chr  Start  End  Size  Gene  Strand  Feature  Biotype  Min depth  Ave depth  Std dev.  W/n 20% of ave  ...
                chrom, s, e, size, symbol, strand, feature, biotype, min_depth, ave_depth = ts[:10]
                ampl = Region(
                    gene_name=symbol, chrom=chrom, strand=strand, feature=feature,
                    start=int(s) + 1, end=int(e), size=int(size), avg_depth=float(ave_depth))
                amplicons.append(ampl)

    if not amplicons:
        critical('No "Capture" record was found in ' + detailed_gene_report_fpath)

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
    # dedup_bam_dirpath = join(cnf.work_dir, source.dedup_bam)
    # safe_mkdir(dedup_bam_dirpath)
    dedupped_bam_by_sample = dict()
    dedup_jobs = []
    for s in bcbio_structure.samples:
        s.dedup_bam = intermediate_fname(cnf, s.bam, source.dedup_bam)
        # s.dedup_bam = add_suffix(s.bam, source.dedup_bam)
        dedupped_bam_by_sample[s.name] = s.dedup_bam
        if verify_bam(s.dedup_bam, silent=True):
            info(s.dedup_bam + ' exists')
        else:
            info('Deduplicating bam file ' + s.dedup_bam)
            dedup_jobs.append(remove_dups(cnf, s.bam, s.dedup_bam, use_grid=True))
    dedup_jobs = wait_for_jobs(cnf, dedup_jobs)

    info('Getting reads and cov stats')
    mapped_read_fpath = join(cnf.output_dir, 'mapped_reads_by_sample.tsv')
    __get_mapped_reads(cnf, bcbio_structure, dedupped_bam_by_sample, mapped_read_fpath)
    info()

    combined_gene_depths_fpath = join(cnf.output_dir, 'cov.tsv')
    # __cov2cnv(cnf, bcbio_structure.sv_bed or bcbio_structure.bed, bcbio_structure.samples, dedupped_bam_by_sample, combined_gene_depths_fpath)
    __simulate_cov2cnv_w_bedtools(cnf, bcbio_structure, bcbio_structure.samples,
                                  dedupped_bam_by_sample, combined_gene_depths_fpath)
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
    exons_bed_fpath = cnf.exons if cnf.exons else cnf.genome.exons  # only for annotation
    _, _, target_bed, seq2c_bed = \
        prepare_beds(cnf, exons_bed=exons_bed_fpath, target_bed=bcbio_structure.bed or cnf.genome.refseq)

    output_dirpath = dirname(output_fpath)
    seq2c_exposed_fpath = join(output_dirpath, 'seq2c_target.bed')
    try:
        copyfile(seq2c_bed, seq2c_exposed_fpath)
    except OSError:
        err(format_exc())
        info()
    else:
        info('Seq2C bed file is saved in ' + seq2c_exposed_fpath)

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

        if not cnf.reuse_intermediate and isfile(seq2cov_output_by_sample[s.name]):
            os.remove(seq2cov_output_by_sample[s.name])

        if cnf.reuse_intermediate and verify_file(seq2cov_output_by_sample[s.name], silent=True):
            info(seq2cov_output_by_sample[s.name] + ' exists, reusing')

        elif target_bed == seq2c_bed and verify_file(s.targetcov_detailed_tsv, silent=True):
            info('Target and Seq2C bed are the same after correction. Using bedcoverage output for Seq2C coverage.')
            info(s.name + ': parsing targetseq output')
            amplicons = _read_amplicons_from_targetcov_report(s.targetcov_detailed_tsv,
                                                              is_wgs=(bcbio_structure.bed is None))
            amplicons = (a for a in amplicons if a.gene_name and a.gene_name != '.')
            save_regions_to_seq2cov_output(cnf, s.name, amplicons, seq2cov_output_by_sample[s.name])

        else:
            if target_bed != seq2c_bed:
                info('target_bed ' + target_bed + ' != seq2c_bed ' + seq2c_bed + ', cannot reuse ' +
                     s.targetcov_detailed_tsv + ' for Seq2C')
            if not verify_file(s.targetcov_detailed_tsv, silent=True):
                info(s.targetcov_detailed_tsv + ' does not exist, regenerating hist for Seq2C')

            info(s.name + ': submitting bedcoverage hist')
            bam_fpath = dedupped_bam_by_sample[s.name]
            # Need to convert BAM to BED to make bedtools histogram
            bedcov_output = join(seq2c_work_dirpath, s.name + '_bedcov' + '.txt')
            bedcov_output_by_sample[s.name] = bedcov_output
            if cnf.reuse_intermediate and verify_file(bedcov_output, silent=True):
                info(bedcov_output + ' exists, reusing')
            else:
                j = launch_bedcoverage_hist(cnf, seq2c_bed, bam_fpath,
                                            bedcov_output_fpath=bedcov_output, qsub=True, sample=s)
                jobs_to_wait.append(j)
        info()
    info('*' * 50)

    jobs_to_wait = wait_for_jobs(cnf, jobs_to_wait)

    sum_jobs_to_wait = []
    info('* Submitting seq2cov output *')
    for j in jobs_to_wait:
        s = j.sample
        if not verify_file(seq2cov_output_by_sample[s.name], silent=True):
            info(s.name + ': summarizing bedcoverage output ' + bedcov_output_by_sample[s.name])

            script = get_script_cmdline(cnf, 'python', join('tools', 'bed_processing', 'find_ave_cov_for_regions.py'),
                                        is_critical=True)
            bedcov_hist_fpath = bedcov_output_by_sample[s.name]
            bed_col_num = count_bed_cols(seq2c_bed)
            cmdline = '{script} {bedcov_hist_fpath} {s.name} {bed_col_num}'.format(**locals())
            j = submit_job(cnf, cmdline, s.name + '_bedcov_2_seq2cov', sample=s,
                           output_fpath=seq2cov_output_by_sample[s.name])
            sum_jobs_to_wait.append(j)

    sum_jobs_to_wait = wait_for_jobs(cnf, sum_jobs_to_wait)

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
    exons_bed_fpath = adjust_path(cnf.exons) or adjust_path(cnf.genome.exons)
    # print any(not verify_file(seq2cov_fpath_by_sample[s.name], silent=True) for s in samples)
    seq2c_bed = None
    if any(not verify_file(seq2cov_fpath_by_sample[s.name], description='seq2cov_fpath for ' + s.name,
                           silent=True) for s in samples) or not cnf.reuse_intermediate:
        _, _, _, seq2c_bed = \
            prepare_beds(cnf, exons_bed=exons_bed_fpath, target_bed=target_bed)

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
    jobs_to_wait = wait_for_jobs(cnf, jobs_to_wait)
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

