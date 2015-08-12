from itertools import dropwhile
from os.path import isfile, join, abspath, basename
import sys
from subprocess import check_output
from source.calling_process import call
from source.file_utils import intermediate_fname, iterate_file
from source.logger import info, critical, warn, err
from source.qsub_utils import submit_job
from source.tools_from_cnf import get_system_path, get_script_cmdline


def index_bam(cnf, bam_fpath, samtools=None):
    indexed_bam = bam_fpath + '.bai'
    if not isfile(bam_fpath + '.bai'):
        info('Indexing to ' + indexed_bam + '...')
        if samtools is None:
            samtools = get_system_path(cnf, 'samtools', is_critical=True)
        cmdline = '{samtools} index {bam_fpath}'.format(**locals())
        call(cnf, cmdline)
    info('Index: ' + indexed_bam)


def index_bam_grid(cnf, bam, samtools=None):
    info('Indexing to ' + bam + '...')
    samtools = samtools or get_system_path(cnf, 'samtools')
    if samtools is None:
        samtools = get_system_path(cnf, 'samtools', is_critical=True)
    cmdline = '{samtools} index {bam}'.format(**locals())  # -F (=not) 1024 (=duplicate)
    j = submit_job(cnf, cmdline, basename(bam) + '_dedup', output_fpath=bam + '.bai', stdout_to_outputfile=False)
    info()
    return j


def count_bed_cols(bed_fpath):
    with open(bed_fpath) as f:
        for l in f:
            if l and l.strip() and not l.startswith('#'):
                return len(l.split('\t'))
    # return len(next(dropwhile(lambda x: x.strip().startswith('#'), open(bed_fpath))).split('\t'))
    err('Empty bed file: ' + bed_fpath)
    return None


def remove_comments(cnf, bed_fpath):
    def f(l, i):
        if not l.startswith('#'):
            return l
        else:
            return None
    return iterate_file(cnf, bed_fpath, f, suffix='rmcmt')


def prepare_beds(cnf, exons_bed=None, target_bed=None):
    if exons_bed is None and target_bed is None:
        warn('No bed specified (WGS?) and no exons in the system config; cannot run on target.')
        return None, None

    if exons_bed and target_bed and abspath(exons_bed) == abspath(target_bed):
        warn('Same file used for exons and amplicons: ' + exons_bed)

    # Exons
    exons_no_genes_bed = None
    if exons_bed:
        info()
        info('Merging regions within genes...')
        exons_bed = group_and_merge_regions_by_gene(cnf, exons_bed, keep_genes=True)

        info()
        info('Sorting exons by (chrom, gene name, start)')
        exons_bed = sort_bed(cnf, exons_bed, cnf.genome.name)

        info()
        info('Filtering exon bed file to have only non-gene records...')
        exons_no_genes_bed = intermediate_fname(cnf, exons_bed, 'no_genes_cut')
        call(cnf, 'grep -vw Gene ' + exons_bed, output_fpath=exons_no_genes_bed)

    if target_bed:
        info()
        info('Remove comments in target...')
        target_bed = remove_comments(cnf, target_bed)

        info()
        info('Cut -f1,2,3,4 target...')
        target_bed = cut(cnf, target_bed, 4)

        info()
        info('Sorting target...')
        target_bed = sort_bed(cnf, target_bed, cnf.genome.name)

    if target_bed:
        cols = count_bed_cols(target_bed)
        if cnf.reannotate or cols < 4:
            info()
            if not exons_bed:
                critical(str(cols) + ' columns (less than 4), and no Ensembl exons to annotate regions.')
            info('cnf.reannotate is ' + str(cnf.reannotate) + ', and cols in amplicons bed is ' + str(cols) +
                 '. Annotating amplicons with gene names from Ensembl...')
            target_bed = annotate_amplicons(cnf, target_bed, exons_bed)

    seq2c_bed = prep_bed_for_seq2c(cnf, target_bed or cut(cnf, exons_no_genes_bed, 4))

    if target_bed:
        info()
        info('Merging amplicons...')
        target_bed = group_and_merge_regions_by_gene(cnf, target_bed, keep_genes=False)

        info('Sorting exons by (chrom, gene name, start)')
        target_bed = sort_bed(cnf, target_bed, cnf.genome.name)

    return exons_bed, exons_no_genes_bed, target_bed, seq2c_bed


def annotate_amplicons(cnf, amplicons_bed, exons_bed):
    output_fpath = intermediate_fname(cnf, amplicons_bed, 'ann')

    annotate_bed_py = get_system_path(cnf, 'python', join('tools', 'bed_processing', 'annotate_bed.py'))
    bedtools = get_system_path(cnf, 'bedtools')

    cmdline = '{annotate_bed_py} {amplicons_bed} {cnf.work_dir} {exons_bed} {bedtools}'.format(**locals())
    call(cnf, cmdline, output_fpath)

    return output_fpath


def group_and_merge_regions_by_gene(cnf, bed_fpath, keep_genes=False):
    output_fpath = intermediate_fname(cnf, bed_fpath, 'grp_mrg')

    group_merge_bed_py = get_system_path(cnf, 'python', join('tools', 'bed_processing', 'group_and_merge_by_gene.py'))

    cmdline = '{group_merge_bed_py} {bed_fpath}'.format(**locals())
    if not keep_genes:
        cmdline += ' | grep -vw Gene'

    call(cnf, cmdline, output_fpath)

    return output_fpath


def cut(cnf, fpath, col_num):
    cut_fpath = intermediate_fname(cnf, fpath, 'cut')
    cmdline = 'cut -f' + ','.join(map(str, range(1, col_num + 1))) + ' ' + fpath
    call(cnf, cmdline, cut_fpath)
    return cut_fpath


def prep_bed_for_seq2c(cnf, bed_fpath):
    info()
    info('Preparing BED file for seq2c...')

    cols = count_bed_cols(bed_fpath)

    seq2c_bed = None
    if 8 > cols > 4:
        seq2c_bed = cut(cnf, bed_fpath, 4)
    elif cols > 8:
        seq2c_bed = cut(cnf, bed_fpath, 8)
    else:
        seq2c_bed = bed_fpath

    # removing regions with no gene annotation
    def f(l, i):
        if l.split('\t')[3].strip() == '.': return None
        else: return l
    seq2c_bed = iterate_file(cnf, seq2c_bed, f, suffix='filt')

    info('Done: ' + seq2c_bed)
    return seq2c_bed


def filter_bed_with_gene_set(cnf, bed_fpath, gene_names_set):
    def fn(l, i):
        if l:
            fs = l.split('\t')
            new_gns = []
            for g in fs[3].split(','):
                if g in gene_names_set:
                    new_gns.append(g)
            if new_gns:
                return l.replace(fs[3], ','.join(new_gns))

    return iterate_file(cnf, bed_fpath, fn, suffix='key', check_result=False)


def sort_bed(cnf, bed_fpath, genome=None):
    output_fpath = intermediate_fname(cnf, bed_fpath, 'sorted')

    cmdl = get_script_cmdline(cnf, 'python', join('tools', 'bed_processing', 'sort_bed.py'), is_critical=True)
    if genome: cmdl += ' ' + genome

    res = call(cnf, cmdl, stdin_fpath=bed_fpath, output_fpath=output_fpath)
    if not res:
        return None

    # sort = get_system_path(cnf, 'sort')
    # cmdline = '{sort} -V -k1,1 -k2,2 -k3,3 {bed_fpath}'.format(**locals())
    # res = call(cnf, cmdline, output_fpath, exit_on_error=False)
    # if not res:
    #     warn('Cannot sort with -V, trying with -n')
    #     cmdline = '{sort} -k1,1n -k2,2n -k3,3n {bed_fpath}'.format(**locals())
    #     res = call(cnf, cmdline, output_fpath, exit_on_error=False)
    #     if not res:
    #         warn('Cannot uniq-sort, trying with bedtools')
    #         bedtools = get_system_path(cnf, 'bedtools')
    #         cmdline = '{bedtools} sort -i {bed_fpath}'.format(**locals())
    #         res = call(cnf, cmdline, output_fpath)
    #
    # if genome != 'mm10':
    #     cmdline = 'grep "^chrM" {output_fpath} > {output_fpath}_1; grep -v "^chrM" {output_fpath} >> {output_fpath}_1; mv {output_fpath}_1 {output_fpath}'.format(**locals())
    #     res = call(cnf, cmdline)
    # else:
    #     cmdline = 'grep "^chrM" {output_fpath} > {output_fpath}_1; grep -v "^chrM" {output_fpath} >> {output_fpath}_1; mv {output_fpath}_1 {output_fpath}'.format(**locals())

    return output_fpath


def total_merge_bed(cnf, bed_fpath):
    bedops = get_system_path(cnf, 'bedops')
    if bedops:
        cmdline = '{bedops} --merge {bed_fpath}'.format(**locals())
        output_fpath = intermediate_fname(cnf, bed_fpath, 'total_merged')
        call(cnf, cmdline, output_fpath)
        return output_fpath
    else:
        bedtools = get_system_path(cnf, 'bedtools')
        cmdline = '{bedtools} merge -i {bed_fpath}'.format(**locals())
        output_fpath = intermediate_fname(cnf, bed_fpath, 'total_merged')
        call(cnf, cmdline, output_fpath)
        return output_fpath


def calc_sum_of_regions(bed_fpath):
    total_bed_size = 0

    with open(bed_fpath) as f:
        for l in f:
            l = l.strip()
            if not l.startswith('#'):
                start, end = [int(f) for f in l.split('\t')[1:3]]
                total_bed_size += end - start

    return total_bed_size


def get_total_bed_size(cnf, bed_fpath):
    merged_bed = total_merge_bed(cnf, bed_fpath)
    return calc_sum_of_regions(merged_bed)


def bedtools_version(bedtools):
    v = check_output([bedtools, '--version'])  # bedtools v2.24.0
    try:
        v = int(v.split(' ')[1].split('.')[1])
    except:
        return None
    else:
        return v
