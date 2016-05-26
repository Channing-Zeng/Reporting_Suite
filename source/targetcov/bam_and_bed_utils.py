import os
import subprocess
from itertools import dropwhile
from os.path import isfile, join, abspath, basename, dirname, getctime
import sys
from subprocess import check_output
from collections import OrderedDict

from source.calling_process import call, call_pipe
from source.config import with_cnf
from source.file_utils import intermediate_fname, iterate_file, splitext_plus, verify_file, adjust_path, add_suffix, \
    safe_mkdir, file_transaction
from source.logger import info, critical, warn, err, debug
from source.qsub_utils import submit_job
from source.targetcov.Region import SortableByChrom
from source.tools_from_cnf import get_system_path, get_script_cmdline
from source.utils import md5, get_chr_lengths_from_seq, get_ext_tools_dirname


def index_bam(cnf, bam_fpath, sambamba=None, samtools=None, use_grid=False):
    if use_grid:
        return index_bam_grid(cnf, bam_fpath, sambamba)

    sambamba = sambamba or get_system_path(cnf, 'sambamba')
    indexed_bam = bam_fpath + '.bai'
    if not isfile(indexed_bam) or getctime(indexed_bam) < getctime(bam_fpath):
        info('Indexing BAM, writing ' + indexed_bam + '...')
        cmdline = '{sambamba} index {bam_fpath}'.format(**locals())
        res = call(cnf, cmdline, exit_on_error=False)
        if not isfile(indexed_bam) or getctime(indexed_bam) < getctime(bam_fpath):
            samtools = samtools or get_system_path(cnf, 'samtools')
            cmdline = '{samtools} index {bam_fpath}'.format(**locals())
            call(cnf, cmdline)
    else:
        debug('Actual "bai" index exists.')


def index_bam_grid(cnf, bam_fpath, sambamba=None):
    indexed_bam = bam_fpath + '.bai'
    if not isfile(indexed_bam) or getctime(indexed_bam) < getctime(bam_fpath):
        info('Indexing BAM, writing ' + indexed_bam + '...')
        sambamba = sambamba or get_system_path(cnf, 'sambamba')
        if sambamba is None:
            sambamba = get_system_path(cnf, 'sambamba', is_critical=True)
        cmdline = '{sambamba} index {bam_fpath}'.format(**locals())  # -F (=not) 1024 (=duplicate)
        j = submit_job(cnf, cmdline, basename(bam_fpath) + '_index', output_fpath=indexed_bam, stdout_to_outputfile=False)
        info()
        return j
    else:
        debug('Actual "bai" index exists.')
        return None


def markdup_sam(cnf, in_sam_fpath, samblaster=None):
    """Perform non-stream based deduplication of SAM input files using samblaster.
    """
    if not samblaster:
        samblaster = get_system_path(cnf, 'samblaster')
        if not samblaster:
            warn('No samblaster, can\'t mark duplicates.')
            return None

    out_sam_fpath = add_suffix(in_sam_fpath, 'markdup')
    tmp_fpath = join(cnf.work_dir, splitext_plus(basename(in_sam_fpath))[0] + '_markdup')
    safe_mkdir(dirname(tmp_fpath))
    cmdline = ('{samblaster} -i {in_sam_fpath} -o {out_sam_fpath}').format(**locals())
    res = call(cnf, cmdline, output_fpath=out_sam_fpath, stdout_to_outputfile=False, exit_on_error=False)
    if res:
        return out_sam_fpath
    else:
        return None


def markdup_bam(cnf, in_bam_fpath, bammarkduplicates=None):
    """Perform non-stream based deduplication of BAM input files using biobambam.
    """
    if not bammarkduplicates:
        bammarkduplicates = get_system_path(cnf, 'bammarkduplicates')
        if not bammarkduplicates:
            warn('No biobambam bammarkduplicates, can\'t mark duplicates.')
            return None

    out_bam_fpath = add_suffix(in_bam_fpath, 'markdup')
    tmp_fpath = join(cnf.work_dir, splitext_plus(basename(in_bam_fpath))[0] + '_markdup')
    safe_mkdir(dirname(tmp_fpath))
    cmdline = ('{bammarkduplicates} tmpfile={tmp_fpath} I={in_bam_fpath} O={out_bam_fpath}'
           ).format(**locals())
    res = call(cnf, cmdline, output_fpath=out_bam_fpath, stdout_to_outputfile=False, exit_on_error=False)
    if res:
        return out_bam_fpath
    else:
        return None


def bam_to_bed(cnf, bam_fpath, to_gzip=True):
    info('Converting the BAM to BED to save some memory.')  # from here: http://davetang.org/muse/2015/08/05/creating-a-coverage-plot-using-bedtools-and-r/
    bam_bed_fpath = splitext_plus(bam_fpath)[0] + ('.bed.gz' if to_gzip else '.bed')
    bedtools = get_system_path(cnf, 'bedtools')
    gzip = get_system_path(cnf, 'gzip')
    cmdline = '{bedtools} bamtobed -i {bam_fpath}'.format(**locals())
    cmdline += ' | {gzip}'.format(**locals()) if to_gzip else ''
    call(cnf, cmdline, output_fpath=bam_bed_fpath)
    return bam_bed_fpath


def bam_to_bed_nocnf(bam_fpath, bedtools='bedtools', gzip='gzip'):
    info('Converting the BAM to BED to save some memory.')  # from here: http://davetang.org/muse/2015/08/05/creating-a-coverage-plot-using-bedtools-and-r/
    bam_bed_fpath = splitext_plus(bam_fpath)[0] + '.bed.gz'
    cmdline = '{bedtools} bamtobed -i {bam_fpath} | {gzip} > {bam_bed_fpath}'.format(**locals())
    info(cmdline)
    os.system(cmdline)
    bam_bed_fpath = verify_file(bam_bed_fpath)
    if bam_bed_fpath:
        info('Done, saved to ' + bam_bed_fpath)
    else:
        err('Error, result is non-existent or empty')
    return bam_bed_fpath


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


def prepare_beds(cnf, features_bed=None, target_bed=None, seq2c_bed=None):
    if features_bed is None and target_bed is None:
        warn('No input target BED, and no features BED in the system config specified. Not making detailed per-gene reports.')
        # return None, None, None, None

    if target_bed:
        target_bed = verify_bed(target_bed, is_critical=True)

    if seq2c_bed:
        seq2c_bed = verify_bed(seq2c_bed, is_critical=True)

    if features_bed:
        features_bed = verify_bed(features_bed, is_critical=True)

    # if features_bed and target_bed and abspath(features_bed) == abspath(target_bed):
    #     warn('Same file used for exons and amplicons: ' + features_bed)

    # Features
    features_no_genes_bed = None
    if features_bed:
        # info()
        # info('Merging regions within genes...')
        # exons_bed = group_and_merge_regions_by_gene(cnf, exons_bed, keep_genes=True)
        #
        # info()
        # info('Sorting exons by (chrom, gene name, start)')
        # exons_bed = sort_bed(cnf, exons_bed)

        info()
        info('Filtering the features bed file to have only non-gene and no-transcript records...')
        features_no_genes_bed = intermediate_fname(cnf, features_bed, 'no_genes')
        call(cnf, 'grep -vw Gene ' + features_bed + ' | grep -vw Transcript', output_fpath=features_no_genes_bed)

    ori_target_bed_path = target_bed
    if target_bed:
        info()
        info('Remove comments in target...')
        target_bed = remove_comments(cnf, target_bed)

        info()
        info('Cut -f1,2,3,4 target...')
        target_bed = cut(cnf, target_bed, 4)

        info()
        info('Sorting target...')
        target_bed = sort_bed(cnf, target_bed)

        cols = count_bed_cols(target_bed)
        if True or cnf.reannotate or cols < 4:
            info()
            if not features_bed:
                critical(str(cols) + ' columns (less than 4), and no features to annotate regions '
                    '(please make sure you have set the "features" key in the corresponding genome section '
                    '(' + cnf.genome.name + ') in ' + cnf.sys_cnf)
            info('cnf.reannotate is ' + str(cnf.reannotate) + ', and cols in the target BED is ' + str(cols) +
                '. Annotating target with the gene names from the "features" file ' + features_bed + '...')
            target_bed = annotate_target(cnf, target_bed)

    def remove_no_anno(l, i):
        if l.split('\t')[3].strip() == '.': return None
        else: return l

    if not seq2c_bed and target_bed or seq2c_bed and seq2c_bed == ori_target_bed_path:
        info('Seq2C bed: remove regions with no gene annotation')
        seq2c_bed = target_bed
        seq2c_bed = iterate_file(cnf, seq2c_bed, remove_no_anno, suffix='filt')

    elif seq2c_bed:
        info()
        info('Remove comments in seq2c bed...')
        seq2c_bed = remove_comments(cnf, seq2c_bed)

        info()
        info('Sorting seq2c bed...')
        seq2c_bed = sort_bed(cnf, seq2c_bed)

        cols = count_bed_cols(seq2c_bed)
        if cols < 4:
            info()
            info('Number columns in SV bed is ' + str(cols) + '. Annotating amplicons with gene names...')
            seq2c_bed = annotate_target(cnf, seq2c_bed)
        elif 8 > cols > 4:
            seq2c_bed = cut(cnf, seq2c_bed, 4)
        elif cols > 8:
            seq2c_bed = cut(cnf, seq2c_bed, 8)
        info('Filtering non-annotated entries in seq2c bed')
        seq2c_bed = iterate_file(cnf, seq2c_bed, remove_no_anno, suffix='filt')

    else:
        seq2c_bed = verify_bed(cnf.genome.cds)

    if target_bed:
        info()
        # info('Merging amplicons...')
        # target_bed = group_and_merge_regions_by_gene(cnf, target_bed, keep_genes=False)

        info('Sorting target by (chrom, gene name, start)')
        target_bed = sort_bed(cnf, target_bed)

    return features_bed, features_no_genes_bed, target_bed, seq2c_bed


def extract_gene_names_and_filter_exons(cnf, target_bed, features_bed, features_no_genes_bed):
    gene_key_set = set()
    gene_key_list = []

    info()
    info('Getting gene list')

    # if genes_fpath:
    #     with open(genes_fpath) as f:
    #         gene_key_list = [g.strip() for g in f.read().split('\n') if g]
    #         gene_key_set = set(gene_key_list)
    #     info('Using genes from ' + genes_fpath + ', filtering exons and amplicons with this genes.')
    #     if target_bed:
    #         target_bed = filter_bed_with_gene_set(cnf, target_bed, gene_key_set)
    #     if exons_bed:
    #         exons_bed = filter_bed_with_gene_set(cnf, exons_bed, gene_key_set)
    #         exons_no_genes_bed = filter_bed_with_gene_set(cnf, exons_no_genes_bed, gene_key_set)
    # else:

    if target_bed:
        info()
        gene_key_set, gene_key_list = get_gene_keys(target_bed)
        info('Using genes from the amplicons list ' + target_bed)
        if features_bed and cnf.prep_bed is not False:
            info('Trying filtering exons with these ' + str(len(gene_key_list)) + ' genes.')
            features_filt_bed = filter_bed_with_gene_set(cnf, features_bed, gene_key_set, suffix='target_genes_1st_round')
            if not verify_file(features_filt_bed):
                info()
                warn('No gene symbols from the capture BED file was found in the features BED file. Re-annotating target...')
                target_bed = annotate_target(cnf, target_bed)
                #info('Merging regions within genes...')
                #target_bed = group_and_merge_regions_by_gene(cnf, target_bed, keep_genes=False)
                info('Sorting amplicons_bed by (chrom, gene_name, start)')
                target_bed = sort_bed(cnf, target_bed)
                info('Getting gene names again...')
                gene_key_set, gene_key_list = get_gene_keys(target_bed)
                info()
                info('Using genes from the new amplicons list, filtering features with this genes again.')
                features_filt_bed = filter_bed_with_gene_set(cnf, features_bed, gene_key_set, suffix='target_genes_2nd_round')
                if not verify_file(features_filt_bed):
                    critical('No gene symbols from the capture BED file was found in the features BED.')
            features_bed = features_filt_bed
            info('Filtering the full features file including gene records.')
            features_no_genes_bed = filter_bed_with_gene_set(cnf, features_no_genes_bed, gene_key_set, suffix='target_genes')
    elif features_no_genes_bed:
        info()
        info('No target (WGS), getting the gene names from the full features list...')
        gene_key_set, gene_key_list = get_gene_keys(features_no_genes_bed)
    info()

    return gene_key_set, gene_key_list, target_bed, features_bed, features_no_genes_bed


def calc_region_number(bed_fpath):
    with open(bed_fpath) as f:
        return sum(1 for l in f if l.strip() and not l.strip().startswith('#'))


def get_gene_keys(bed_fpath, chrom_index=0, gene_index=3):
    gene_keys_set = set()
    gene_keys_list = list()
    with open(bed_fpath) as f:
        for line in f:
            if not line or not line.strip() or line.startswith('#'):
                continue

            tokens = line.replace('\n', '').split('\t')
            if len(tokens) <= gene_index or len(tokens) <= chrom_index:
                continue

            chrom = tokens[chrom_index]
            for gn in tokens[gene_index].split(','):
                if (gn, chrom) not in gene_keys_set:
                    gene_keys_set.add((gn, chrom))
                    gene_keys_list.append((gn, chrom))

    return gene_keys_set, gene_keys_list


def annotate_target(cnf, target_bed):
    output_fpath = intermediate_fname(cnf, target_bed, 'ann')
    features_bed = verify_bed(cnf.genome.bed_annotation_features, is_critical=True, description='bed_annotation_features in system config')

    annotate_bed_py = get_system_path(cnf, 'python', join('tools', 'bed_processing', 'annotate_bed.py'))
    bedtools = get_system_path(cnf, 'bedtools')

    cmdline = '{annotate_bed_py} {target_bed} --work-dir {cnf.work_dir} --reference {features_bed} ' \
              '--genome {cnf.genome.name} --sys-cnf {cnf.sys_cnf} --run-cnf {cnf.run_cnf} ' \
              '-o {output_fpath}'.format(**locals())
    call(cnf, cmdline, output_fpath, stdout_to_outputfile=False)

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
    info('Doing some Seq2C specific preparation of the BED file...')

    cols = count_bed_cols(bed_fpath)

    seq2c_bed = None
    if 8 > cols > 4:
        seq2c_bed = cut(cnf, bed_fpath, 4)
    elif cols > 8:
        seq2c_bed = cut(cnf, bed_fpath, 8)
    else:
        seq2c_bed = bed_fpath

    if cols >= 4:
        # removing regions with no gene annotation
        def f(l, i):
            if l.split('\t')[3].strip() == '.': return None
            else: return l
        seq2c_bed = iterate_file(cnf, seq2c_bed, f, suffix='filt')

    info('Done: ' + seq2c_bed)
    return seq2c_bed


def filter_bed_with_gene_set(cnf, bed_fpath, gene_keys_set, suffix=None):
    def fn(l, i):
        if l:
            fs = l.split('\t')
            if len(fs) < 4:
                return None
            new_gns = []
            c = fs[0]
            for g in fs[3].split(','):
                if (g, c) in gene_keys_set:
                    new_gns.append(g)
            if new_gns:
                return l.replace(fs[3], ','.join(new_gns))

    return iterate_file(cnf, bed_fpath, fn, suffix=suffix or 'filt_genes', check_result=False)


def sort_bed(cnf, input_bed_fpath, output_bed_fpath=None):
    input_bed_fpath = verify_bed(input_bed_fpath)
    output_bed_fpath = adjust_path(output_bed_fpath) if output_bed_fpath else intermediate_fname(cnf, input_bed_fpath, 'sorted')

    class Region(SortableByChrom):
        def __init__(self, chrom, start, end, other_fields, chrom_ref_order):
            SortableByChrom.__init__(self, chrom, chrom_ref_order)
            self.start = start
            self.end = end
            self.chrom_ref_order = chrom_ref_order
            self.other_fields = tuple(other_fields)

        def get_key(self):
            return self.chrom_ref_order, self.start, self.end, self.other_fields

    regions = []
    chr_lengths = get_chr_lengths_from_seq(cnf.genome.seq)
    chr_order = {c: i for i, (c, l) in enumerate(chr_lengths)}

    info('Sorting regions in ' + input_bed_fpath)
    if cnf.reuse_intermediate and isfile(output_bed_fpath) and verify_bed(output_bed_fpath):
        info(output_bed_fpath + ' exists, reusing')
        return output_bed_fpath

    with open(input_bed_fpath) as f:
        with file_transaction(cnf.work_dir, output_bed_fpath) as tx:
            with open(tx, 'w') as out:
                for l in f:
                    if not l.strip():
                        continue
                    if l.strip().startswith('#'):
                        out.write(l)
                        continue

                    fs = l.strip().split('\t')
                    chrom = fs[0]
                    start = int(fs[1])
                    end = int(fs[2])
                    other_fields = fs[3:]
                    order = chr_order.get(chrom, -1)
                    regions.append(Region(chrom, start, end, other_fields, order))

                for region in sorted(regions, key=lambda r: r.get_key()):
                    fs = [region.chrom, str(region.start), str(region.end)]
                    fs.extend(region.other_fields)
                    out.write('\t'.join(fs) + '\n')

    info('Sorted ' + str(len(regions)) + ' regions, saved to ' + output_bed_fpath + '\n')
    return output_bed_fpath


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


def fix_bed_for_qualimap(bed_fpath, qualimap_bed_fpath):
    with open(qualimap_bed_fpath, 'w') as out, open(bed_fpath) as inn:
        for l in inn:
            fields = l.strip().split('\t')

            if len(fields) < 3:
                continue
            try:
                int(fields[1]), int(fields[2])
            except ValueError:
                continue

            if len(fields) < 4:
                fields.append('.')

            if len(fields) < 5:
                fields.append('0')

            if len(fields) < 6:
                fields.append('+')

            out.write('\t'.join(fields) + '\n')


def call_sambamba(cnf, cmdl, bam_fpath, output_fpath=None, sambamba=None, use_grid=False,
                  command_name='', sample_name=None):
    sambamba = sambamba or get_system_path(cnf, 'sambamba', is_critical=True)
    sample_name = sample_name or basename(bam_fpath).split('.')[0]
    if use_grid:
        grid_sambabma = get_script_cmdline(cnf, 'python', join('tools', 'bed_processing', 'sambamba.py'))
        cmdl = cmdl.replace(' "', ' \'\"__QUOTE__')
        cmdl = cmdl.replace('" ', '__QUOTE__\"\' ')
        grid_cmdl = grid_sambabma + ' ' + bam_fpath + ' ' + sambamba + ' ' + cmdl
        job_name = command_name + '_' + sample_name
        j = submit_job(cnf, grid_cmdl, job_name=job_name, output_fpath=output_fpath)
        info()
        return j
    else:
        index_bam(cnf, bam_fpath, sambamba=sambamba)
        cmdl = sambamba + ' ' + cmdl
        stderr_dump = []
        res = call(cnf, cmdl, output_fpath=output_fpath, exit_on_error=False, stderr_dump=stderr_dump)
        if not res:
            for l in stderr_dump:
                if 'sambamba-view: BAM index file (.bai) must be provided' in l:
                    if isfile(isfile(bam_fpath + '.bai')):
                        info('Removing .bai and re-indexing...')
                        os.remove(bam_fpath + '.bai')
                    index_bam(cnf, bam_fpath, sambamba)
                    res = call(cnf, cmdl, output_fpath=output_fpath)
        return res


def sambamba_depth(cnf, bed, bam, output_fpath=None, use_grid=False, depth_thresholds=None, sample_name=None):
    sample_name = sample_name or splitext_plus(basename(bam))[0]

    if not output_fpath:
        output_fpath = join(cnf.work_dir,
            splitext_plus(basename(bed))[0] + '_' + sample_name + '_sambamba_depth.txt')

    if cnf.reuse_intermediate and verify_file(output_fpath, silent=True):
        info(output_fpath + ' exists, reusing.')
        if use_grid:
            return None
        else:
            return output_fpath
    sambamba = get_system_path(cnf, join(get_ext_tools_dirname(), 'sambamba'), is_critical=True)
    depth_thresholds = depth_thresholds or cnf.coverage_reports.depth_thresholds
    thresholds_str = ' -T'.join([str(d) for d in depth_thresholds])
    cmdline = 'depth region -F "not duplicate and not failed_quality_control" -L {bed} -T {thresholds_str} {bam}'.format(**locals())

    return call_sambamba(cnf, cmdline, output_fpath=output_fpath, bam_fpath=bam,
        sambamba=sambamba, use_grid=use_grid, command_name='depth_' + splitext_plus(basename(bed))[0],
                         sample_name=sample_name)


def remove_dups(cnf, bam, output_fpath, sambamba=None, use_grid=False):
    cmdline = 'view --format=bam -F "not duplicate" {bam}'.format(**locals())  # -F (=not) 1024 (=duplicate)
    sambamba = sambamba or get_system_path(cnf, 'sambamba', is_critical=True)
    # if use_grid:
    #     j = submit_job(cnf, sambamba + ' ' + cmdline,
    #                    'DEDUP__' + cnf.project_name + '__' + basename(bam).split('.')[0],
    #                    output_fpath=output_fpath)
    #     info()
    #     return j
    # else:
    return call_sambamba(cnf, cmdline, output_fpath=output_fpath, bam_fpath=bam,
                         sambamba=sambamba, use_grid=use_grid, command_name='not_duplicate')


def remove_dups_picard(cnf, bam_fpath):
    picard = get_system_path(cnf, 'java', 'picard')
    if not picard:
        critical('No picard in the system')

    info('Running picard dedup for "' + basename(bam_fpath) + '"')

    dup_metrics_txt = join(cnf.work_dir, 'picard_dup_metrics.txt')
    output_fpath = intermediate_fname(cnf, bam_fpath, 'pcd_dedup')

    cmdline = '{picard} MarkDuplicates' \
              ' I={bam_fpath}' \
              ' O={output_fpath}' \
              ' METRICS_FILE={dup_metrics_txt}' \
              ' REMOVE_DUPLICATES=True' \
              ' VALIDATION_STRINGENCY=LENIENT'
    res = call(cnf, cmdline.format(**locals()), output_fpath=output_fpath,
        stdout_to_outputfile=False, exit_on_error=False)

    if res != output_fpath:  # error occurred, try to correct BAM and restart
        warn('Picard deduplication failed for "' + basename(bam_fpath) + '". Fixing BAM and restarting Picard...')
        bam_fpath = _fix_bam_for_picard(cnf, bam_fpath)
        res = call(cnf, cmdline.format(**locals()), output_fpath=output_fpath,
            stdout_to_outputfile=False, exit_on_error=False)

    if res == output_fpath:
        dup_rate = _parse_picard_dup_report(dup_metrics_txt)
        assert dup_rate <= 1.0 or dup_rate is None, str(dup_rate)
        info('Duplication rate (picard): ' + str(dup_rate))
        return output_fpath
    else:
        return None


def _fix_bam_for_picard(cnf, bam_fpath):
    def __process_problem_read_aligns(read_aligns):
        # each alignment: 0:NAME 1:FLAG 2:CHR 3:COORD 4:MAPQUAL 5:CIGAR 6:MATE_CHR 7:MATE_COORD TLEN SEQ ...
        def __get_key(align):
            return align.split('\t')[2] + '@' + align.split('\t')[3]

        def __get_mate_key(align):
            return (align.split('\t')[6] if align.split('\t')[2] != '=' else align.split('\t')[2]) \
                   + '@' + align.split('\t')[7]

        chr_coord = OrderedDict()
        for align in read_aligns:
            key = __get_key(align)
            if key not in chr_coord:
                chr_coord[key] = []
            chr_coord[key].append(align)
        correct_pairs = []
        for align in read_aligns:
            mate_key = __get_mate_key(align)
            if mate_key in chr_coord:
                for pair_align in chr_coord[mate_key]:
                    if read_aligns.index(pair_align) <= read_aligns.index(align):
                        continue
                    if __get_mate_key(pair_align) == __get_key(align):
                        correct_pairs.append((align, pair_align))
        if not correct_pairs:
            return []
        if len(correct_pairs) > 1:
            # sort by sum of mapping quality of both alignments
            correct_pairs.sort(key=lambda pair: pair[0].split('\t')[4] + pair[1].split('\t')[4], reverse=True)
        return [correct_pairs[0][0], correct_pairs[0][1]]

    samtools = get_system_path(cnf, 'samtools')
    try:
        import pysam
        without_pysam = False
    except ImportError:
        without_pysam = True

    # find reads presented more than twice in input BAM
    if without_pysam:
        qname_sorted_sam_fpath = intermediate_fname(cnf, bam_fpath, 'qname_sorted')[:-len('bam')] + 'sam'
        # queryname sorting; output is SAM
        cmdline = '{samtools} view {bam_fpath} | sort '.format(**locals())
        call(cnf, cmdline, qname_sorted_sam_fpath)
        qname_sorted_file = open(qname_sorted_sam_fpath, 'r')
    else:
        qname_sorted_bam_fpath = intermediate_fname(cnf, bam_fpath, 'qname_sorted')
        # queryname sorting (-n), to stdout (-o), 'prefix' is not used; output is BAM
        cmdline = '{samtools} sort -n -o {bam_fpath} prefix'.format(**locals())
        call(cnf, cmdline, qname_sorted_bam_fpath)
        qname_sorted_file = pysam.Samfile(qname_sorted_bam_fpath, 'rb')
    problem_reads = dict()
    cur_read_aligns = []
    for line in qname_sorted_file:
        line = str(line)
        if cur_read_aligns:
            if line.split('\t')[0] != cur_read_aligns[0].split('\t')[0]:
                if len(cur_read_aligns) > 2:
                    problem_reads[cur_read_aligns[0].split('\t')[0]] = cur_read_aligns
                cur_read_aligns = []
        flag = int(line.split('\t')[1])
        cur_read_aligns.append(line)
    if len(cur_read_aligns) > 2:
        problem_reads[cur_read_aligns[0].split('\t')[0]] = cur_read_aligns
    qname_sorted_file.close()

    for read_id, read_aligns in problem_reads.items():
        problem_reads[read_id] = __process_problem_read_aligns(read_aligns)

    # correct input BAM
    fixed_bam_fpath = intermediate_fname(cnf, bam_fpath, 'fixed_for_picard')
    fixed_sam_fpath = fixed_bam_fpath[:-len('bam')] + 'sam'
    if without_pysam:
        sam_fpath = intermediate_fname(cnf, bam_fpath, 'tmp')[:-len('bam')] + 'sam'
        cmdline = '{samtools} view -h {bam_fpath}'.format(**locals())
        call(cnf, cmdline, sam_fpath)
        input_file = open(sam_fpath, 'r')
        fixed_file = open(fixed_sam_fpath, 'w')
    else:
        input_file = pysam.Samfile(bam_fpath, 'rb')
        fixed_file = pysam.Samfile(fixed_bam_fpath, 'wb', template=input_file)
    for line in input_file:
        if without_pysam and line.startswith('@'):  # header
            fixed_file.write(line)
            continue
        read_name = str(line).split('\t')[0]
        if read_name in problem_reads and str(line) not in problem_reads[read_name]:
            continue
        fixed_file.write(line)
    input_file.close()
    fixed_file.close()
    if without_pysam:
        cmdline = '{samtools} view -bS {fixed_sam_fpath}'.format(**locals())
        call(cnf, cmdline, fixed_bam_fpath)

    return fixed_bam_fpath


def _parse_picard_dup_report(dup_report_fpath):
    with open(dup_report_fpath) as f:
        for l in f:
            if l.startswith('## METRICS CLASS'):
                l_NEXT = None
                ind = None
                try:
                    l_LIBRARY = next(f)
                    if l_LIBRARY.startswith('LIBRARY'):
                        ind = l_LIBRARY.strip().split().index('PERCENT_DUPLICATION')
                        l_NEXT = next(f)
                        while l_NEXT.startswith(' ') or l_NEXT.startswith('\t'):
                            l_NEXT = next(f)
                except StopIteration:
                    pass
                else:
                    if l_NEXT and ind:
                        fields = l_NEXT.split()
                        if fields[0] == 'Unknown':
                            ind += 1
                        if len(fields) > ind:
                            dup_rate = 1.0 * float(fields[ind])
                            return dup_rate
    err('Error: cannot read duplication rate from ' + dup_report_fpath)


def count_in_bam(cnf, bam, query, dedup=False, bed=None, use_grid=False, sample_name=None):
    if dedup:
        query += ' and not duplicate'
    name = 'num_' + (query.replace(' ', '_') or 'reads')
    if bed:
        name += '_on_target_' + basename(bed)
    sample_name = sample_name or basename(bam)
    output_fpath = join(cnf.work_dir, sample_name + '_' + name)

    cmdline = 'view -c -F "{query}" {bam}'.format(**locals())
    if bed:
        cmdline += ' -L ' + bed

    res = call_sambamba(cnf, cmdline, output_fpath=output_fpath, bam_fpath=bam,
                        use_grid=use_grid, command_name=name, sample_name=sample_name)
    if not use_grid:
        with open(output_fpath) as f:
            return int(f.read().strip())
    else:
        return res  # Job object


def number_of_reads(cnf, bam, dedup=False, use_grid=False, sample_name=None):
    return count_in_bam(cnf, bam, '', dedup, use_grid=use_grid, sample_name=sample_name)


def number_of_mapped_reads(cnf, bam, dedup=False, use_grid=False, sample_name=None):
    return count_in_bam(cnf, bam, 'not unmapped', dedup, use_grid=use_grid, sample_name=sample_name)


def number_of_properly_paired_reads(cnf, bam, dedup=False, use_grid=False, sample_name=None):
    return count_in_bam(cnf, bam, 'proper_pair', dedup, use_grid=use_grid, sample_name=sample_name)


def number_of_dup_reads(cnf, bam, use_grid=False, sample_name=None):
    return count_in_bam(cnf, bam, 'duplicate', use_grid=use_grid, sample_name=sample_name)


def number_mapped_reads_on_target(cnf, bed, bam, dedup=False, use_grid=False, sample_name=None):
    return count_in_bam(cnf, bam, 'not unmapped', dedup, bed=bed, use_grid=use_grid, sample_name=sample_name)


def flag_stat(cnf, bam):
    output_fpath = join(cnf.work_dir, basename(bam) + '_flag_stats')
    cmdline = 'flagstat {bam}'.format(**locals())
    call_sambamba(cnf, cmdline, output_fpath=output_fpath, bam_fpath=bam, command_name='flagstat')
    stats = dict()
    with open(output_fpath) as f:
        lines = f.readlines()
        for stat, fun in [('total', number_of_reads),
                          ('duplicates', number_of_dup_reads),  # '-f 1024'
                          ('mapped', number_of_mapped_reads),   # '-F 4'
                          ('properly paired', number_of_properly_paired_reads)]:  # '-f 2'
            try:
                val = next(l.split()[0] for l in lines if stat in l)
            except StopIteration:
                warn('Cannot extract ' + stat + ' from flagstat output ' + output_fpath + '. Trying samtools view -c...')
                val = None
            else:
                try:
                    val = int(val)
                except ValueError:
                    warn('Cannot parse value ' + str(val) + ' from ' + stat + ' from flagstat output ' + output_fpath + '. Trying samtools view -c...')
                    val = None
            if val is not None:
                stats[stat] = val
            else:
                stats[stat] = fun(cnf, bam)
    return stats


def get_padded_bed_file(cnf, bed, genome, padding):
    info('Making bed file for padded regions...')
    bedtools = get_system_path(cnf, 'bedtools')
    cmdline = '{bedtools} slop -i {bed} -g {genome} -b {padding}'.format(**locals())
    output_fpath = intermediate_fname(cnf, bed, 'padded')
    call(cnf, cmdline, output_fpath)
    return output_fpath


def intersect_bed(cnf, bed1, bed2):
    bed1_fname, _ = splitext_plus(basename(bed1))
    bed2_fname, _ = splitext_plus(basename(bed2))
    output_fpath = join(cnf['work_dir'], bed1_fname + '__' + bed2_fname + '.bed')
    bedtools = get_system_path(cnf, 'bedtools')
    cmdline = '{bedtools} intersect -u -a {bed1} -b {bed2}'.format(**locals())
    call(cnf, cmdline, output_fpath, verify_output_not_empty=False)
    return output_fpath


def verify_bam(fpath, description='', is_critical=False, silent=False):
    if not verify_file(fpath, description, is_critical=is_critical, silent=silent):
        return None

    fpath = adjust_path(fpath)

    logfn = critical if is_critical else err
    if not fpath.endswith('.bam'):
        logfn('The file ' + fpath + ' is supposed to be BAM but does not have the .bam '
            'extension. Please, make sure you pass proper file.')
        return None

    textchars = ''.join(map(chr, [7, 8, 9, 10, 12, 13, 27] + range(0x20, 0x100)))
    is_binary_string = lambda baitiki: bool(baitiki.translate(None, textchars))
    if not is_binary_string(open(fpath).read(3)):
        logfn('The BAM file ' + fpath + ' must be a binary file.')
        return None

    return fpath


def verify_bed(fpath, description='', is_critical=False, silent=False):
    if not verify_file(fpath, description, is_critical=is_critical, silent=silent):
        return None

    fpath = adjust_path(fpath)

    error = BedFile(fpath).checkformat()
    if error:
        fn = critical if is_critical else err
        fn('Error: incorrect bed file format (' + fpath + '): ' + str(error) + '\n')
        return None

    return fpath


def check_md5(work_dir, fpath, file_type, silent=False):
    md5_fpath = join(work_dir, file_type + '_md5.txt')
    new_md5 = md5(fpath)
    info('md5 of ' + fpath + ' is ' + str(new_md5))
    prev_md5 = None
    if isfile(md5_fpath):
        with open(md5_fpath) as f:
            prev_md5 = f.read()
    else:
        info('Previous md5 file ' + md5_fpath + ' does not exist')
    info('Checking previous md5 from ' + md5_fpath + ': ' + str(prev_md5))

    if prev_md5 == new_md5:
        if not silent:
            info('Reusing previous ' + file_type.upper() + ' files.')
        return True
    else:
        if not silent:
            info('Pre-processing input ' + file_type.upper() + ' file')
        if prev_md5:
            if not silent:
                info('Prev ' + file_type.upper() + ' md5: ' + str(prev_md5))
                info('New ' + file_type.upper() + ' md5: ' + str(new_md5))

        with open(md5_fpath, 'w') as f:
            f.write(str(new_md5))
        return False


class BedFile:
    def __init__(self, _filename):
        self.filename = _filename
        self.chrs = None
        self.nregions = None

    def checkformat(self):
        """************************************************************************************************************************************************************
        Task: checks the format of the bed file. The only requirements checked are that each line presents at least 3 tab separated columns, the
            two on the right must present integer values indicating the start/end position respectively. Right value must be greater than the
            left value.
        Outputs:
            err: string containing the detected error. Empty string in case of a correct format.
        ************************************************************************************************************************************************************"""

        fd = file(self.filename)

        line = fd.readline()
        fields = line.split('\t')
        lc = 1
        err = ''

        # Checks that the two columns on the right contain integer values
        try:
            # Parses each line and checks that there are at least 3 fields, the two on the right containing integer values and being the right one
            # greater than the left one
            while (line <> '' and len(fields) > 2 and int(fields[1]) < int(fields[2])):
                lc += 1
                line = fd.readline()
                fields = line.split('\t')
        except ValueError:
            err += 'Incorrect start/end values at line ' + str(lc) + '\n'
            err += 'Start/End coordinates must be indicated with integer values. The right value must be greater than the left value.\n'
            err += 'Line found: ' + line
            fd.close()

            return err

        # If it get to this point means that either the file ended or there is a line with less than 3 fields
        if (line <> ''):
            err += 'Incorrect line format at line ' + str(lc) + '\n'
            err += 'At least three columns are expected in each line\n'
            err += 'The right value must be greater than the left value.\n'
            err += 'Line found: ' + line
            fd.close()

        return err

    def count_lines(self, filename=None):
        if filename is None:
            filename = self.filename
        return len(open(filename).readlines())