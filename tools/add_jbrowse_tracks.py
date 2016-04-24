#!/usr/bin/env python
import os
import traceback
from genericpath import isfile
from os.path import join, basename, splitext, isdir

from source.calling_process import call
from source.file_utils import verify_dir, safe_mkdir, verify_file, splitext_plus
from source.logger import warn, info
from source.targetcov.bam_and_bed_utils import index_bam
from source.utils import is_uk
from ext_modules.joblib import Parallel, delayed


def add_project_files_to_jbrowse(cnf, bcbio_structure):
    genome = cnf.genome.name
    jbrowse_data_path, _, _ = set_folders(genome)

    jbrowse_dirpath = join(jbrowse_data_path, 'tracks')
    jbrowse_project_dirpath = join(jbrowse_dirpath, bcbio_structure.project_name)

    safe_mkdir(jbrowse_project_dirpath)
    jbrowse_tracks_fpath = join(jbrowse_data_path, 'tracks.conf')

    vcf_fpath_by_sample = None
    caller = bcbio_structure.variant_callers.get('vardict') or \
             bcbio_structure.variant_callers.get('vardict-java')
    if caller:
        vcf_fpath_by_sample = caller.get_filt_vcf_by_sample()

    for sample in bcbio_structure.samples:
        if sample.bam:
            index_bam(cnf, sample.bam, use_grid=True)

    for sample in bcbio_structure.samples:
        if all(isfile(join(jbrowse_project_dirpath, sample.name + ext)) for ext in ['.bam', '.bam.bai', '.vcf.gz', '.vcf.gz.tbi', '.bigwig']):
            continue
        vcf_link = None
        if vcf_fpath_by_sample:
            vcf_fpath = vcf_fpath_by_sample[sample.name] if sample.name in vcf_fpath_by_sample else None
            if vcf_fpath and verify_file(vcf_fpath):
                vcf_link = create_jbrowse_symlink(genome, bcbio_structure.project_name, sample.name, vcf_fpath)
                if not verify_file(vcf_fpath + '.tbi'):
                    cmdline = '{tabix} {vcf_fpath}'.format(**locals())
                    call(cnf, cmdline, exit_on_error=False)
                create_jbrowse_symlink(genome, bcbio_structure.project_name, sample.name, vcf_fpath + '.tbi')

        if sample.bam:
            bam_link = create_jbrowse_symlink(genome, bcbio_structure.project_name, sample.name, sample.bam)
            create_jbrowse_symlink(genome, bcbio_structure.project_name, sample.name, sample.bam + '.bai')
            bigwig_link = create_jbrowse_symlink(genome, bcbio_structure.project_name, sample.name, splitext(sample.bam)[0] + '.bigwig')
            print_sample_tracks_info(sample.name, bcbio_structure.project_name, trunc_symlink(bam_link),
                                     trunc_symlink(bigwig_link), trunc_symlink(vcf_link), jbrowse_tracks_fpath)


def print_sample_tracks_info(sample, project_name, bam_link, bigwig_link, vcf_link, jbrowse_tracks_fpath):
    with open(jbrowse_tracks_fpath, 'a') as tracks:
        print >> tracks, '\n[ tracks.{sample} ]\n' \
                         '\nstoreClass     = JBrowse/Store/SeqFeature/BAM' \
                         '\nurlTemplate    = {bam_link}' \
                         '\nbaiUrlTemplate = {bam_link}.bai' \
                         '\nchunkSizeLimit = 100000000' \
                         '\nmaxHeight      = 10000' \
                         '\ncategory = {project_name}' \
                         '\ntype = JBrowse/View/Track/Alignments2' \
                         '\nkey  = {sample}\n'.format(**locals())
        print >> tracks, '\n[ tracks.{sample}_cov ]\n' \
                         '\nstoreClass     = JBrowse/Store/SeqFeature/BAM' \
                         '\nurlTemplate    = {bam_link}' \
                         '\nbaiUrlTemplate = {bam_link}.bai' \
                         '\nchunkSizeLimit = 100000000' \
                         '\ncategory = {project_name}' \
                         '\ntype = SNPCoverage' \
                         '\nkey  = {sample}_coverage_bam\n'.format(**locals())
        print >> tracks, '\n[ tracks.{sample}_bigwig ]\n' \
                         '\nstoreClass     = JBrowse/Store/SeqFeature/BigWig' \
                         '\nurlTemplate    = {bigwig_link}' \
                         '\ncategory = {project_name}' \
                         '\ntype = JBrowse/View/Track/Wiggle/XYPlot' \
                         '\nautoscale = local' \
                         '\nkey  = {sample}_coverage\n'.format(**locals())
        if vcf_link:
            print >> tracks, '\n[ tracks.{sample}_vcf ]\n' \
                         '\nstoreClass     = JBrowse/Store/SeqFeature/VCFTabix' \
                         '\nurlTemplate    = {vcf_link}' \
                         '\ncategory = {project_name}' \
                         '\ntype = JBrowse/View/Track/CanvasVariants' \
                         '\nkey  = {sample}_variants\n'.format(**locals())

def trunc_symlink(link):
    if not link:
        return None
    return 'tracks' + link.split('tracks', 1)[1]


def set_folders(genome):
    jbrowse_basepath = '/home/klpf990/public_html/JBrowse-1.11.6'
    jbrowse_browser_path = 'http://blue.usbod.astrazeneca.net/~klpf990/JBrowse-1.11.6'
    if is_uk():
        jbrowse_basepath = '/ngs/web_content/reports/JBrowse'
        jbrowse_browser_path = 'http://ukapdlnx115.ukapd.astrazeneca.net/ngs/reports/JBrowse'
    data_dirname = 'data_hg19'
    if 'hg38' in genome:
        data_dirname = 'data_hg38'
    elif 'hg' not in genome:
        return None, None, None
    jbrowse_data_path = join(jbrowse_basepath, data_dirname)
    return jbrowse_data_path, data_dirname, jbrowse_browser_path


def get_jbrowser_link(genome, sample, bed_fname=None):
    jbrowse_data_path, data_dirname, jbrowse_browser_path = set_folders(genome)
    bed = ''
    if bed_fname:
        bed = ',' + bed_fname
    return '{jbrowse_browser_path}/?data={data_dirname}&tracks=DNA,' \
           '{sample}_bigwig,{sample}_vcf{bed},{sample}&highlight='.format(**locals())


def create_jbrowse_symlink(genome, project_name, sample, file_fpath):
    jbrowse_data_path, _, _ = set_folders(genome)
    jbrowse_dirpath = join(jbrowse_data_path, 'tracks')
    jbrowse_project_dirpath = join(jbrowse_dirpath, project_name)
    base, ext = splitext_plus(file_fpath)
    if ext in ['.tbi', '.bai']:
        base, ext2 = splitext_plus(base)
        ext = ext2 + ext
    sym_link = join(jbrowse_project_dirpath, sample + ext)
    if not verify_dir(jbrowse_project_dirpath):
        safe_mkdir(jbrowse_project_dirpath)
    if isfile(file_fpath) and not isfile(sym_link):
        try:
            os.symlink(file_fpath, sym_link)
        except OSError:
            warn(traceback.format_exc())
    if isfile(sym_link):
        change_permissions(sym_link)
    return sym_link



def change_permissions(path):
    try:
        os.system('chmod -R g+w ' + path)
    except:
        pass
