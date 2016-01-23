from os.path import join, splitext
from collections import OrderedDict
from source.file_utils import verify_file, add_suffix
from source.logger import info
from source import targetcov


project_kind_by_prefix = {
    'Bio_': 'bioscience',
    'Dev_': 'dev',
    'EXT_': 'external',
    'TS_': 'translation',
}

detail_gene_report_baseending = '.details.gene'

varfilter_name           = 'varFilter'
varannotate_name         = 'varAnnotate'
targetseq_name           = 'targetSeq'
flag_regions_name        = 'flaggedRegions'
cnv_dir                  = 'cnv'

varqc_name               = 'varQC'
varqc_after_name         = 'varQC_postVarFilter'
ngscat_name              = 'ngscat'
qualimap_name            = 'qualimap'
picard_name              = 'picard'
targqc_name              = 'targQC'
fastqc_name              = 'fastqc'

varqc_repr               = 'Var QC'
varqc_after_repr         = 'Var QC after filtering'
ngscat_repr              = 'Ngscat'
qualimap_repr            = 'Qualimap'
targqc_repr              = 'Target QC'
fastqc_repr              = 'FastQC'

clinreport_name = clinreport_dir = 'oncologyNGSReport'

seq2c_name               = 'Seq2C'
seq2c_seq2cov_ending     = 'seq2c_seq2cov.txt'

dedup_bam                = 'dedup'  # work/post_processing/dedup and -ready.dedup.bam

ngscat_report_fname             = 'captureQC.html'
qualimap_report_fname           = 'qualimapReport.html'
qualimap_genome_results_fname   = 'genome_results.txt'
qualimap_raw_data_dirname       = 'raw_data_qualimapReport'
qualimap_ishist_fsubpath        = join(qualimap_raw_data_dirname, 'insert_size_histogram.txt')
qualimap_covhist_fsubpath       = join(qualimap_raw_data_dirname, 'coverage_histogram.txt')

picard_ishist_fname = 'picard_ishist.txt'
fastqc_report_fname = 'fastqc_report.html'

mut_file_ext = 'txt'
mut_fname_template = '{caller_name}.' + mut_file_ext
mut_single_suffix = 'single'
mut_paired_suffix = 'paired'
mut_pass_suffix = 'PASS'


class BaseSample:
    def __init__(self, name, dirpath, bam=None, bed=None, vcf=None, genome=None,
                 targqc_dirpath=None, ngscat_dirpath=None, qualimap_dirpath=None,
                 fastqc_dirpath=None, picard_dirpath=None, clinical_report_dirpath=None,
                 flagged_regions_dirpath=None, normal_match=None, sv_fpath=None):
        self.name = name
        self.bam = bam
        self.dedup_bam = None
        self.bed = bed
        self.qualimap_bed = None
        self.vcf_by_callername = OrderedDict()  # string -> vcf_fpath
        self.vcf = vcf
        self.dirpath = dirpath
        self.phenotype = None
        self.gender = None
        self.genome = None
        self.var_dirpath = None
        self.normal_match = normal_match
        self.min_af = None
        self.sv_fpath = sv_fpath

        self.targqc_dirpath                  = None
        self.targetcov_html_fpath            = None
        self.targetcov_json_fpath            = None
        self.targetcov_detailed_txt          = None
        self.targetcov_detailed_tsv          = None
        self.targetcov_norm_depth_vcf_txt    = None
        self.targetcov_norm_depth_vcf_tsv    = None

        self.flagged_regions_dirpath         = None
        self.flagged_tsv                     = None
        self.flagged_txt                     = None

        self.clinical_report_dirpath         = None
        self.clinical_targqc_tsv             = None
        self.clinical_mutation_tsv           = None
        self.clinical_target_tsv             = None
        self.clinical_html                   = None

        self.ngscat_dirpath                  = None
        self.ngscat_html_fpath               = None

        self.qualimap_dirpath                = None
        self.qualimap_html_fpath             = None
        self.qualimap_genome_results_fpath   = None
        self.qualimap_ins_size_hist_fpath    = None
        self.qualimap_cov_hist_fpath         = None

        self.fastqc_dirpath                  = None
        self.fastqc_html_fpath               = None

        self.picard_dirpath                  = None
        self.picard_ins_size_hist_txt_fpath  = None
        self.picard_ins_size_hist_pdf_fpath  = None

        if targqc_dirpath:
            self.targqc_dirpath = targqc_dirpath
            self.targetcov_html_fpath           = join(self.targqc_dirpath, name + '.' + targetseq_name +  '.html')
            self.targetcov_json_fpath           = join(self.targqc_dirpath, name + '.' + targetseq_name +  '.json')
            self.targetcov_detailed_txt         = join(self.targqc_dirpath, name + '.' + targetseq_name +  detail_gene_report_baseending + '.txt')
            self.targetcov_detailed_tsv         = join(self.targqc_dirpath, name + '.' + targetseq_name +  detail_gene_report_baseending + '.tsv')
            self.targetcov_norm_depth_vcf_txt   = None
            self.targetcov_norm_depth_vcf_tsv   = None

        if flagged_regions_dirpath:
            self.flagged_regions_dirpath        = flagged_regions_dirpath
            self.flagged_tsv                    = join(self.flagged_regions_dirpath, name + '.' + flag_regions_name + '.tsv')
            self.flagged_txt                    = join(self.flagged_regions_dirpath, name + '.' + flag_regions_name + '.txt')

        if clinical_report_dirpath:
            self.clinical_report_dirpath        = clinical_report_dirpath
            self.clinical_targqc_tsv            = join(self.clinical_report_dirpath, name + '.coverage.tsv')
            self.clinical_mutation_tsv          = join(self.clinical_report_dirpath, name + '.mutations.tsv')
            self.clinical_target_tsv            = join(self.clinical_report_dirpath, name + '.target.tsv')
            self.clinical_html                  = join(self.clinical_report_dirpath, name + '.html')

        if ngscat_dirpath:
            self.ngscat_dirpath                 = ngscat_dirpath
            self.ngscat_html_fpath              = join(self.ngscat_dirpath, ngscat_report_fname)
 
        if qualimap_dirpath:
            self.qualimap_dirpath               = qualimap_dirpath
            self.qualimap_html_fpath            = join(self.qualimap_dirpath, qualimap_report_fname)
            self.qualimap_genome_results_fpath  = join(self.qualimap_dirpath, qualimap_report_fname)
            self.qualimap_ins_size_hist_fpath   = join(self.qualimap_dirpath, qualimap_ishist_fsubpath)
            self.qualimap_cov_hist_fpath        = join(self.qualimap_dirpath, qualimap_covhist_fsubpath)

        if fastqc_dirpath:
            self.fastqc_dirpath                 = fastqc_dirpath
            self.fastqc_html_fpath              = join(self.fastqc_dirpath, fastqc_report_fname)

        if picard_dirpath:
            self.picard_dirpath                 = picard_dirpath
            self.picard_ins_size_hist_txt_fpath = join(self.picard_dirpath, picard_ishist_fname)
            self.picard_ins_size_hist_pdf_fpath = join(self.picard_dirpath, splitext(picard_ishist_fname)[0] + '.pdf')

        # # Only for Bcbio and TargQC:
        # self.ngscat_html_fpath                 = self.make_fpath(path_base + 'captureQC.html', name=ngscat_name)
        # self.qualimap_html_fpath               = self.make_fpath(path_base + 'qualimapReport.html', name=qualimap_name)
        # self.qualimap_genome_results_fpath     = self.make_fpath(path_base + 'genome_results.txt', name=qualimap_name)
        # self.qualimap_ins_size_hist_fpath      = self.make_fpath(path_base + 'raw_data_qualimapReport/insert_size_histogram.txt', name=qualimap_name)
        # self.qualimap_cov_hist_fpath           = self.make_fpath(path_base + 'raw_data_qualimapReport/coverage_histogram.txt', name=qualimap_name)
        # # self.qualimap_postproc_cov_hist_fpath  = self.make_fpath(path_base + 'raw_data_qualimapReport/coverage_histogram.txt', name=qualimap_name)
        # self.fastqc_html_fpath                 = self.make_fpath(path_base + 'fastqc_report.html', name=fastqc_name)
        # # self.seq2cov_output_fpath              = self.make_fpath(path_base + '{sample}.{name}_' + seq2c_seq2cov_ending, name=targetseq_name)

    # def make_fpath(self, path_template, **kwargs):
    #     keys = dict(dict(dirpath=self.dirpath, sample=self.name).items() + kwargs.items())
    #     return join(*path_template.split('/')).format(**keys)

    def __cmp__(self, other):
        return cmp(self.key_to_sort(), other.key_to_sort())

    def key_to_sort(self):
        parts = []

        cur_part = []
        prev_was_num = False

        for c in self.name:
            if prev_was_num == c.isdigit() and c not in ['-', '.']:  # same type of symbol, but not - or .
                cur_part.append(c)
            else:
                if cur_part:
                    part = ''.join(cur_part)
                    if prev_was_num:
                        part = int(part)
                    parts.append(part)
                    cur_part = []

                if c in ['-', '.']:
                    pass
                else:
                    if c.isdigit():
                        prev_was_num = True
                    else:
                        prev_was_num = False
                    cur_part.append(c)
        if cur_part:
            part = ''.join(cur_part)
            if prev_was_num:
                part = int(part)
            parts.append(part)

        return tuple(parts)

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name


class TargQC_Sample(BaseSample):
    def __init__(self, name, dirpath, bam=None, bed=None, vcf=None, genome=None):
        BaseSample.__init__(self, name, dirpath,
            targqc_dirpath=dirpath,
            ngscat_dirpath=join(dirpath, ngscat_name),
            qualimap_dirpath=join(dirpath, qualimap_name),
            picard_dirpath=join(dirpath, picard_name),
            flagged_regions_dirpath=join(dirpath, flag_regions_name))
        self.bam = bam
        self.bed = bed
        self.vcf = vcf
        self.genome = genome

    def targetcov_done(self):
        if verify_file(self.targetcov_json_fpath) \
           and verify_file(self.targetcov_html_fpath):
           # and verify_file(self.targetcov_detailed_tsv):
            # info(self.targetcov_json_fpath + ', ' +
            #      self.targetcov_html_fpath + ', and ' +
            #      self.targetcov_detailed_tsv + ' exist.')
            return True
        return False

    # def ngscat_done(self):
    #     if verify_file(self.ngscat_html_fpath):
    #         # info(self.ngscat_html_fpath + ' exists.')
    #         return True
    #     return False
    #
    # def qualimap_done(self):
    #     if verify_file(self.qualimap_html_fpath):
    #         # info(self.qualimap_html_fpath + ' exists.')
    #         return True
    #     return False


class VarSample(BaseSample):
    def __init__(self, name, dirpath, **kwargs):
        BaseSample.__init__(self, name, dirpath, **kwargs)

        self.anno_vcf_fpath = None
        self.filt_vcf_fpath = None
        self.pass_filt_vcf_fpath = None
        self.filt_tsv_fpath = None
        self.varfilter_dirpath = None
        self.prep_vcf_fpath = None


# class TargQCStandaloneSample(BaseSample):
#     def __init__(self, name, dirpath, *args, **kwargs):
#         BaseSample.__init__(self, name, join(dirpath, '{sample}_{name}'), *args, **kwargs)
