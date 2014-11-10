from dircache import listdir
from genericpath import isdir, isfile
import os
import shutil
import sys
from collections import defaultdict, OrderedDict
from os.path import join, abspath, exists, pardir, splitext, basename, islink, dirname
import re
from source import logger
from source.logger import info, err, critical, warn
from source.calling_process import call
from source.config import load_yaml_config
from source.file_utils import verify_dir, verify_file, adjust_path, remove_quotes
from source.ngscat.bed_file import verify_bed, verify_bam
from source.targetcov.bam_file import index_bam
from source.tools_from_cnf import get_system_path
from source.file_utils import file_exists, safe_mkdir
from source.utils import OrderedDefaultDict


class Sample:
    def __init__(self, name, bcbio_structure=None, bam=None, bed=None, vcf=None):
        self.name = name
        self.bam = bam
        self.bed = bed
        self.vcf_by_callername = OrderedDict()  # string -> vcf_fpath
        self.vcf = vcf
        # self.filtered_vcf_by_callername = OrderedDict()
        # self.filtered_tsv_by_callername = OrderedDict()
        # self.filtered_maf_by_callername = OrderedDict()
        self.phenotype = None
        self.dirpath = None
        self.var_dirpath = None
        self.normal_match = None
        self.min_af = None

    # def get_filtered_vcf_fpath_by_callername(callername):
    #     fpath = join(
    #         self.bcbio_structure.final_dirpath,
    #         s.name,
    #         dirname,
    #         s.name + '-' + self.suf + ending)
    #     return verify_file(fpath)
    #
    # def find_filtered_vcf_by_callername(callername):
    #     return


    def find_vcf_fpath_by_callername(self, callername):
        return verify_file(self.get_vcf_fpath_by_callername(callername))

    def get_vcf_fpath_by_callername(self, callername):
        return join(self.dirpath, BCBioStructure.var_dir, self.name + '-' + callername + '.vcf')


    def find_anno_vcf_fpath_by_callername(self, callername):
        return verify_file(self.get_anno_vcf_fpath_by_callername(callername))

    def get_anno_vcf_fpath_by_callername(self, callername):
        return join(self.dirpath, BCBioStructure.varannotate_dir,
                    self.name + '-' + callername + BCBioStructure.anno_vcf_ending)


    def find_filt_vcf_by_callername(self, callername):
        return verify_file(self.get_filt_vcf_fpath_by_callername(callername))

    def find_filt_tsv_by_callername(self, callername):
        return verify_file(self.get_filt_tsv_fpath_by_callername(callername))

    def find_filt_maf_by_callername(self, callername):
        return verify_file(self.get_filt_maf_fpath_by_callername(callername))

    def find_pass_filt_vcf_by_callername(self, callername):
        return verify_file(self.get_pass_filt_vcf_fpath_by_callername(callername))

    def get_filt_vcf_fpath_by_callername(self, callername):
        return join(self.dirpath, BCBioStructure.varfilter_dir,
                    self.name + '-' + callername + BCBioStructure.filt_vcf_ending)

    def get_filt_tsv_fpath_by_callername(self, callername):
        return join(self.dirpath, BCBioStructure.varfilter_dir,
                    self.name + '-' + callername + BCBioStructure.filt_tsv_ending)

    def get_filt_maf_fpath_by_callername(self, callername):
        return join(self.dirpath, BCBioStructure.varfilter_dir,
                    self.name + '-' + callername + BCBioStructure.filt_maf_ending)

    def get_pass_filt_vcf_fpath_by_callername(self, callername):
        return join(self.dirpath, BCBioStructure.varfilter_dir,
                    self.name + '-' + callername + BCBioStructure.pass_filt_vcf_ending)

    def __str__(self):
        return self.name

    def for_json(self):
        return dict((k, v) for k, v in self.__dict__.items() if k != 'bcbio_structure')
    #     return dict(
    #         (k, (v if k != 'vcf_by_caller' else (dict((c.name, v) for c, v in v.items()))))
    #         for k, v in self.__dict__.items())

    def key_to_sort(self):
        m = re.search(r'\d+$', self.name)  # split_name_and_number
        if m:
            num = m.group()
            return self.name.split(num)[0], int(num)
        else:
            return self.name, 0

    @staticmethod
    def load(data, bcbio_structure):
        data['bcbio_structure'] = bcbio_structure
        sample = Sample(**data)
        sample.__dict__ = data
        return sample


class VariantCaller:
    def __init__(self, bcbio_structure, name):
        self.name = self.suf = name
        self.bcbio_structure = bcbio_structure
        self.samples = []

        self.summary_qc_report = None
        self.summary_qc_rep_fpaths = []

        self.combined_filt_maf_fpath = None

    def find_filt_maf_by_sample(self):
        return self._find_files_by_sample(BCBioStructure.varfilter_dir, BCBioStructure.filt_maf_ending)

    def find_fpaths_by_sample(self, dir_name, name, ext):
        return self._find_files_by_sample(dir_name, '.' + name + '.' + ext)

    def find_anno_vcf_by_sample(self):
        return self._find_files_by_sample(BCBioStructure.varannotate_dir, BCBioStructure.anno_vcf_ending)

    def get_filt_vcf_by_sample(self):
        return self._find_files_by_sample(BCBioStructure.varfilter_dir, BCBioStructure.filt_vcf_ending)

    def find_pass_filt_vcf_by_sample(self):
        return self._find_files_by_sample(BCBioStructure.varfilter_dir, BCBioStructure.pass_filt_vcf_ending)

    def _find_files_by_sample(self, dir_name, ending):
        files_by_sample = OrderedDict()

        for s in self.samples:
            fpath = join(
                self.bcbio_structure.final_dirpath,
                s.name,
                dir_name,
                s.name + '-' + self.suf + ending)

            if isfile(fpath):
                if verify_file(fpath):
                    files_by_sample[s.name] = fpath
            elif s.phenotype != 'normal':
                info('Warning: no ' + fpath + ' for ' + s.name + ', ' + self.name)

        return files_by_sample

    def __str__(self):
        return self.name

    def for_json(self):
        return {k: v for k, v in self.__dict__
                if k not in ['bcbio_structure', 'samples']}


class Batch:
    def __init__(self, name=None):
        self.name = name
        self.normal = None
        self.tumor = []

    def __str__(self):
        return self.name


class BCBioStructure:
    varfilter_name   = varfilter_dir                           = 'varFilter'
    varannotate_name = varannotate_dir                         = 'varAnnotate'
    targetseq_name   = targetseq_dir = targetseq_summary_dir   = 'targetSeq'
    cnv_dir                          = cnv_summary_dir         = 'cnv'
    varqc_name               = 'varQC'
    varqc_summary_name       = 'varQC_summary'
    varqc_after_name         = 'varQC_postVarFilter'
    varqc_after_summary_name = 'varQC_postVarFilter_summary'
    ngscat_name              = 'ngscat'
    qualimap_name            = 'qualimap'
    targqc_name              = 'targQC'
    fastqc_name              = 'fastqc'
    varqc_dir        = varqc_summary_dir       = join('qc', varqc_name)
    varqc_after_dir  = varqc_after_summary_dir = join('qc', varqc_after_name)
    ngscat_dir       = ngscat_summary_dir      = join('qc', ngscat_name)
    qualimap_dir     = qualimap_summary_dir    = join('qc', qualimap_name)
    targqc_summary_dir                         = join('qc', targqc_name)
    fastqc_dir       = fastqc_summary_dir      = join('qc', fastqc_name)
    seq2c_name       = 'Seq2C'
    combined_report_name = combined_report_dir = 'combined_report'
    detail_gene_report_baseending = '.details.gene'
    detail_sorted_gene_report_baseending = '.sorted.details.gene'
    detail_lowcov_gene_report_baseending = '.details.low_cov.gene'
    detail_highcov_gene_report_baseending = '.details.high_cov.gene'
    detail_gene_report_ending = detail_gene_report_baseending + '.txt'
    detail_gene_report_tsv_ending = detail_gene_report_baseending + '.tsv'
    seq2c_seq2cov_ending   = 'seq2c_seq2cov.txt'
    anno_vcf_ending  = '.anno.vcf'
    filt_vcf_ending  = '.anno.filt.vcf'
    pass_filt_vcf_ending = '.anno.filt.pass.vcf'
    filt_tsv_ending  = '.anno.filt.tsv'
    filt_maf_ending  = '.anno.filt.maf'
    var_dir          = 'var'

    def __init__(self, cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath=None, proc_name=None):
        self._set_final_dir(bcbio_cnf, bcbio_project_dirpath, final_dirpath)

        self.bcbio_cnf = cnf.bcbio_cnf
        self.cnf = cnf
        self.batches = OrderedDefaultDict(Batch)
        self.samples = []
        self.variant_callers = OrderedDict()

        # Date dirpath is from bcbio and named after fc_name, not our own project name
        self.date_dirpath = join(self.final_dirpath, bcbio_cnf['fc_date'] + '_' + bcbio_cnf['fc_name'])
        if not verify_dir(self.date_dirpath):
            err('Warning: no project directory of format {fc_date}_{fc_name}, creating ' + self.date_dirpath)
        safe_mkdir(self.date_dirpath)

        bcbio_project_dirname = basename(dirname(self.final_dirpath))
        bcbio_project_parent_dirname = basename(dirname(dirname(self.final_dirpath)))
        self.project_name = cnf.project_name or bcbio_project_parent_dirname + '_' + bcbio_project_dirname
        info('Project name: ' + self.project_name)
        self.cnf.name = proc_name or self.project_name

        self.set_up_log(cnf, proc_name, self.project_name, self.final_dirpath)

        self.work_dir = join(self.final_dirpath, pardir, 'work', 'post_processing')
        self.cnf.work_dir = self.work_dir
        if not isdir(self.work_dir):
            safe_mkdir(self.work_dir)

        self.var_dirpath = join(self.date_dirpath, BCBioStructure.var_dir)

        # Moving raw variants in the date dir to var_raw
        raw_var_dirpath = join(self.date_dirpath, BCBioStructure.var_dir + '_raw')
        for fname in os.listdir(self.date_dirpath):
            if '.vcf' in fname:
                if not isdir(raw_var_dirpath):
                    safe_mkdir(raw_var_dirpath)
                src_fpath = join(self.date_dirpath, fname)
                dst_fpath = join(raw_var_dirpath, fname)
                info('Moving ' + src_fpath + ' to ' + raw_var_dirpath)
                try:
                    os.rename(src_fpath, dst_fpath)
                except OSError:
                    pass

        info(' '.join(sys.argv))
        info()
        info('-' * 70)

        for sample in [self._read_sample_details(sample_info) for sample_info in bcbio_cnf['details']]:
            if sample.dirpath is None:
                err('For sample ' + sample.name + ', directory does not exist. Thus, skipping that sample.')
            else:
                self.samples.append(sample)

        if not self.samples:
            critical('No directory for any sample. Exiting.')

        for b in self.batches.values():
            for t_sample in b.tumor:
                t_sample.normal_match = b.normal

        self.samples.sort(key=lambda s: s.key_to_sort())
        for caller in self.variant_callers.values():
            caller.samples.sort(key=lambda s: s.key_to_sort())

        if not self.cnf.verbose:
            info('', ending='')

        # all_variantcallers = set()
        # for s_info in self.bcbio_cnf.details:
        #     all_variantcallers |= set(s_info['algorithm'].get('variantcaller')) or set()

        # samples_fpath = abspath(join(self.cnf.work_dir, 'samples.txt'))
        # with open(samples_fpath, 'w') as f:
        #     for sample_info in self.bcbio_cnf.details:
        #         sample = sample_info['description']
        #         f.write(sample + '\n')

        if not self.cnf.verbose:
            print ''
        else:
            info('Done loading BCBio structure.')

    def set_up_log(self, cnf, proc_name, project_name, project_fpath):
        logger.proc_name = proc_name
        logger.project_name = project_name
        logger.project_fpath = project_fpath
        logger.address = remove_quotes(cnf.email) if cnf.email else ''
        logger.smtp_host = cnf.smtp_host

        self.log_dirpath = join(self.date_dirpath, 'log')
        info('log_dirpath: ' + self.log_dirpath)
        safe_mkdir(self.log_dirpath)

        if not proc_name:
            info('self.cnf.name: ' + self.cnf.name)
            self.cnf.log = join(self.log_dirpath, self.cnf.name + '.log')
            info('log_dirpath: ' + self.cnf.log)
            i = 1
            if file_exists(self.cnf.log):
                bak_fpath = self.cnf.log + '.' + str(i)
                while isfile(bak_fpath):
                    bak_fpath = self.cnf.log + '.' + str(i)
                    i += 1
                os.rename(self.cnf.log, bak_fpath)
            elif isfile(self.cnf.log):
                try:
                    os.remove(self.cnf.log)
                except OSError:
                    pass

            logger.log_fpath = self.cnf.log

    @staticmethod
    def move_vcfs_to_var(sample):
        fpaths = []
        for fname in os.listdir(sample.dirpath):
            if any(fname.endswith(ending) for ending in
                   [BCBioStructure.filt_maf_ending,
                    BCBioStructure.filt_tsv_ending,
                    BCBioStructure.filt_vcf_ending,
                    BCBioStructure.filt_vcf_ending + '.idx']):
                continue

            if 'vcf' in fname.split('.') and \
                    not (islink(fname) and '.anno.filt' in fname):
                fpaths.append([sample, fname])

        if fpaths:
            if not exists(sample.var_dirpath):
                info('Creating "var" directory ' + sample.var_dirpath)
                safe_mkdir(sample.var_dirpath)

        for sample, fname in fpaths:
            src_fpath = join(sample.dirpath, fname)
            dst_fpath = join(sample.var_dirpath, fname)
            if exists(dst_fpath):
                try:
                    os.remove(dst_fpath)
                except OSError:
                    info('Cannot move ' + src_fpath + ' to ' + dst_fpath + ': dst exists, and permissions denied to remove it.')
                    continue
            safe_mkdir(sample.var_dirpath)
            info('Moving ' + src_fpath + ' to ' + dst_fpath)
            os.rename(src_fpath, dst_fpath)

    def _read_sample_details(self, sample_info):
        sample = Sample(name=sample_info['description'], bcbio_structure=self)

        info('Sample "' + sample.name + '"')
        if not self.cnf.verbose: info(ending='')

        sample.dirpath = adjust_path(join(self.final_dirpath, sample.name))
        if not verify_dir(sample.dirpath):
            sample.dirpath = None
            return sample

        self._set_bed_file(sample, sample_info)

        self._set_bam_file(sample)

        if 'min_allele_fraction' in sample_info['algorithm']:
            sample.min_af = float(sample_info['algorithm']['min_allele_fraction']) / 100

        sample.phenotype = None

        if 'metadata' in sample_info:
            sample.phenotype = sample_info['metadata'].get('phenotype') or 'tumor'
            info('Phenotype: ' + str(sample.phenotype))

            batch_names = sample_info['metadata']['batch']
            if isinstance(batch_names, basestring):
                batch_names = [batch_names]

            for batch_name in batch_names:
                if sample.phenotype == 'normal':
                    if self.batches[batch_name].normal:
                        critical('Multiple normal samples for batch ' + batch_name)
                    self.batches[batch_name].normal = sample

                elif sample.phenotype == 'tumor':
                    self.batches[batch_name].tumor.append(sample)

        sample.var_dirpath = adjust_path(join(sample.dirpath, BCBioStructure.var_dir))
        # self.move_vcfs_to_var(sample)

        to_exit = False
        variantcallers = sample_info['algorithm'].get('variantcaller') or []
        if isinstance(variantcallers, basestring):
            variantcallers = [variantcallers]

        for caller_name in variantcallers:
            info(caller_name)
            caller = self.variant_callers.get(caller_name)
            if not caller:
                self.variant_callers[caller_name] = VariantCaller(self, caller_name)

            to_exit, vcf_fpath = self._set_vcf_file(caller_name, sample, to_exit)
            self.variant_callers[caller_name].samples.append(sample)
            sample.vcf_by_callername[caller_name] = vcf_fpath

        if to_exit:
            sys.exit(1)

        info()
        return sample

    def _set_final_dir(self, bcbio_cnf, bcbio_project_dirpath, final_dirpath=None):
        if final_dirpath:
            self.final_dirpath = final_dirpath
        elif 'upload' in bcbio_cnf and 'dir' in bcbio_cnf['upload']:
            final_dirname = bcbio_cnf['upload']['dir']
            self.final_dirpath = adjust_path(join(bcbio_project_dirpath, 'config', final_dirname))
            if not verify_dir(self.final_dirpath, 'upload directory specified in the bcbio config'):
                sys.exit(1)
        else:
            self.final_dirpath = join(bcbio_project_dirpath, 'final')
            if not verify_dir(self.final_dirpath):
                critical('If final directory it is not named "final", please, specify it in the bcbio config.')
        info('Final dirpath: ' + self.final_dirpath)

    def _set_bed_file(self, sample, sample_info):
        bed = None
        if self.cnf.bed:  # Custom BED provided in command line?
            bed = adjust_path(self.cnf.bed)
            if not verify_bed(bed):
                sys.exit(1)
        elif sample_info['algorithm'].get('variant_regions'):  # Variant regions?
            bed = adjust_path(sample_info['algorithm']['variant_regions'])
            if not verify_bed(bed):
                sys.exit(1)
        elif self.cnf.genome.exons:
            warn('Warning: no amplicon BED file provided, using exons instead.')
            bed = self.cnf.genome.exons
            if not verify_bed(bed):
                sys.exit(1)
        else:
            err('No BED file for sample, no default BED file and exons (or cannot read them)'
                ' - skipping targetSeq reproting.')
        sample.bed = bed
        if sample.bed:
            info('BED file for ' + sample.name + ': ' + sample.bed)
        else:
            err('No BED file for ' + sample.name)

    def _set_bam_file(self, sample):
        bam = adjust_path(join(sample.dirpath, sample.name + '.bam'))
        if verify_bam(bam):
            sample.bam = bam
            info('BAM file for ' + sample.name + ': ' + sample.bam)
            index_bam(self.cnf, bam)
        else:
            sample.bam = None
            err('No BAM file for ' + sample.name)

    def _set_vcf_file(self, caller_name, sample, to_exit):
        vcf_fname = sample.name + '-' + caller_name + '.vcf'
        vcf_fpath = adjust_path(join(sample.var_dirpath, vcf_fname))  # in var
        if not isfile(vcf_fpath):  # not in var, looking in sample dir
            vcf_fpath = adjust_path(join(sample.dirpath, vcf_fname))  # in sample dir
        _ungzip_if_needed(self.cnf, vcf_fpath)
        if isfile(vcf_fpath) and not verify_file(vcf_fpath):  # bad file, error :(
            err('Error: Phenotype is ' + str(sample.phenotype) + ', and VCF file is empty.')
            to_exit = True
            vcf_fpath = None
        if not isfile(vcf_fpath):
            if sample.phenotype != 'normal':  # no VCF file is OK if phenotype is normal, otherwise - warning
                err('Warning: Phenotype is ' + str(sample.phenotype) + ', and no VCF file.')
            vcf_fpath = None
        if vcf_fpath:
            info(vcf_fpath)
        return to_exit, vcf_fpath

    def find_gene_reports_by_sample(self):
        return dict((sname, verify_file(fpath))
                    for sname, fpath in self.get_gene_report_fpaths_by_sample().items())

    def get_gene_report_fpaths_by_sample(self):
        return self._get_fpaths_per_sample(
            BCBioStructure.targetseq_dir,
            lambda sample: sample.name + '.' +
                           BCBioStructure.targetseq_dir +
                           BCBioStructure.detail_gene_report_ending)

    def get_gene_report_tsv_fpaths_by_sample(self):
        return self._get_fpaths_per_sample(
            BCBioStructure.targetseq_dir,
            lambda sample: sample.name + '.' +
                           BCBioStructure.targetseq_dir +
                           BCBioStructure.detail_gene_report_tsv_ending)

    def find_targetcov_reports_by_sample(self, ext='json'):
        return dict((sname, verify_file(fpath))
                    for sname, fpath in self.get_targetcov_report_fpaths_by_sample(ext).items())

    def get_targetcov_report_fpaths_by_sample(self, ext='json'):
        return self._get_fpaths_per_sample(
            BCBioStructure.targetseq_dir,
            lambda sample: sample.name + '.' +
                           BCBioStructure.targetseq_dir + '.' + ext)

    def find_ngscat_reports_by_sample(self):
        return dict((sname, verify_file(fpath))
                    for sname, fpath in self.get_ngscat_report_fpaths_by_sample().items())

    def get_ngscat_report_fpaths_by_sample(self):
        return self._get_fpaths_per_sample(
            BCBioStructure.ngscat_dir,
            lambda sample: 'captureQC.html')

    def find_qualimap_reports_by_sample(self):
        return dict((sname, verify_file(fpath))
                    for sname, fpath in self.get_qualimap_report_fpaths_by_sample().items())

    def get_qualimap_report_fpaths_by_sample(self):
        return self._get_fpaths_per_sample(
            BCBioStructure.qualimap_dir,
            lambda sample: 'qualimapReport.html')

    def get_fastqc_report_fpaths_by_sample(self):
        return self._get_fpaths_per_sample(
            BCBioStructure.fastqc_dir,
            lambda sample: 'fastqc_report.html')

    def _get_fpaths_per_sample(self, base_dir, get_name_fn):
        fpaths_by_sample = OrderedDict()

        for sample in self.samples:
            report_fpath = join(self.final_dirpath, sample.name, base_dir, get_name_fn(sample))
            info(report_fpath)

            if verify_file(report_fpath):
                fpaths_by_sample[sample.name] = report_fpath

        # if len(fpaths_by_sample) < len(self.samples):
        #     raise RuntimeError('No ')

        return fpaths_by_sample

    def _find_files_per_sample(self, base_dir, get_name_fn):
        return dict((sname, verify_file(fpath))
                    for sname, fpath in self._get_fpaths_per_sample(base_dir, get_name_fn).items())

    def clean(self):
        for sample in self.samples:
            info('Sample ' + sample.name)

            for dic in [sample.filtered_vcf_by_callername,
                        sample.filtered_tsv_by_callername,
                        sample.filtered_maf_by_callername]:
                for c, fpath in dic.items():
                    try:
                        os.unlink(fpath)
                        info('Removed symlink ' + fpath)
                    except OSError:
                        pass

            for fname in listdir(sample.var_dirpath):
                if not fname.startswith('.'):
                    fpath = join(sample.var_dirpath, fname)
                    os.rename(fpath, join(sample.dirpath, fname))

            for dirname in [BCBioStructure.varannotate_dir,
                            BCBioStructure.varfilter_dir,
                            BCBioStructure.varqc_dir,
                            BCBioStructure.varqc_after_dir,
                            BCBioStructure.ngscat_dir,
                            BCBioStructure.qualimap_dir,
                            BCBioStructure.targetseq_dir,
                            BCBioStructure.var_dir]:
                dirpath = join(sample.dirpath, dirname)
                if isdir(dirpath):
                    info('  removing ' + dirpath)
                    shutil.rmtree(dirpath)
            info()

        for dirname in [BCBioStructure.targetseq_summary_dir,
                        BCBioStructure.cnv_summary_dir,
                        BCBioStructure.varqc_summary_dir,
                        BCBioStructure.varqc_after_summary_dir,
                        BCBioStructure.ngscat_summary_dir,
                        BCBioStructure.qualimap_summary_dir]:
            dirpath = join(self.date_dirpath, dirname)
            if isdir(dirpath):
                info('  removing ' + dirpath)
                shutil.rmtree(dirpath)


def load_bcbio_cnf(cnf, config_dirpath):
    yaml_files = [join(config_dirpath, fname)
                  for fname in listdir(config_dirpath)
                  if fname.endswith('.yaml')]

    if len(yaml_files) == 0:
        critical('No YAML file in the config directory.')

    config_fpaths = [fpath for fpath in yaml_files if not any(n in fpath for n in ['run_info', 'system_info'])]
    if not config_fpaths:
        critical('No BCBio YAMLs in the config directory (only ' + ', '.join(map(basename, yaml_files)) + ')')

    yaml_fpath = config_fpaths[0]
    if len(config_fpaths) > 1:
        some_yaml_files = [f for f in config_fpaths if splitext(basename(f))[0] in config_dirpath]
        if len(some_yaml_files) == 0:
            critical('More than one YAML file in the config directory ' + ' '.join(config_fpaths) +
                     ', and no YAML file named after the project.')
        yaml_fpath = some_yaml_files[0]

    yaml_file = abspath(yaml_fpath)

    info('Using bcbio YAML config ' + yaml_file)

    return load_yaml_config(yaml_file)


def _normalize(name):
    return name.lower().replace('_', '').replace('-', '')


def _ungzip_if_needed(cnf, fpath):
    if fpath.endswith('.gz'):
        fpath = fpath[:-3]
    if not file_exists(fpath) and file_exists(fpath + '.gz'):
        gz_fpath = fpath + '.gz'
        gunzip = get_system_path(cnf, 'gunzip')
        cmdline = '{gunzip} -c {gz_fpath}'.format(**locals())
        res = call(cnf, cmdline, output_fpath=fpath, exit_on_error=False)
        info()
        if not res:
            return None
    return fpath


# def get_trailing_number(string):
#     m = re.search(r'\d+$', string)
#     return int(m.group()) if m else None

