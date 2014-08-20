from dircache import listdir
from genericpath import isdir, isfile
import os
import shutil
import sys
from collections import defaultdict, OrderedDict
from os.path import join, abspath, exists, pardir, splitext, basename, islink
import re
from source import logger
from source.logger import info, err, critical
from source.calling_process import call
from source.config import load_yaml_config
from source.file_utils import verify_dir, verify_file, adjust_path
from source.ngscat.bed_file import verify_bed, verify_bam
from source.tools_from_cnf import get_tool_cmdline
from source.file_utils import file_exists, safe_mkdir
from source.utils import OrderedDefaultDict


class Sample:
    def __init__(self, name, bam=None, bed=None, vcf=None):
        self.name = name
        self.bam = bam
        self.bed = bed
        self.vcf_by_callername = OrderedDict()  # string -> vcf_fpath
        self.filtered_vcf_by_callername = OrderedDict()
        self.filtered_tsv_by_callername = OrderedDict()
        self.filtered_maf_by_callername = OrderedDict()
        self.phenotype = None
        self.dirpath = None
        self.var_dirpath = None
        self.normal_match = None

    def __str__(self):
        return self.name

    def for_json(self):
        return self.__dict__
    #     return dict(
    #         (k, (v if k != 'vcf_by_caller' else (dict((c.name, v) for c, v in v.items()))))
    #         for k, v in self.__dict__.items())


class VariantCaller:
    def __init__(self, bcbio_structure, name):
        self.name = self.suf = name
        self.summary_qc_report = None
        self.summary_qc_rep_fpaths = []
        self.anno_vcf_fpaths = OrderedDict()
        self.anno_filt_vcf_fpaths = OrderedDict()

        self.bcbio_structure = bcbio_structure
        self.samples = []

    def _get_files_by_sample(self, dirname, ending):
        files_by_sample = OrderedDict()

        to_exit = False
        for s in self.samples:
            if self.name in s.vcf_by_callername:
                fpath = join(
                    self.bcbio_structure.final_dirpath,
                    s.name,
                    dirname,
                    s.name + '-' + self.suf + ending)

                if verify_file(fpath):
                    files_by_sample[s] = fpath
                else:
                    to_exit = True

        if to_exit:
            sys.exit(1)

        return files_by_sample

    def get_qc_reports_by_samples(self):
        return self._get_files_by_sample(
            BCBioStructure.varqc_dir, '.' + BCBioStructure.varqc_name + '.json')

    def get_anno_vcf_by_samples(self):
        return self._get_files_by_sample(
            BCBioStructure.varannotate_dir, BCBioStructure.anno_vcf_ending)

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
    varfilter_name      = varfilter_dir                           = 'varFilter'
    varannotate_name    = varannotate_dir                         = 'varAnnotate'
    targetseq_name      = targetseq_dir = targetseq_summary_dir   = 'targetSeq'
    cnv_dir                             = cnv_summary_dir         = 'cnv'
    varqc_name                          = varqc_summary_dir       = 'varQC'
    varqc_after_name                    = varqc_after_summary_dir = 'varQC_postVarFilter'
    ngscat_name                         = ngscat_summary_dir      = 'ngscat'
    qualimap_name                       = qualimap_summary_dir    = 'qualimap'
    varqc_dir           = join('qc', varqc_name)
    varqc_after_dir     = join('qc', varqc_after_name)
    ngscat_dir          = join('qc', ngscat_name)
    qualimap_dir        = join('qc', qualimap_name)
    seq2c_name          = 'Seq2C'
    detail_gene_report_ending = '.details.gene.txt'
    anno_vcf_ending     = '.anno.vcf'
    filt_vcf_ending     = '.anno.filt.vcf'
    filt_tsv_ending     = '.anno.filt.tsv'
    filt_maf_ending     = '.anno.filt.passed.maf'
    var_dir             = 'var'


    def __init__(self, cnf, bcbio_final_dirpath, bcbio_cnf, proc_name=None):
        self.final_dirpath = bcbio_final_dirpath
        self.bcbio_cnf = bcbio_cnf
        self.cnf = cnf
        self.batches = OrderedDefaultDict(Batch)
        self.samples = []
        self.variant_callers = OrderedDict()
        self.project_name = bcbio_cnf.fc_name
        self.cnf.name = proc_name or self.project_name or critical('No fc_name in bcbio YAML file.')

        self.date_dirpath = join(bcbio_final_dirpath, bcbio_cnf.fc_date + '_' + bcbio_cnf.fc_name)
        if not verify_dir(self.date_dirpath): err('Warning: no project directory of format {fc_date}_{fc_name}, creating ' + self.date_dirpath)
        safe_mkdir(self.date_dirpath)

        self.set_up_log(proc_name)

        self.work_dir = join(cnf.bcbio_final_dir, pardir, 'work', 'post_processing')
        self.cnf.work_dir = self.work_dir
        if not isdir(self.work_dir):
            safe_mkdir(self.work_dir)

        info(' '.join(sys.argv))
        info()
        info('-' * 70)

        self.samples = [self._read_sample_details(sample_info) for sample_info in self.bcbio_cnf.details]
        if any(s is None for s in self.samples):
            sys.exit(1)

        for b in self.batches.values():
            for t_sample in b.tumor:
                t_sample.normal_match = b.normal

        if all(get_trailing_number(s.name) for s in self.samples):
            self.samples.sort(key=lambda s: split_name_and_number(s.name))

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
            info('Done.')

    def set_up_log(self, proc_name):
        self.log_dirpath = join(self.date_dirpath, 'log')
        safe_mkdir(self.log_dirpath)

        if not proc_name:
            self.cnf.log = join(self.log_dirpath, self.cnf.name + '.log')
            if file_exists(self.cnf.log):
                bak_fpath = self.cnf.log + '.bak'
                if isfile(bak_fpath):
                    os.remove(bak_fpath)
                os.rename(self.cnf.log, bak_fpath)
            elif isfile(self.cnf.log):
                try:
                    os.remove(self.cnf.log)
                except OSError:
                    pass

            logger.log_fpath = self.cnf.log

    @staticmethod
    def _ungzip_if_needed(cnf, fpath):
        if not file_exists(fpath) and file_exists(fpath + '.gz'):
            gz_fpath = fpath + '.gz'
            gunzip = get_tool_cmdline(cnf, 'gunzip')
            cmdline = '{gunzip} -c {gz_fpath}'.format(**locals())
            call(cnf, cmdline, output_fpath=fpath)
            info()

    @staticmethod
    def _move_vcfs_to_var(sample):
        if not exists(sample.var_dirpath):
            info('Creating "var" directory ' + sample.var_dirpath)
            safe_mkdir(sample.var_dirpath)

        for fname in os.listdir(sample.dirpath):
            if any(fname.endswith(ending) for ending in
                   [BCBioStructure.filt_maf_ending,
                    BCBioStructure.filt_vcf_ending,
                    BCBioStructure.filt_tsv_ending]):
                continue

            if 'vcf' in fname.split('.') and \
                    not (islink(fname) and fname.endswith('.anno.filt.vcf')):
                src_fpath = join(sample.dirpath, fname)
                dst_fpath = join(sample.var_dirpath, fname)
                if exists(dst_fpath):
                    os.remove(dst_fpath)
                safe_mkdir(sample.var_dirpath)
                info('Moving ' + src_fpath + ' to ' + dst_fpath)
                os.rename(src_fpath, dst_fpath)

    def _read_sample_details(self, sample_info):
        sample = Sample(name=sample_info['description'])
        self.samples.append(sample)

        info('Sample "' + sample.name + '"')
        if not self.cnf.verbose: info(ending='')

        sample.dirpath = adjust_path(join(self.final_dirpath, sample.name))
        if not verify_dir(sample.dirpath): sys.exit(1)

        bed = adjust_path(sample_info['algorithm'].get('variant_regions'))
        bam = adjust_path(join(sample.dirpath, sample.name + '-ready.bam'))
        sample.bed = bed if verify_bed(bed) else None
        sample.bam = bam if verify_bam(bam) else None

        safe_mkdir(join(sample.dirpath, 'qc'))

        sample.phenotype = None

        if 'metadata' in sample_info:
            sample.phenotype = sample_info['metadata']['phenotype']
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
        self._move_vcfs_to_var(sample)

        to_exit = False
        for caller_name in sample_info['algorithm'].get('variantcaller') or []:
            vcf_fname = sample.name + '-' + caller_name + '.vcf'
            vcf_fpath = adjust_path(join(sample.var_dirpath, vcf_fname))
            self._ungzip_if_needed(self.cnf, vcf_fpath)

            if isfile(vcf_fpath) and not verify_file(vcf_fpath):
                to_exit = True
                continue

            if not file_exists(vcf_fpath):
                if sample.phenotype != 'normal':
                    err('Phenotype is ' + str(sample.phenotype) + ', and VCF does not exist.')
                vcf_fpath = None

            if vcf_fpath:
                if caller_name not in self.variant_callers:
                    self.variant_callers[caller_name] = VariantCaller(self, caller_name)
                self.variant_callers[caller_name].samples.append(sample)
                info(vcf_fpath)

            sample.vcf_by_callername[caller_name] = vcf_fpath  # could be None, that's OK

            # And filtered symlinks
            for ending, dic in zip([BCBioStructure.filt_vcf_ending,
                                    BCBioStructure.filt_tsv_ending,
                                    BCBioStructure.filt_maf_ending],
                                   [sample.filtered_vcf_by_callername,
                                    sample.filtered_tsv_by_callername,
                                    sample.filtered_maf_by_callername]):

                fpath = join(sample.dirpath, sample.name + '-' + caller_name + ending)
                if islink(fpath) or (isfile(fpath) and verify_file(fpath)):
                    dic[caller_name] = fpath

        if to_exit:
            sys.exit(1)

        info()

        return sample

    def get_gene_reports_by_sample(self):
        return self.get_per_sample_fpaths_for_bcbio_final_dir(
            BCBioStructure.targetseq_dir,
            lambda sample: sample.name + '.' +
                           BCBioStructure.targetseq_dir +
                           BCBioStructure.detail_gene_report_ending)

    def get_targetcov_json_by_sample(self):
        return self.get_per_sample_fpaths_for_bcbio_final_dir(
            BCBioStructure.targetseq_dir,
            lambda sample: sample.name + '.' + BCBioStructure.targetseq_dir + '.json')

    def get_ngscat_html_by_sample(self):
        return self.get_per_sample_fpaths_for_bcbio_final_dir(
            BCBioStructure.ngscat_dir,
            lambda sample: 'captureQC.html')

    def get_qualimap_html_by_sample(self):
        return self.get_per_sample_fpaths_for_bcbio_final_dir(
            BCBioStructure.qualimap_dir,
            lambda sample: 'qualimapReport.html')

    def get_per_sample_fpaths_for_bcbio_final_dir(self, base_dir, get_name_fn):
        fpaths = OrderedDict()

        for sample in self.samples:
            report_fpath = join(self.final_dirpath, sample.name, base_dir, get_name_fn(sample))
            info(report_fpath)

            if verify_file(report_fpath):
                fpaths[sample] = report_fpath

        if len(fpaths) < len(self.samples):
            sys.exit(1)

        return fpaths

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


def load_bcbio_cnf(cnf):
    bcbio_config_dirpath = join(cnf.bcbio_final_dir, pardir, 'config')
    yaml_files = [join(bcbio_config_dirpath, fname)
                  for fname in listdir(bcbio_config_dirpath)
                  if fname.endswith('.yaml')]

    if len(yaml_files) == 0:
        critical('No YAML file in the config directory.')

    config_fpaths = [fpath for fpath in yaml_files if not any(n in fpath for n in ['run_info', 'system_info'])]
    if not config_fpaths:
        critical('No BCBio YAMLs in the config directory (only ' + ', '.join(map(basename, yaml_files)) + ')')

    yaml_fpath = config_fpaths[0]
    if len(config_fpaths) > 1:
        some_yaml_files = [f for f in config_fpaths if splitext(basename(f))[0] in cnf.bcbio_final_dir]
        if len(some_yaml_files) == 0:
            critical('More than one YAML file in the config directory ' + ' '.join(config_fpaths) +
                     ', and no YAML file named after the project.')
        yaml_fpath = some_yaml_files[0]

    yaml_file = abspath(yaml_fpath)

    info('Using bcbio YAML config ' + yaml_file)

    cnf.bcbio_cnf = load_yaml_config(yaml_file)


def _normalize(name):
    return name.lower().replace('_', '').replace('-', '')


def split_name_and_number(string):
    m = re.search(r'\d+$', string)
    if m:
        num = m.group()
        return string.split(num)[0], int(num)
    else:
        return string, None


def get_trailing_number(string):
    m = re.search(r'\d+$', string)
    return int(m.group()) if m else None
