from collections import OrderedDict
from genericpath import isdir, exists
import os
import shutil
import pickle
import operator

from os.path import basename, join, isfile, dirname, splitext, islink, pardir, abspath
from joblib import Parallel, delayed
##from memory_profiler import profile
import sys

from source.bcbio_structure import BCBioStructure
from source.calling_process import call
from source.config import defaults
from source.tools_from_cnf import get_script_cmdline, get_system_path
from source.variants.Effect import Effect
from source.logger import step_greetings, info, critical, err, warn
from source.variants.anno import _snpsift_annotate
from source.variants.vcf_processing import iterate_vcf, vcf_one_per_line, \
    get_main_sample_index, tabix_vcf, vcf_merge, leave_main_sample
from source.utils import mean
from source.file_utils import safe_mkdir, add_suffix, verify_file
from source.variants.tsv import make_tsv
from source.variants.vcf_processing import remove_rejected, convert_to_maf, vcf_is_empty, igvtools_index
from source.logger import info


# TODO:
# add filters to VCF header

class Filter:
    filt_cnf = None
    main_sample = None

    def __init__(self, word, check, required=True):
        self.check = check
        self.required = required
        self.word = word.upper()

        # self.num_passed = 0
        # self.num_rejected = 0

    @staticmethod
    def __remove_pass(rec):
        if rec.FILTER == ['PASS']:
            rec.FILTER = None

    # Call check function. If False, add value to the FILTER field (unless only_check)
    def apply(self, rec, only_check=False, *args, **kvargs):
        # if rec.is_rejected():
        #     return False

        if only_check:
            return self.check(rec, *args, **kvargs)

        if self.check(rec, *args, **kvargs):  # True if passed, False otherwise
            # self.num_passed += 1
            return True
        else:
            # self.num_rejected += 1

            if rec.FILTER:
                self.__remove_pass(rec)
            rec.add_filter(self.word)
            return False


class CnfFilter(Filter):  # get value from config by key, compare with provided value
    def __init__(self, key, check, cnf=None, *args, **kvargs):
        self.cnf = cnf or Filter.filt_cnf
        self.key = key

        self.cnf_val = self.cnf.get(key)
        if self.cnf_val:
            try:
                cnf_val = float(self.cnf_val)
            except ValueError:
                try:
                    cnf_val = int(self.cnf_val)
                except ValueError:
                    pass
                else:
                    self.cnf_val = cnf_val
            else:
                self.cnf_val = cnf_val

        Filter.__init__(self, key.upper(), check, *args, **kvargs)

    def apply(self, rec, only_check=False, *args, **kvargs):
        if self.cnf_val is None:
            return True

        return Filter.apply(self, rec, only_check, *args, **kvargs)


class InfoFilter(CnfFilter):  # cnf filter, but compare with value from INFO by key, using the op
    def __init__(self, cnf_key, info_key, op=operator.gt, *args, **kwargs):
        # True if PASS
        def check(rec):
            anno_val = rec.get_val(info_key)
            if anno_val is None:
                if self.required:
                    critical(
                        'Error: no field ' + info_key + ' in INFO or SAMPLE for variant ' +
                        rec.get_variant() + ' - required to test ' + cnf_key)
                else:
                    return True  # PASS

            return op(anno_val, self.cnf_val)

        CnfFilter.__init__(self, cnf_key, check, *args, **kwargs)


class EffectFilter(CnfFilter):  # cnf filter, but compare with value from EFF
    def __init__(self, cnf_key, *args, **kwargs):
        def check(rec):
            if 'EFF' not in rec.INFO:
                # err('Warning: EFF field is missing for variant ' + str(rec.get_variant()))
                return False
            else:
                return any(eff.impact.upper() in self.important_impacts
                       for eff in map(Effect, rec.INFO['EFF']))

        CnfFilter.__init__(self, cnf_key, check, *args, **kwargs)

        self.important_impacts = ''
        if self.cnf_val:
            if not isinstance(self.cnf_val, basestring):
                critical('The value of filter impact should be string, like MODERATE|HIGH')
                return
            self.important_impacts = [s.upper() for s in self.cnf_val.split('|')]


class VariantInfo:  # collects information (here, AFs) for all (chrom, pos, ref, alt) tuples
    def __init__(self, var_str):  # var_str="CHROM:POS:REF:ALT"
        self.var_str = var_str
        self.afs = []
        self._avg_af = None

    def occurence_num(self):
        return len(self.afs)

    def frac(self):
        return float(self.occurence_num()) / len(filtering.sample_names)

    def avg_af(self):
        if self._avg_af is None:
            self._avg_af = mean(self.afs)
        return self._avg_af


def process_vcf(vcf_fpath, fun, suffix, cnf=None, *args, **kwargs):
    return iterate_vcf(cnf, vcf_fpath, fun, suffix, cnf_=cnf, self_=filtering, *args, **kwargs)


def rm_prev_round(vcf_fpath, sample_name):
    return process_vcf(vcf_fpath, proc_line_remove_prev_filter, 'rm_prev', cnfs_for_sample_names[sample_name])


def first_round(vcf_fpath, sample_name):
    variant_dict = dict()  # tuples of (chrom, pos, ref, alt)
    control_variants = set()  # variants from controls samples
    res = process_vcf(
        vcf_fpath, proc_line_1st_round, 'r1', cnfs_for_sample_names[sample_name],
        variant_dict=variant_dict, control_vars=control_variants)

    return res, variant_dict, control_variants


def second_round(vcf_fpath, sample_name):
    cnf = cnfs_for_sample_names[sample_name]
    main_sample_index = get_main_sample_index(vcf_fpath, sample_name)
    cnf.main_sample_index = main_sample_index

    #TODO: remove this since it is made in annotation
    res = _snpsift_annotate(cnf, cnf.get('clinvar'), 'clinvar', vcf_fpath)
    if not res:
        err('Could not annotate with clinvar.')
    else:
        vcf_fpath = res

    return process_vcf(vcf_fpath, proc_line_2nd_round, 'r2', cnf)


def impact_round(vcf_fpath, sample_name):
    return process_vcf(vcf_fpath, proc_line_impact, 'impact', cnfs_for_sample_names[sample_name])


def one_per_line(vcf_fpath, sample_name):
    return vcf_one_per_line(cnfs_for_sample_names[sample_name], vcf_fpath)


class Filtering:
    cnf = None
    filt_cnf = None

    def __init__(self, cnf, bcbio_structure, caller):
        Filtering.cnf = cnf
        Filtering.filt_cnf = Filter.filt_cnf = filt_cnf = cnf.variant_filtering.__dict__

        self.caller = caller
        self.control_vars = set()
        self.sample_names = set([s.name for s in caller.samples])
        self.variant_dict = dict()  # "CHROM:POS:REF:ALT" -> VariantInfo("CHROM:POS:REF:ALT", afs)
        self.polymorphic_variants = None

        self.round1_filters = []
        if filt_cnf.get('filt_depth') is not None:
            self.round1_filters.append(InfoFilter('filt_depth', 'DP', required=False))
        if filt_cnf.get('filt_q_mean') is not None:
            self.round1_filters.append(InfoFilter('filt_q_mean', 'QUAL', required=False))
        if filt_cnf.get('filt_p_mean') is not None:
            self.round1_filters.append(InfoFilter('filt_p_mean', 'PMEAN', required=False))

        self.control = filt_cnf.get('control')

        self.impact_filter = EffectFilter('impact')

        self.polymorphic_filter = Filter('POLYMORPHIC', lambda rec: rec)

        self.min_freq_filters = dict()
        for sample in caller.samples:
            def min_freq_filter_check(rec):
                af = rec.get_af(cnf.main_sample_index, caller.name)
                if af and af < sample.min_af:
                    return False
                return True
            self.min_freq_filters[sample.name] = CnfFilter('min_freq', min_freq_filter_check,
                cnfs_for_sample_names[sample.name], required=True)

        self.round2_filters = [
            InfoFilter('min_p_mean', 'PMEAN', required=False),
            InfoFilter('min_q_mean', 'QMEAN', required=False),
            InfoFilter('mean_mq', 'MQ', required=False),
            InfoFilter('signal_noise', 'SN', required=False),
            InfoFilter('mean_vd', 'VD', required=False)]

        def dup_filter_check(rec, main_sample_index):
            pstd = rec.get_val('PSTD')
            bias = rec.get_bias(main_sample_index)

            # all variants from one position in reads
            if pstd is not None and bias is not None:
                return not (pstd == 0 and bias[-1] != '0' and bias[-1] != '1')
            return True
        self.dup_filter = Filter('DUP', dup_filter_check)

        def bias_filter_check(rec, cls, main_sample_index, caller_name):  # Filter novel variants with strand bias.
            if not Filter.filt_cnf['bias'] is True:
                return True  # I we don't need to check bias, just keep the variant

            if not cls in ['Novel', 'dbSNP']:
                return True  # Check only for novel and dnSNP

            if not (rec.get_bias(main_sample_index) and rec.get_bias(main_sample_index) in ['2:1', '2:0']):
                return True  # Variant has to have bias and with proper values

            return rec.get_af(main_sample_index, caller_name) >= 0.3
        self.bias_filter = CnfFilter('bias', bias_filter_check)

        def nonclnsnp_filter_check(rec, cls):
            if cls in ['COSMIC']:
                return True

            return rec.check_clnsig() != 'not_significant'
        self.nonclnsnp_filter = Filter('NonClnSNP', nonclnsnp_filter_check)

        def multi_filter_check(rec, var_n, frac, avg_af):
            if (rec.ID is None and      # reject if novel and present in [fraction] samples
                frac > Filter.filt_cnf['fraction'] and
                var_n >= Filter.filt_cnf['sample_cnt'] and
                avg_af < Filter.filt_cnf['freq']):
                return False
            else:
                return True
        self.multi_filter = Filter('MULTI', multi_filter_check)

        def max_rate_filter_check(rec, frac, af):  # reject if present in [max_ratio] samples
            if frac >= Filter.filt_cnf['max_ratio'] and af and af < 0.3:
                return False
            return True
        self.max_rate_filter = CnfFilter('max_ratio', max_rate_filter_check)

        # self.undet_sample_filter = Filter('UNDET_SAMPLE', lambda rec: False)

        def control_filter_check(rec):
            return not (filt_cnf['control'] and rec.get_variant() in self.control_vars)
        self.control_filter = CnfFilter('control', control_filter_check)

    def run_regengineered_filtering(self, vcf_fpaths, sample_names, n_threads=1):
        step_greetings('Filtering')

        global filtering
        filtering = self

        info('Fixing previous . values to PASS')
        vcf_fpaths = Parallel(n_jobs=n_threads)(
            delayed(rm_prev_round)(vcf_fpath, s)
            for vcf_fpath, s in zip(vcf_fpaths, sample_names))
        info()

        info('First round')
        results = Parallel(n_jobs=n_threads)(
            delayed(first_round)(vcf_fpath, s)
            for vcf_fpath, s in zip(vcf_fpaths, sample_names))
        info()

        control_vars_dump_fpath = join(self.cnf.work_dir, self.caller.name + '_control_vars.txt')
        variants_dump_fpath = join(self.cnf.work_dir, self.caller.name + 'variants.txt')

        if not [all([variant_dict, control_vars]) for (_, variant_dict, control_vars) in results]:
            info('Skipped first round, restoring varks and controls vars.')
            if not verify_file(control_vars_dump_fpath) or not verify_file(variants_dump_fpath):
                critical('Cannot restore varks and control_vars, please, run without the --reuse flag.')

            with open(control_vars_dump_fpath) as f:
                info('Loading control vars...')
                self.control_vars = pickle.load(f)
                info('Loaded control vars: ' + str(len(self.control_vars)))

            with open(variants_dump_fpath) as f:
                info('Loading varks...')
                self.variant_dict = pickle.load(f)
                info('Loaded varks: ' + str(len(self.variant_dict)))

        else:
            for vcf_fpath, variant_dict, control_vars in results:
                for var_str, variant_info in variant_dict.items():
                    if var_str not in self.variant_dict:
                        self.variant_dict[var_str] = variant_info
                    else:
                        self.variant_dict[var_str].afs.extend(variant_info.afs)

                self.control_vars.update(control_vars)

            with open(control_vars_dump_fpath, 'w') as f:
                info('Saving control vars...')
                pickle.dump(self.control_vars, f)
                info('Saved control vars: ' + str(len(self.control_vars)))

            with open(variants_dump_fpath, 'w') as f:
                info('Saving varks...')
                pickle.dump(self.variant_dict, f)
                info('Saved varks: ' + str(len(self.variant_dict)))

        vcf_fpaths = [vcf_fpath for vcf_fpath, _, _ in results]
        filtering = self
        info()

        info('One effect per line')
        vcf_fpaths = Parallel(n_jobs=n_threads)(delayed(one_per_line)(vcf_fpath, s)
                     for vcf_fpath, s in zip(vcf_fpaths, sample_names))
        info()

        info('Second round')
        vcf_fpaths = Parallel(n_jobs=n_threads)(delayed(second_round)(vcf_fpath, s)
                     for vcf_fpath, s in zip(vcf_fpaths, sample_names))
        info()

        # if 'polymorphic_variants' in self.filt_cnf:
        #     polymorphic_variants_path = self.filt_cnf['polymorphic_variants']
        #     if verify_file(polymorphic_variants_path):
        #         info('Filtering out polymorphic variants')
        #         with open(polymorphic_variants_path) as f:
        #

        return vcf_fpaths


def proc_line_remove_prev_filter(rec, cnf_, self_):
    if not rec.FILTER:
        rec.FILTER = ['PASS']
    return rec


# Counting samples, variants and AF_by_vark
def proc_line_1st_round(rec, cnf_, self_, variant_dict, control_vars):
    # Strict filter of DP, QUAL, PMEAN
    for f in self_.round1_filters:
        if not f.check(rec):
            return None

    # For those who passed, collect controls, samples and af_by_varid
    samplename = cnf_.name

    var_str = rec.get_variant()  # = tuple (chrom, pos, ref, alt)

    [f.apply(rec) for f in self_.round2_filters]
    min_freq_filter = self_.min_freq_filters[samplename]
    min_freq_filter.apply(rec)

    if samplename and self_.control and samplename == self_.control:
        if not rec.is_rejected() and rec.get_cls() == 'Novel':
            # So that any novel variants showed up in control won't be filtered:
            control_vars.add(var_str)

    if samplename and not rec.is_rejected():
        # Undetermined won't count toward samples
        if 'undetermined' not in samplename.lower() or Filtering.filt_cnf['count_undetermined']:
            if var_str not in variant_dict:
                variant_dict[var_str] = VariantInfo(var_str)
            af = rec.get_af(cnf_.main_sample_index, self_.caller.name)
            variant_dict[var_str].afs.append(af or .0)
        else:
            if not self_.undet_sample_filter.apply(rec):
                err('Undetermined sample for rec ' + rec.get_variant() + ', sample ' + str(samplename))
                return None

    return rec


# Based on counted samples, variants and AF_by_variant:
#   var_n    = number of occurences of variant       must be >= [sample_cnt]
#   fraction = variant_n / number of total samples   must be > [fraction]
#   avg_af   = avg AF for this variant               must be < [freq]
def proc_line_2nd_round(rec, cnf_, self_):
    if rec.is_rejected():
        return rec

    sample_name = rec.sample_field()
    if sample_name:
        var_info = self_.variant_dict.get(rec.get_variant())
        if not var_info:
            return None

        var_n = var_info.occurence_num()
        frac = var_info.frac()
        avg_af = var_info.avg_af()
        rec.INFO['Num_samples'] = len(filtering.sample_names)
        rec.INFO['Num_samples_with_the_same_variant'] = var_n
        rec.INFO['Ave_AF_in_samples_with_the_same_variant'] = avg_af

        if not self_.multi_filter.apply(rec, var_n=var_n, frac=frac, avg_af=avg_af):
            # info('Multi filter: vark = ' + rec.get_variant() +
            #      ', var_n = ' + str(var_n) + ', n_sample = ' +
            #      str(len(self_.sample_names)) + ', avg_af = ' + str(avg_af))
            pass

        if not self_.dup_filter.apply(rec, main_sample_index=cnf_.main_sample_index):
            # info('Dup filter: variant = ' + rec.get_variant())
            pass

        af = rec.get_af(cnf_.main_sample_index, self_.caller.name)
        if af is not None:
            if not self_.max_rate_filter.apply(rec, frac=frac, af=af):
                # info('Max rate filter: variant = ' + rec.get_variant() + ', frac = ' + str(frac) + ', af = ' +
                #      str(rec.af(cnf_.main_sample_index)))
                pass

        self_.control_filter.apply(rec)

        cls = rec.get_cls(req_maf=Filtering.filt_cnf.get('maf'), max_frac_for_del=0.95, sample_name=sample_name)
        rec.INFO['Class'] = cls

        self_.bias_filter.apply(
            rec,
            cls=cls,
            main_sample_index=cnf_.main_sample_index,
            caller_name=self_.caller.name)  # dbSNP, Novel, and bias satisfied - keep
        self_.nonclnsnp_filter.apply(rec, cls=cls)

    return rec


def proc_line_impact(rec, cnf_, self_):
    self_.impact_filter.apply(rec)
    return rec


def proc_line_polymorphic(rec, cnf_, self_):
    self_.polymorphic_filter.apply(rec)
    return rec


# def run_snpsift(cnf, vcf_cnf, vcf_fpath):
#     expression = vcf_cnf.get('expression')
#     if not expression:
#         return vcf_fpath
#
#     step_greetings('Running SnpSift filter.')
#
#     executable = get_java_tool_cmdline(cnf, 'snpsift')
#     cmdline = '{executable} filter -a EXPR -n -p -f ' \
#               '{vcf_fpath} "{expression}"'.format(**locals())
#     filtered_fpath = intermediate_fname(cnf, vcf_fpath, 'snpsift')
#     call(cnf, cmdline, filtered_fpath)
#
#     info('Done.')
#
#     return filtered_fpath


cnfs_for_sample_names = dict()


def prep_vcf(cnf, vcf_fpath, sample_name, caller_name):
    main_sample_index = get_main_sample_index(vcf_fpath, sample_name)

    vcf_fpath = leave_main_sample(cnf, vcf_fpath, sample_name)

    def set_af(rec):
        af, t_ref_count, t_alt_count = None, None, None

        ads = rec.get_val('AD', main_sample_index)
        if ads:
            try:
                t_ref_count, t_alt_count = ads[0], ads[1]
            except:
                t_ref_count, t_alt_count = ads, None

            if caller_name == 'vardict':
                af = rec.get_val('AF', main_sample_index)

            elif caller_name == 'mutect':
                af = rec.get_val('FREQ', main_sample_index)

            elif caller_name == 'freebayes':
                af = rec.get_val('AB', main_sample_index)

            else:
                dp = rec.get_val('DP', main_sample_index=main_sample_index)
                if t_alt_count is not None and dp:
                    af = float(t_alt_count) / dp

            rec.INFO['AF'] = af
        return rec

    vcf_fpath = iterate_vcf(cnf, vcf_fpath, set_af, 'vcf2txt')
    return vcf_fpath


def filter_with_vcf2txt(cnf, bcbio_structure, vcf_fpaths, sample_names, caller, n_threads, sample_min_freq):
    # info()
    # info('Keeping only first sample info in VCFs...')
    # vcf_fpaths = Parallel(n_jobs=n_threads)(delayed(leave_main_sample)(cnf, fpath, s) for fpath, s in zip(vcf_fpaths, sample_names))

    # info()
    # info('Indexing with tabix ')
    # vcf_fpaths = Parallel(n_jobs=n_threads)(delayed(tabix_vcf)(cnf, vcf_fpath) for vcf_fpath in vcf_fpaths)

    safe_mkdir(bcbio_structure.var_dirpath)

    # combined_vcf_fpath = join(bcbio_structure.var_dirpath, caller.name + '.vardict.vcf')
    #
    # info()
    # info('Merging VCFs...')
    # vcf_merge(cnf, vcf_fpaths, combined_vcf_fpath)

    info()
    info('Preparing VCFs for vcf2txt')
    vcf_fpaths = Parallel(n_jobs=n_threads)(delayed(prep_vcf)(cnf, fpath, s, caller.name)
        for fpath, s in zip(vcf_fpaths, sample_names))

    caller.vcf2txt_res_fpath = join(bcbio_structure.var_dirpath, caller.name + '.txt')
    res = run_vcf2txt(cnf, vcf_fpaths, caller.vcf2txt_res_fpath, bcbio_structure.paired, sample_min_freq)
    if not res:
        err('vcf2txt run returned non-0')
        return None

    res = run_pickline(cnf, caller, caller.vcf2txt_res_fpath)
    if not res:
        err('pickLine run returned non-0')
        return None
    return res


def run_pickline(cnf, caller, vcf2txt_res_fpath):
    pick_line = get_script_cmdline(cnf, 'perl', join('external', 'pickLine.pl'))
    if not pick_line:
        sys.exit(1)
        return None

    caller.pickline_res_fpath = add_suffix(vcf2txt_res_fpath, 'PASS')

    cmdline = '{pick_line} -l PASS:TRUE -c 45 {vcf2txt_res_fpath} | grep -vw dbSNP | ' \
              'grep -v UTR_ | grep -vw SILENT | grep -v intron_variant | grep -v upstream_gene_variant | ' \
              'grep -v downstream_gene_variant | grep -v intergenic_region | grep -v intragenic_variant | ' \
              'grep -v NON_CODING'

# grep -E "(^Sample)|(HIGH)|(MODERATE)" variants.txt | pickLine -l PASS:TRUE -c 45 | grep -vw dbSNP | pickLine -v -i 12:3 -c 14:11

    if cnf.genome.polymorphic_variants:
        poly_vars = abspath(cnf.genome.polymorphic_variants)
        cmdline += ' | {pick_line} -v -i 12:3 -c 14:11 {poly_vars}'
    cmdline = cmdline.format(**locals())
    res = call(cnf, cmdline, caller.pickline_res_fpath, exit_on_error=False)
    if not res:
        return None
    else:
        return res


def filter_for_variant_caller(caller, cnf, bcbio_structure):
    IN_PARALLEL = False

    info('Running for ' + caller.name)

    vcf_by_sample = caller.find_anno_vcf_by_sample(optional_ext='.gz')

    if len(vcf_by_sample) == 0:
        err('No vcfs for ' + caller.name + '. Skipping.')
        return caller

    info('*' * 70)
    info('Samples (total ' + str(len(vcf_by_sample)) + '):')

    # Set up cnfs for each sample
    # global cnfs_for_sample_names
    # for sample in caller.samples:
    #     info(sample.name)
    #     cnf_copy = cnf.copy()
    #     cnf_copy['name'] = sample.name
    #     if cnf_copy['variant_filtering'].get('min_freq') is not None:
    #         sample.min_af = cnf_copy['variant_filtering']['min_freq']
    #     elif sample.min_af is None:
    #         sample.min_af = 0
    #     cnfs_for_sample_names[sample.name] = cnf_copy

    n_threads = cnf.threads if cnf.threads else min(10, len(vcf_by_sample) + 1)
    if not IN_PARALLEL: n_threads = 1
    info('Number of threads for filtering: ' + str(n_threads))

    vcf_fpaths = vcf_by_sample.values()
    sample_names = vcf_by_sample.keys()

    # global filtering
    # filtering = Filtering(cnf, bcbio_structure, caller)

    # info('Filtering by impact')
    # vcf_fpaths = Parallel(n_jobs=n_threads)(delayed(impact_round)(vcf_fpath, s) for vcf_fpath, s in zip(vcf_fpaths, sample_names))

    info()
    info('-' * 70)
    info('Filtering using vcf2txt...')
    res = filter_with_vcf2txt(cnf, bcbio_structure, vcf_fpaths, sample_names, caller, n_threads, caller.samples[0].min_af)
    if not res:
        return None

    # symlinking
    pass_comb_basefname = basename(caller.pickline_res_fpath)
    comb_pass_maf_fpath_symlink = join(bcbio_structure.date_dirpath, pass_comb_basefname)
    if not exists(comb_pass_maf_fpath_symlink) \
            and not islink(comb_pass_maf_fpath_symlink) \
            and caller.pickline_res_fpath != comb_pass_maf_fpath_symlink:
        os.symlink(caller.pickline_res_fpath, comb_pass_maf_fpath_symlink)

    info('-' * 70)
    info()


    # info()
    # info('-' * 70)
    # info('Filtering using reengineered stuff...')
    # vcf_fpaths__intrinsic = filtering.run_regengineered_filtering(vcf_fpaths, sample_names, n_threads)
    # postprocess_filtered_vcfs(cnf, bcbio_structure, vcf_fpaths__intrinsic, sample_names, caller, n_threads)

    info('-' * 70)
    info()

    return caller


def combine_mafs(cnf, maf_fpaths, output_basename):
    output_fpath = output_basename + '.maf'
    output_pass_fpath = output_basename + '.pass.maf'

    if isfile(output_fpath): os.remove(output_fpath)
    if isfile(output_pass_fpath): os.remove(output_pass_fpath)

    if not maf_fpaths:
        warn('No MAFs - no combined MAF will be made.')
        return None, None

    if not isdir(dirname(output_fpath)): safe_mkdir(dirname(output_fpath))

    with open(output_fpath, 'w') as out, \
         open(output_pass_fpath, 'w') as out_pass:

        for i, fpath in enumerate(maf_fpaths):
            with open(fpath) as inp:
                for j, line in enumerate(inp):
                    if i > 0 and j in [0, 1]:
                        continue
                    out.write(line)
                    if '\tInvalid\t' not in line:
                        out_pass.write(line)
    return output_fpath, output_pass_fpath


def run_vcf2txt(cnf, vcf_fpaths, final_maf_fpath, paired, sample_min_freq=None):
    info()
    info('Running VarDict vcf2txt...')

    vcf2txt = None
    # if paired:
    #    vcf2txt = get_script_cmdline(cnf, 'perl', join('external', 'vcf2txt_paired.pl'))
    # else:
    #    vcf2txt = get_script_cmdline(cnf, 'perl', join('external', 'vcf2txt.pl'))
    vcf2txt = get_script_cmdline(cnf, 'perl', join('external', 'vcf2txt.pl'))

    if not vcf2txt:
        sys.exit(1)

    vcf_fpaths = ' '.join(vcf_fpaths)

    c = cnf.variant_filtering
    min_freq = c.min_freq or sample_min_freq or defaults.default_min_freq

    cmdline = '{vcf2txt} ' \
        '-F {min_freq} -n {c.sample_cnt} -f {c.freq} -p {c.min_p_mean} -q {c.min_q_mean} ' \
        '-r {c.fraction} -R {c.max_ratio} -P {c.filt_p_mean} -Q {c.filt_q_mean} -D {c.filt_depth} ' \
        '-M {c.min_mq} -V {c.min_vd} -G {c.maf} -o {c.signal_noise} ' \
        '{vcf_fpaths}'.format(**locals())

    if cnf.genome.dbsnp_multi_mafs:
        cmdline += ' -A ' + cnf.genome.dbsnp_multi_mafs

    res = call(cnf, cmdline, final_maf_fpath, exit_on_error=False)
    return res


def postprocess_vcf(sample, caller_name, work_filt_vcf_fpath):
    if not work_filt_vcf_fpath:
        err(sample.name + ': work_filt_vcf_fpath is None')
        return None, None, None

    cnf = cnfs_for_sample_names.get(sample.name)
    if cnf is None:
        err('Error: for ' + sample.name + ': cnf is None')
        return None, None, None
    # work_filt_vcf_fpath = leave_first_sample(cnf, work_filt_vcf_fpath)

    safe_mkdir(sample.get_filtered_vcfs_dirpath())
    final_vcf_fpath = join(sample.get_filtered_vcfs_dirpath(), sample.name + '-' + caller_name + '.anno.filt.vcf')

    file_basefpath = splitext(final_vcf_fpath)[0]
    final_vcf_fpath = file_basefpath + '.vcf'
    final_tsv_fpath = file_basefpath + '.tsv'
    final_maf_fpath = file_basefpath + '.maf'
    final_pass_vcf_fpath = join(dirname(final_vcf_fpath), sample.name + '-' + caller_name + \
        BCBioStructure.pass_filt_vcf_ending)  # for futrher processing

    BCBioStructure.move_vcfs_to_var(sample)

    # Moving final VCF
    if isfile(final_vcf_fpath): os.remove(final_vcf_fpath)
    # shutil.move(work_filt_vcf_fpath, final_vcf_fpath)
    # os.symlink(final_vcf_fpath, work_filt_vcf_fpath)
    shutil.copy(work_filt_vcf_fpath, final_vcf_fpath)
    igvtools_index(cnf, final_vcf_fpath)

    # Cleaning rejected variants
    pass_filtered_vcf_fpath = remove_rejected(cnf, work_filt_vcf_fpath)
    if vcf_is_empty(cnf, pass_filtered_vcf_fpath):
        info('All variants are rejected.')
    if isfile(final_pass_vcf_fpath): os.remove(final_pass_vcf_fpath)
    if islink(pass_filtered_vcf_fpath): os.unlink(pass_filtered_vcf_fpath)
    shutil.copy(pass_filtered_vcf_fpath, final_pass_vcf_fpath)
    # os.symlink(final_clean_vcf_fpath, clean_filtered_vcf_fpath)
    igvtools_index(cnf, final_pass_vcf_fpath)

    # Converting to TSV
    if work_filt_vcf_fpath and 'tsv_fields' in cnf:
        tsv_fpath = make_tsv(cnf, work_filt_vcf_fpath, sample.name)
        if isfile(final_tsv_fpath): os.remove(final_tsv_fpath)
        if not tsv_fpath:
            err('TSV convertion didn\'t work for ' + sample.name)
            final_tsv_fpath = None
        else:
            shutil.copy(tsv_fpath, final_tsv_fpath)
    else:
        if isfile(final_tsv_fpath): os.remove(final_tsv_fpath)
        final_tsv_fpath = None

    # Converting clean VCF to TSV
    # if clean_filtered_vcf_fpath and 'tsv_fields' in cnf:
    #     clean_tsv_fpath = make_tsv(cnf, clean_filtered_vcf_fpath)
    #
    #     if isfile(final_clean_tsv_fpath):
    #         os.remove(final_clean_tsv_fpath)
    #     shutil.move(clean_tsv_fpath, final_clean_tsv_fpath)
    # else:
    #     final_clean_tsv_fpath = None

    # Converting to MAF
    if cnf.make_maf and work_filt_vcf_fpath:
        maf_fpath = convert_to_maf(
            cnf, work_filt_vcf_fpath,
            tumor_sample_name=sample.name,
            bam_fpath=sample.bam,
            transcripts_fpath=cnf.transcripts_fpath,
            normal_sample_name=sample.normal_match.name if sample.normal_match else None)
        if isfile(final_maf_fpath): os.remove(final_maf_fpath)
        if maf_fpath:
            shutil.copy(maf_fpath, final_maf_fpath)
        else:
            final_maf_fpath = None
        info('-' * 70)
        info()
    else:
        if isfile(final_maf_fpath): os.remove(final_maf_fpath)
        final_maf_fpath = None

    return [final_vcf_fpath, final_tsv_fpath, final_maf_fpath]
