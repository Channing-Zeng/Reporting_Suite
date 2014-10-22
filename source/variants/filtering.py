from genericpath import isdir
import os
import shutil
import pickle
import operator

from os.path import basename, join, isfile, dirname, splitext, islink
from joblib import Parallel, delayed
##from memory_profiler import profile

from source.bcbio_structure import BCBioStructure
from source.variants.Effect import Effect
from source.logger import step_greetings, info, critical, err
from source.variants.anno import _snpsift_annotate
from source.variants.vcf_processing import iterate_vcf, vcf_one_per_line, leave_main_sample, get_trasncripts_fpath, \
    get_main_sample_index
from source.utils import mean
from source.file_utils import safe_mkdir, add_suffix, verify_file
from source.variants.tsv import make_tsv
from source.variants.vcf_processing import remove_rejected, convert_to_maf, vcf_is_empty, igvtools_index
from source.logger import info


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
                err('Warning: EFF field is missing for variant ' + str(rec.get_variant()))
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


cnf_for_samples = dict()


# @profile
def process_vcf(vcf_fpath, fun, suffix, cnf=None, *args, **kwargs):
    return iterate_vcf(cnf, vcf_fpath, fun, suffix, cnf_=cnf, self_=filtering, *args, **kwargs)


def rm_prev_round(vcf_fpath, sample_name):
    return process_vcf(vcf_fpath, proc_line_remove_prev_filter, 'rm_prev', cnf_for_samples[sample_name])


def first_round(vcf_fpath, sample_name):
    variant_dict = dict()  # tuples of (chrom, pos, ref, alt)
    control_variants = set()  # variants from controls samples
    res = process_vcf(
        vcf_fpath, proc_line_1st_round, 'r1', cnf_for_samples[sample_name],
        variant_dict=variant_dict, control_vars=control_variants)

    return res, variant_dict, control_variants


def second_round(vcf_fpath, sample_name):
    cnf = cnf_for_samples[sample_name]
    main_sample_index = get_main_sample_index(cnf, vcf_fpath, sample_name)
    cnf.main_sample_index = main_sample_index

    #TODO: tmp
    res = _snpsift_annotate(cnf, cnf.get('clinvar'), 'clinvar', vcf_fpath)
    if not res:
        err('Could not annotate with clinvar.')

    return process_vcf(vcf_fpath, proc_line_2nd_round, 'r2', cnf)


def impact_round(vcf_fpath, sample_name):
    return process_vcf(vcf_fpath, proc_line_impact, 'impact', cnf_for_samples[sample_name])


def one_per_line(vcf_fpath, sample_name):
    return vcf_one_per_line(cnf_for_samples[sample_name], vcf_fpath)


class Filtering:
    cnf = None
    filt_cnf = None

    def __init__(self, cnf, bcbio_structure, caller, cnf_for_samples):
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
            self.min_freq_filters[sample.name] = \
                InfoFilter('min_freq', 'AF', cnf=cnf_for_samples[sample.name], required=False)

        self.round2_filters = [
            InfoFilter('min_p_mean', 'PMEAN', required=False),
            InfoFilter('min_q_mean', 'QMEAN', required=False),
            InfoFilter('mean_mq', 'MQ', required=False),
            InfoFilter('signal_noise', 'SN', required=False),
            InfoFilter('mean_vd', 'VD', required=False)]

        def dup_filter_check(rec, main_sample_index):
            pstd = rec.get_val('PSTD')
            bias = rec.bias(main_sample_index)

            # all variants from one position in reads
            if pstd is not None and bias is not None:
                return not (pstd == 0 and bias[-1] != '0' and bias[-1] != '1')
            return True

        def bias_filter_check(rec, cls, main_sample_index):  # Filter novel variants with strand bias.
            if not Filter.filt_cnf['bias'] is True:
                return True  # I we don't need to check bias, just keep the variant

            if not cls in ['Novel', 'dbSNP']:
                return True  # Check only for novel and dnSNP

            if not (rec.bias(main_sample_index) and rec.bias(main_sample_index) in ['2:1', '2:0']):
                return True  # Variant has to have bias and with proper values

            return rec.af(main_sample_index) >= 0.3

        def nonclnsnp_filter_check(rec, cls):
            if cls in ['COSMIC']:
                return True

            return rec.check_clnsig() != 'not_significant'

        def multi_filter_check(rec, var_n, frac, avg_af):
            if (rec.ID is None and      # reject if novel and present in [fraction] samples
                frac > Filter.filt_cnf['fraction'] and
                var_n >= Filter.filt_cnf['sample_cnt'] and
                avg_af < Filter.filt_cnf['freq']):
                return False

        def max_rate_filter_check(rec, frac):  # reject if present in [max_ratio] samples
            if frac >= Filter.filt_cnf['max_ratio'] and rec.get_val('AF') < 0.3:
                return False

        # self.undet_sample_filter = Filter('UNDET_SAMPLE', lambda rec: False)
        self.control_filter = CnfFilter('control', lambda rec:
            not (filt_cnf['control'] and rec.get_variant() in self.control_vars))
        self.dup_filter = Filter('DUP', dup_filter_check)
        self.bias_filter = CnfFilter('bias', bias_filter_check)
        self.nonclnsnp_filter = Filter('NonClnSNP', nonclnsnp_filter_check)
        self.multi_filter = Filter('MULTI', multi_filter_check)
        self.max_rate_filter = CnfFilter('max_ratio', max_rate_filter_check)


    # @profile
    def run_filtering(self, sample_names, vcf_fpaths, n_jobs=1):
        step_greetings('Filtering')

        global filtering
        filtering = self

        info('Fixing previous . values to PASS')
        vcf_fpaths = Parallel(n_jobs=n_jobs)(delayed(rm_prev_round)(vcf_fpath, s)
                                             for vcf_fpath, s in zip(vcf_fpaths, sample_names))
        info()

        info('First round')
        results = Parallel(n_jobs=n_jobs)(delayed(first_round)(vcf_fpath, s)
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

            # if isfile(varks_dump_fpath):
            #     try:
            #         with open(varks_dump_fpath) as f:
            #             varks_2 = pickle.load(f)
            #         print 'Varks restored (to check): len=', str(len(varks_2.keys()))
            #         print 'Varks new                : len=', str(len(self.varks.keys()))
            #     except:
            #         pass

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
        vcf_fpaths = Parallel(n_jobs=n_jobs)(delayed(one_per_line)(vcf_fpath, s)
                     for vcf_fpath, s in zip(vcf_fpaths, sample_names))
        info()

        info('Second round')
        vcf_fpaths = Parallel(n_jobs=n_jobs)(delayed(second_round)(vcf_fpath, s)
                     for vcf_fpath, s in zip(vcf_fpaths, sample_names))
        info()

        info('Filtering by impact')
        vcf_fpaths = Parallel(n_jobs=n_jobs)(delayed(impact_round)(vcf_fpath, s)
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

    all_passed = all([f.apply(rec) for f in self_.round1_filters])
    if not all_passed:
        return None

    # For those who passed, collect controls, samples and af_by_varid
    samplename = cnf_.name

    var_str = rec.get_variant()  # = tuple (chrom, pos, ref, alt)

    [f.apply(rec) for f in self_.round2_filters]

    min_freq_filter = self_.min_freq_filters[samplename]
    min_freq_filter.apply(rec)

    if samplename and self_.control and samplename == self_.control:
        if not rec.is_rejected() or rec.cls() == 'Novel':
            # So that any novel variants showed up in control won't be filtered:
            control_vars.add(var_str)

    if samplename and not rec.is_rejected():
        # Undetermined won't count toward samples
        if 'undetermined' not in samplename.lower() or Filtering.filt_cnf['count_undetermined']:
            if var_str not in variant_dict:
                variant_dict[var_str] = VariantInfo(var_str)
            variant_dict[var_str].afs.append(rec.get_val('AF', .0))
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

        if not self_.multi_filter.apply(rec, var_n=var_n, frac=frac, avg_af=avg_af):
            info('Multi filter: vark = ' + rec.get_variant() +
                 ', var_n = ' + str(var_info.occurence_num()) + ', n_sample = ' +
                 str(len(self_.sample_names)) + ', avg_af = ' + str(var_info.avg_af()))

        if not self_.dup_filter.apply(rec, main_sample_index=cnf_.main_sample_index):
            info('Dup filter: variant = ' + rec.get_variant())

        if rec.af(cnf_.main_sample_index) is not None:
            if not self_.max_rate_filter.apply(rec, frac=frac):
                info('Max rate filter: variant = ' + rec.get_variant() + ', frac = ' + str(frac) + ', af = ' +
                     str(rec.af(cnf_.main_sample_index)))

        self_.control_filter.apply(rec)

        cls = rec.cls()

        # We check for caf - and if it is above the required value, we assume that the variant class is dnSNP regardless of Cosmic. The class is needed for further CLNSIG and BIAS filters.
        # Then we check if the calss is deleterious dnSNP
        # Then we filter for strand bias only variants with Novel and dbSNP class
        # And then filter out all dbSNP variants that are significant according to Clinvar
        if 'CAF' in rec.INFO:
            vals = rec.INFO['CAF']
            
            cafs = [''.join(c for c in v if c not in '[]') for v in vals if v]
            print cafs
            if len(cafs) == 0:
                print 'cafs = ' + str(cafs) + ', vals = ' + str(vals) + ' for ' + rec.get_variant() + ' in ' + sample_name
            else:
                allele_cafs = [c for c in cafs[1:] if c]
                if len(allele_cafs) == 0:
                    print 'cafs = ' + str(cafs) + ', allele_cafs = ' + str(allele_cafs) + ' for ' + rec.get_variant() + ' in ' + sample_name
                else:
                    allele_cafs = map(float, allele_cafs)
                    min_allele_caf = min(allele_cafs)
                    req_maf = Filtering.filt_cnf.get('maf')

                    # if there's MAF with frequency, it'll be considered
                    # dbSNP regardless of COSMIC
                    if req_maf is not None and min_allele_caf > req_maf:
                        info('min_allele_caf = ' + str(min_allele_caf) + ', req_maf = ' + str(req_maf) + ', class was ' + cls + ', becomes dbSNP')
                        cls = 'dbSNP'
                    else:
                        info('min_allele_caf = ' + str(min_allele_caf) + ', req_maf = ' + str(req_maf) + ', class stays ' + cls)

        # Rescue deleterious dbSNP, such as rs80357372 (BRCA1 Q139) that is in dbSNP,
        # but not in ClnSNP or COSMIC.
        for eff in map(Effect, rec.INFO.get('EFF') or []):
            if eff.efftype in ['STOP_GAINED', 'FRAME_SHIFT'] and cls == 'dbSNP':
                if eff.pos / int(eff.aal) < 0.95:
                    cls = 'dbSNP_del'



        self_.bias_filter.apply(rec, cls=cls, main_sample_index=cnf_.main_sample_index)  # dbSNP, Novel, and bias satisfied - keep
        self_.nonclnsnp_filter.apply(rec, cls=cls)  # significant and not Cosmic - keep

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


# @profile
def filter_for_variant_caller(caller, cnf, bcbio_structure):
    IN_PARALLEL = True

    info('Running for ' + caller.name)

    anno_vcf_by_sample = caller.find_anno_vcf_by_sample()
    sample_names = anno_vcf_by_sample.keys()
    anno_vcf_fpaths = anno_vcf_by_sample.values()

    if len(anno_vcf_fpaths) == 0:
        err('No vcfs for ' + caller.name + '. Skipping.')
        return caller

    global cnf_for_samples
    for sample in caller.samples:
        cnf_for_samples[sample.name] = cnf.copy()
        cnf_for_samples[sample.name].name = sample.name

    f = Filtering(cnf, bcbio_structure, caller, cnf_for_samples)

    n_jobs = min(len(anno_vcf_fpaths), 20) if IN_PARALLEL else 1
    filt_anno_vcf_fpaths = f.run_filtering(sample_names, anno_vcf_fpaths, n_jobs)

    global cnfs_for_sample_names
    info('*' * 70)
    info('Processed samples (' + str(len(sample_names)) + '):')
    for sample in caller.samples:
        info(sample.name)
        cnf_copy = cnf.copy()
        cnf_copy['name'] = sample.name
        if sample.min_af is not None:
            cnf_copy['variant_filtering']['min_freq'] = sample.min_af
        cnfs_for_sample_names[sample.name] = cnf_copy
        # TODO: pass min_freq to filters

    results = [r for r in Parallel(n_jobs) \
        (delayed(postprocess_vcf)
         (sample, anno_vcf_fpath, work_filt_vcf_fpath)
              for sample, anno_vcf_fpath, work_filt_vcf_fpath in
              zip(caller.samples, anno_vcf_fpaths, filt_anno_vcf_fpaths)
         ) if r is not None and None not in r]
    info('Results: ' + str(len(results)))
    info('*' * 70)

    for sample, [vcf, tsv, maf] in zip(caller.samples, results):
        info('Sample ' + sample.name + ': ' + vcf + ', ' + tsv + ', ' + maf)
        # sample.filtered_vcf_by_callername[caller.name] = vcf
        # sample.filtered_tsv_by_callername[caller.name] = tsv
        # sample.filtered_maf_by_callername[caller.name] = maf

    comb_maf_fpath = join(bcbio_structure.var_dirpath, caller.name + '.maf')
    caller.combined_filt_maf_fpath = combine_mafs(cnf, caller.find_filt_maf_by_sample().values(), comb_maf_fpath)
    info('-' * 70)
    info()

    return caller


def combine_mafs(cnf, maf_fpaths, output_fpath):
    if isfile(output_fpath):
        os.remove(output_fpath)

    if not isdir(dirname(output_fpath)):
        safe_mkdir(dirname(output_fpath))

    with open(output_fpath, 'w') as out:
        for i, fpath in enumerate(maf_fpaths):
            with open(fpath) as inp:
                for j, line in enumerate(inp):
                    if i > 0 and j in [0, 1]:
                        continue
                    out.write(line)
    return output_fpath


def postprocess_vcf(sample, original_anno_vcf_fpath, work_filt_vcf_fpath):
    if not original_anno_vcf_fpath or not work_filt_vcf_fpath:
        if not original_anno_vcf_fpath:
            err(sample.name + ': original_anno_vcf_fpath is None')
        if not work_filt_vcf_fpath:
            err(sample.name + ': work_filt_vcf_fpath is None')
        return None, None, None

    cnf = cnfs_for_sample_names.get(sample.name)
    if cnf is None:
        err('Error: for ' + sample.name + ': cnf is None')
        return None, None, None
    # work_filt_vcf_fpath = leave_first_sample(cnf, work_filt_vcf_fpath)

    final_vcf_fpath = add_suffix(original_anno_vcf_fpath, 'filt').replace('varAnnotate', 'varFilter')
    safe_mkdir(dirname(final_vcf_fpath))

    file_basepath = splitext(final_vcf_fpath)[0]
    final_vcf_fpath = file_basepath + '.vcf'
    final_tsv_fpath = file_basepath + '.tsv'
    final_maf_fpath = file_basepath + '.maf'
    final_clean_vcf_fpath = file_basepath + '.passed.vcf'  # for futrher processing

    BCBioStructure.move_vcfs_to_var(sample)

    # Moving final VCF
    if isfile(final_vcf_fpath): os.remove(final_vcf_fpath)
    # shutil.move(work_filt_vcf_fpath, final_vcf_fpath)
    # os.symlink(final_vcf_fpath, work_filt_vcf_fpath)
    shutil.copy(work_filt_vcf_fpath, final_vcf_fpath)
    igvtools_index(cnf, final_vcf_fpath)

    # Cleaning rejected variants
    clean_filtered_vcf_fpath = remove_rejected(cnf, work_filt_vcf_fpath)
    if vcf_is_empty(cnf, clean_filtered_vcf_fpath):
        info('All variants are rejected.')
    if isfile(final_clean_vcf_fpath): os.remove(final_clean_vcf_fpath)
    if islink(clean_filtered_vcf_fpath): os.unlink(clean_filtered_vcf_fpath)
    shutil.copy(clean_filtered_vcf_fpath, final_clean_vcf_fpath)
    # os.symlink(final_clean_vcf_fpath, clean_filtered_vcf_fpath)
    igvtools_index(cnf, final_clean_vcf_fpath)

    # Converting to TSV
    if work_filt_vcf_fpath and 'tsv_fields' in cnf:
        tsv_fpath = make_tsv(cnf, work_filt_vcf_fpath, sample.name)
        if not tsv_fpath:
            err('TSV convertion didn\'t work for ' + sample.name)
            final_tsv_fpath = None
        else:
            if isfile(final_tsv_fpath): os.remove(final_tsv_fpath)
        shutil.copy(tsv_fpath, final_tsv_fpath)
    else:
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
    if work_filt_vcf_fpath and cnf.make_maf:
        maf_fpath = convert_to_maf(
            cnf, work_filt_vcf_fpath,
            tumor_sample_name=sample.name,
            bam_fpath=sample.bam,
            transcripts_fpath=cnf.transcripts_fpath,
            normal_sample_name=sample.normal_match.name if sample.normal_match else None)
        # if isfile(final_maf_fpath): os.remove(final_maf_fpath)
        shutil.copy(maf_fpath, final_maf_fpath)
        info('-' * 70)
        info()
    else:
        final_maf_fpath = None

    return [final_vcf_fpath, final_tsv_fpath, final_maf_fpath]
