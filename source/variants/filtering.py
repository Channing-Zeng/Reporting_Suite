from genericpath import isdir
import os
import shutil
import pickle
import operator

from os.path import basename, join, isfile, dirname, splitext, islink
from joblib import Parallel, delayed
##from memory_profiler import profile

from source.variants.Effect import Effect
from source.logger import step_greetings, info, critical, err
from source.variants.vcf_processing import iterate_vcf, vcf_one_per_line, leave_first_sample
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

        self.num_passed = 0
        self.num_rejected = 0

    @staticmethod
    def __remove_pass(rec):
        if rec.FILTER == ['PASS']:
            rec.FILTER = None

    def apply(self, rec, only_check=False, *args, **kvargs):
        if only_check:
            return self.check(rec, *args, **kvargs)

        if self.check(rec, *args, **kvargs):  # True if PASS, False otherwise
            self.num_passed += 1
            return True
        else:
            self.num_rejected += 1

            if self.word not in rec.FILTER:
                self.__remove_pass(rec)
                rec.add_filter(self.word)
            return False


class CnfFilter(Filter):
    def __init__(self, key, check, *args, **kvargs):
        self.key = key
        Filter.__init__(self, key.upper(), check, *args, **kvargs)

    def apply(self, rec, only_check=False, *args, **kvargs):
        if Filter.filt_cnf.get(self.key) is None:
            return True

        return Filter.apply(self, rec, only_check, *args, **kvargs)


class InfoFilter(CnfFilter):
    def __init__(self, cnf_key, info_key, op=operator.gt, *args, **kwargs):
        # True if PASS
        def check(rec):
            anno_val = rec.get_val(info_key)
            if anno_val is None:
                if self.required:
                    critical('Error: no field ' + info_key + ' in INFO or SAMPLE for variant ' +
                             rec.var_id() + ' - required to test ' + cnf_key)
                else:
                    return True  # PASS
            try:
                cnf_val = float(Filter.filt_cnf[cnf_key])
            except ValueError:
                cnf_val = int(Filter.filt_cnf[cnf_key])

            return op(anno_val, cnf_val)

        CnfFilter.__init__(self, cnf_key, check, *args, **kwargs)


class EffectFilter(CnfFilter):
    def __init__(self, cnf_key, *args, **kwargs):
        def check(rec):
            if 'EFF' not in rec.INFO:
                err('Warning: EFF field is missing for variant ' + str(rec.var_id()))
                return False

            impact_filter_line = Filter.filt_cnf[cnf_key]
            if impact_filter_line:
                important_impacts = [s.upper() for s in impact_filter_line.split('|')]

                return any(eff.impact.upper() in important_impacts
                       for eff in map(Effect, rec.INFO['EFF']))

        CnfFilter.__init__(self, cnf_key, check, *args, **kwargs)


class VarkInfo:
    def __init__(self, vark):
        self.vark = vark
        self.afs = []
        self._avg_af = None

    def var_n(self):
        return len(self.afs)

    def frac(self):
        return float(self.var_n()) / len(filtering.sample_names)

    def avg_af(self):
        if self._avg_af is None:
            self._avg_af = mean(self.afs)
        return self._avg_af


cnf_for_samples = dict()


# @profile
def process_vcf(vcf_fpath, fun, suffix, *args, **kwargs):
    cnf = cnf_for_samples[basename(vcf_fpath).split('.')[0]]
    return iterate_vcf(cnf, vcf_fpath, fun, suffix, self_=filtering, *args, **kwargs)


def rm_prev_round(vcf_fpath):
    return process_vcf(vcf_fpath, proc_line_remove_prev_filter, 'rm_prev')


def first_round(vcf_fpath):
    varks = dict()
    control_vars = set()
    res = process_vcf(vcf_fpath, proc_line_1st_round, 'r1',
                      varks=varks, control_vars=control_vars)

    return res, varks, control_vars


def second_round(vcf_fpath):
    return process_vcf(vcf_fpath, proc_line_2nd_round, 'r2')


def impact_round(vcf_fpath):
    return process_vcf(vcf_fpath, proc_line_impact, 'impact')


def one_per_line(vcf_fpath):
    cnf = cnf_for_samples[basename(vcf_fpath).split('.')[0]]
    return vcf_one_per_line(cnf, vcf_fpath)


class Filtering:
    cnf = None
    filt_cnf = None

    def __init__(self, cnf, bcbio_structure, caller):
        Filtering.cnf = cnf
        Filtering.filt_cnf = Filter.filt_cnf = filt_cnf = cnf.variant_filtering.__dict__

        self.caller = caller
        self.control_vars = set()
        self.sample_names = set([s.name for s in caller.samples])
        self.varks = dict()  # vark -> VarkInfo(vark, afs)
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

        self.round2_filters = [
            InfoFilter('min_p_mean', 'PMEAN'),
            InfoFilter('min_q_mean', 'QUAL'),
            InfoFilter('min_freq', 'AF'),
            InfoFilter('mean_mq', 'MQ'),
            InfoFilter('signal_noise', 'SN'),
            InfoFilter('mean_vd', 'VD', required=False)]

        def dup_filter_check(rec):
            pstd = rec.get_val('PSTD')
            bias = rec.bias()

            # all variants from one position in reads
            if pstd is not None and bias is not None:
                return not (pstd == 0 and bias[-1] != '0' and bias[-1] != '1')
            return True

        def bias_filter_check(rec, cls):
            return not (  # Filter novel variants with strand bias.
                Filter.filt_cnf['bias'] is True and
                cls in ['Novel', 'dbSNP'] and
                rec.bias() and rec.bias() in ['2:1', '2:0'] and rec.af() < 0.3)

        def nonclnsnp_filter_check(rec, cls):
            return rec.check_clnsig() != -1 or cls == 'COSMIC'

        def multi_filter_check(rec, var_n, frac, avg_af):
            return not (  # reject if novel and present in [fraction] samples
                rec.ID is None and
                frac > Filter.filt_cnf['fraction'] and
                var_n >= Filter.filt_cnf['sample_cnt'] and
                avg_af < Filter.filt_cnf['freq'])

        def max_rate_filter_check(rec, frac):  # reject if present in [max_ratio] samples
            return not (
                frac >= Filter.filt_cnf['max_ratio'] and
                rec.get_val('AF') < 0.3)

        self.undet_sample_filter = Filter('UNDET_SAMPLE', lambda rec: rec.var_id() in self.varks)
        self.control_filter = CnfFilter('control', lambda rec: not (filt_cnf['control'] and rec.var_id() in self.control_vars))
        self.dup_filter = Filter('DUP', dup_filter_check)
        self.bias_filter = CnfFilter('bias', bias_filter_check)
        self.nonclnsnp_filter = Filter('NonClnSNP', nonclnsnp_filter_check)
        self.multi_filter = Filter('MULTI', multi_filter_check)
        self.max_rate_filter = CnfFilter('max_ratio', max_rate_filter_check)


    # @profile
    def run_filtering(self, vcf_fpaths):
        step_greetings('Filtering')

        info('Removing previous FILTER values')

        n_jobs = len(vcf_fpaths)
        # n_jobs = 1

        global cnf_for_samples, filtering
        filtering = self
        for vcf_fpath in vcf_fpaths:
            cnf_for_samples[basename(vcf_fpath).split('.')[0]] = Filtering.cnf.copy()

        vcf_fpaths = Parallel(n_jobs=n_jobs)(delayed(rm_prev_round)(vcf_fpath) for vcf_fpath in vcf_fpaths)
        info()

        info('First round')
        results = Parallel(n_jobs=n_jobs)(delayed(first_round)(vcf_fpath) for vcf_fpath in vcf_fpaths)

        control_vars_dump_fpath = join(self.cnf.work_dir, self.caller.name + '_control_vars.txt')
        varks_dump_fpath = join(self.cnf.work_dir, self.caller.name + 'varks.txt')

        if not [all([varks, control_vars]) for (_, varks, control_vars) in results]:
            info('Skipped first round, restoring varks and controls vars.')
            if not verify_file(control_vars_dump_fpath) or not verify_file(varks_dump_fpath):
                critical('Cannot restore varks and control_vars, please, run without the --reuse flag.')

            with open(control_vars_dump_fpath) as f:
                info('Loading control vars...')
                self.control_vars = pickle.load(f)
                info('Loaded control vars: ' + str(len(self.control_vars)))

            with open(varks_dump_fpath) as f:
                info('Loading varks...')
                self.varks = pickle.load(f)
                info('Loaded varks: ' + str(len(self.varks)))

        else:
            for vcf_fpath, varks, control_vars in results:
                for vark, vark_info in varks.items():
                    if vark not in self.varks:
                        self.varks[vark] = vark_info
                    else:
                        self.varks[vark].afs.extend(vark_info.afs)

                self.control_vars.update(control_vars)

            if isfile(varks_dump_fpath):
                try:
                    with open(varks_dump_fpath) as f:
                        varks_2 = pickle.load(f)
                    print 'Varks restored (to check): len=', str(len(varks_2.keys()))
                    print 'Varks new                : len=', str(len(self.varks.keys()))
                except:
                    pass

            with open(control_vars_dump_fpath, 'w') as f:
                info('Saving control vars...')
                pickle.dump(self.control_vars, f)
                info('Saved control vars: ' + str(len(self.control_vars)))
            with open(varks_dump_fpath, 'w') as f:
                info('Saving varks...')
                pickle.dump(self.varks, f)
                info('Saved varks: ' + str(len(self.varks)))

        vcf_fpaths = [vcf_fpath for vcf_fpath, _, _ in results]
        filtering = self
        info()

        info('One effect per line')
        vcf_fpaths = Parallel(n_jobs=n_jobs)(delayed(one_per_line)(vcf_fpath) for vcf_fpath in vcf_fpaths)

        info('Second round')
        vcf_fpaths = Parallel(n_jobs=n_jobs)(delayed(second_round)(vcf_fpath) for vcf_fpath in vcf_fpaths)
        info()

        info('Filtering by impact')
        vcf_fpaths = Parallel(n_jobs=n_jobs)(delayed(impact_round)(vcf_fpath) for vcf_fpath in vcf_fpaths)
        info()

        # if 'polymorphic_variants' in self.filt_cnf:
        #     polymorphic_variants_path = self.filt_cnf['polymorphic_variants']
        #     if verify_file(polymorphic_variants_path):
        #         info('Filtering out polymorphic variants')
        #         with open(polymorphic_variants_path) as f:
        #

        return vcf_fpaths


def proc_line_remove_prev_filter(rec, self_):
    rec.FILTER = ['PASS']
    return rec


# Counting samples, variants and AF_by_vark
def proc_line_1st_round(rec, self_, varks, control_vars):
    # Strict filter of DP, QUAL, PMEAN
    [f.apply(rec) for f in self_.round1_filters]
    if rec.is_rejected():
        return None

    # For those who passed, collect controls, samples and af_by_varid
    sample = rec.sample_field()
    try:
        sample = sample[0]
    except:
        pass

    vark = rec.var_id()

    if sample and self_.control and sample == self_.control:
        all_passed = all(f.apply(rec, only_check=True) for f in self_.round2_filters)

        if all_passed or rec.cls() == 'Novel':
            # So that any novel variants showed up in control won't be filtered:
            control_vars.add(vark)

    if sample:
        # Undetermined won't count toward samples
        if 'undetermined' not in sample.lower() or Filtering.filt_cnf['count_undetermined']:
            if vark not in varks:
                varks[vark] = VarkInfo(vark)
            varks[vark].afs.append(rec.get_val('AF', .0))
    return rec


# Based on counted samples, variants and AF_by_vark:
#   var_n    = number of variants for vark       must be >= [sample_cnt]
#   fraction = var_n / number of total samples   must be > [fraction]
#   avg_af   = avg AF for this vark              must be < [freq]
def proc_line_2nd_round(rec, self_):
    sample = rec.sample_field()
    if sample:
        if not self_.undet_sample_filter.apply(rec):
            err('Undetermined sample for rec ' + rec.var_id() + ', sample ' + str(sample))
            return None

        vark_info = self_.varks[rec.var_id()]
        var_n = vark_info.var_n()
        frac = vark_info.frac()
        avg_af = vark_info.avg_af()

        if not self_.multi_filter.apply(rec, var_n=var_n, frac=frac, avg_af=avg_af):
            info('Multi filter: vark = ' + rec.var_id() + ', var_n = ' + str(vark_info.var_n()) + ', n_sample = ' +
                 str(len(self_.sample_names)) + ', avg_af = ' + str(vark_info.avg_af()))

        if not self_.dup_filter.apply(rec):
            info('Dup filter: vark = ' + rec.var_id())

        if rec.af() is not None:
            if not self_.max_rate_filter.apply(rec, frac=frac):
                info('Max rate filter: vark = ' + rec.var_id() + ', frac = ' + str(frac) + ', af = ' + str(rec.af()))

        self_.control_filter.apply(rec)

        gmaf = rec.get_val('GMAF')
        req_maf = Filtering.filt_cnf['maf']

        # if there's MAF with frequency, it'll be considered
        # dbSNP regardless of COSMIC
        cls = 'dbSNP' if req_maf and gmaf and gmaf > req_maf else rec.cls()

        # Rescue deleterious dbSNP, such as rs80357372 (BRCA1 Q139) that is in dbSNP,
        # but not in ClnSNP or COSMIC.
        for eff in map(Effect, rec.INFO.get('EFF', [])):
            if eff.efftype in ['STOP_GAINED', 'FRAME_SHIFT'] and cls == 'dbSNP':
                if eff.pos / int(eff.aal) < 0.95:
                    cls = 'dbSNP_del'

        self_.bias_filter.apply(rec, cls=cls)
        self_.nonclnsnp_filter.apply(rec, cls=cls)

    return rec


def proc_line_impact(rec, self_):
    self_.impact_filter.apply(rec)
    return rec


def proc_line_polymorphic(rec, self_):
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
    info('Running for ' + caller.name)

    anno_vcf_by_sample = caller.get_anno_vcf_by_samples()
    anno_vcf_fpaths = anno_vcf_by_sample.values()

    f = Filtering(cnf, bcbio_structure, caller)

    filt_anno_vcf_fpaths = f.run_filtering(anno_vcf_fpaths)

    samples = anno_vcf_by_sample.keys()

    global cnfs_for_sample_names
    info('*' * 70)
    info('Processed samples (' + str(len(samples)) + '):')
    for sample in samples:
        info(sample.name)
        cnf_copy = cnf.copy()
        cnf_copy['name'] = sample.name
        cnfs_for_sample_names[sample.name] = cnf_copy

    results = [r for r in Parallel(n_jobs=len(samples)) \
        (delayed(postprocess_vcf)
         (sample, anno_vcf_fpath, work_filt_vcf_fpath)
              for sample, anno_vcf_fpath, work_filt_vcf_fpath in
              zip(samples, anno_vcf_fpaths, filt_anno_vcf_fpaths)
         ) if r is not None and None not in r]
    info('Results: ' + str(len(results)))
    info('*' * 70)

    for sample, [vcf, tsv, maf] in zip(samples, results):
        info('Sample ' + sample.name + ': ' + vcf + ', ' + tsv + ', ' + maf)
        sample.filtered_vcf_by_callername[caller.name] = vcf
        sample.filtered_tsv_by_callername[caller.name] = tsv
        sample.filtered_maf_by_callername[caller.name] = maf

    comb_maf_fpath = join(bcbio_structure.var_dirpath, caller.name + '.maf')
    caller.combined_filt_maf_fpath = combine_mafs(cnf, caller.get_filtered_mafs(), comb_maf_fpath)
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
        info('Error: for ' + sample.name + ': cnf is None')
        return None, None, None
    work_filt_vcf_fpath = leave_first_sample(cnf, work_filt_vcf_fpath)

    final_vcf_fpath = add_suffix(original_anno_vcf_fpath, 'filt').replace('varAnnotate', 'varFilter')

    safe_mkdir(dirname(final_vcf_fpath))

    file_basepath = splitext(final_vcf_fpath)[0]
    final_vcf_fpath = file_basepath + '.vcf'
    final_tsv_fpath = file_basepath + '.tsv'
    final_maf_fpath = file_basepath + '.maf'
    final_clean_vcf_fpath = file_basepath + '.passed.vcf'  # for futrher processing

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
    shutil.move(clean_filtered_vcf_fpath, final_clean_vcf_fpath)
    # os.symlink(final_clean_vcf_fpath, clean_filtered_vcf_fpath)
    igvtools_index(cnf, final_clean_vcf_fpath)

    # Converting to TSV
    if work_filt_vcf_fpath and 'tsv_fields' in cnf:
        tsv_fpath = make_tsv(cnf, work_filt_vcf_fpath)

        if isfile(final_tsv_fpath): os.remove(final_tsv_fpath)
        shutil.move(tsv_fpath, final_tsv_fpath)
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
            normal_sample_name=sample.normal_match.name if sample.normal_match else None)
        if isfile(final_maf_fpath): os.remove(final_maf_fpath)
        shutil.move(maf_fpath, final_maf_fpath)
        info('-' * 70)
        info()
    else:
        final_maf_fpath = None

    return [final_vcf_fpath, final_tsv_fpath, final_maf_fpath]
