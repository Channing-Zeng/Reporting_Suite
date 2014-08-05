from os.path import basename
import sys

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

import operator
from collections import defaultdict, namedtuple
from joblib import Parallel, delayed

# try:
#     import dill as pickle
#     dill_loaded = True
# except ImportError:
#     import pickle
#     dill_loaded = False

from source.variants.Effect import Effect
from source.logger import step_greetings, info, critical, err
from source.variants.vcf_processing import iterate_vcf, vcf_one_per_line
from source.utils import mean


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

    def __init__(self, cnf, filt_cnf, sample_names):
        Filtering.cnf = cnf
        filt_cnf = filt_cnf.__dict__
        Filter.filt_cnf = filt_cnf
        Filtering.filt_cnf = filt_cnf

        self.control_vars = set()
        self.sample_names = set(sample_names)
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

        def polymorphic_filter_check(rec):
            pass
        self.polymorphic_filter = Filter('POLYMORPHIC', polymorphic_filter_check)

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


    def run_filtering(self, vcf_fpaths):
        step_greetings('Filtering')

        info('Removing previous FILTER values')

        n_jobs = len(vcf_fpaths)
        # n_jobs = 1

        global cnf_for_samples, filtering
        filtering = self
        for vcf_fpath in vcf_fpaths:
            cnf_for_samples[basename(vcf_fpath).split('.')[0]] = Filtering.cnf.copy()

        vcf_fpaths = Parallel(n_jobs=n_jobs)(delayed(rm_prev_round)(vcf_fpath)
                                             for vcf_fpath in vcf_fpaths)
        info()

        info('First round')
        results = Parallel(n_jobs=n_jobs)(delayed(first_round)(vcf_fpath)
                                          for vcf_fpath in vcf_fpaths)

        for vcf_fpath, varks, control_vars in results:
            for vark, vark_info in varks.items():
                if vark == 'chrM:150:T:[C]':
                    pass
                if vark not in self.varks:
                    self.varks[vark] = vark_info
                else:
                    self.varks[vark].afs.extend(vark_info.afs)

            self.control_vars.update(control_vars)
        filtering = self
        info()

        info('One effect per line')
        vcf_fpaths = Parallel(n_jobs=n_jobs)(delayed(one_per_line)(vcf_fpath)
                                             for vcf_fpath in vcf_fpaths)
        info()

        info('Second round')
        vcf_fpaths = Parallel(n_jobs=n_jobs)(delayed(second_round)(vcf_fpath)
                                             for vcf_fpath in vcf_fpaths)
        info()

        info('Filtering by impact')
        vcf_fpaths = Parallel(n_jobs=n_jobs)(delayed(impact_round)(vcf_fpath)
                                             for vcf_fpath in vcf_fpaths)
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
        return rec

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
            return rec

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