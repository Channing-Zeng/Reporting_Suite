import sys
from source.variants.vcf_processing import iterate_vcf

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

import operator
from collections import defaultdict

from source.variants.Effect import Effect
from source.logger import step_greetings, info, critical, err
from source.file_utils import convert_file
from source.utils import mean


class Filter:
    filt_cnf = None

    def __init__(self, word, check=None, required=True):
        self.check = check
        self.required = required
        self.word = word.upper()

        self.num_passed = 0
        self.num_rejected = 0

    @staticmethod
    def __remove_pass(rec):
        if rec.FILTER == ['PASS']:
            rec.FILTER = None

    def apply(self, rec, only_check=False):
        assert self.check, 'check function must be specified for filter before applying'

        if only_check:
            return self.check(rec)

        if self.check(rec):
            self.num_passed += 1
        else:
            self.num_rejected += 1

            if self.word not in rec.FILTER:
                self.__remove_pass(rec)
                rec.add_filter(self.word)


class CnfFilter(Filter):
    def __init__(self, key, *args, **kvargs):
        self.key = key
        Filter.__init__(self, key.upper(), *args, **kvargs)

    def apply(self, rec, only_check=False):
        if Filter.filt_cnf.get(self.key) is None:
            return True

        return Filter.apply(self, rec, only_check)


class InfoFilter(CnfFilter):
    def __init__(self, cnf_key, info_key, op=operator.ge, *args, **kwargs):
        def check(rec):
            if info_key not in rec.INFO:
                if self.required:
                    critical('Error: no field ' + info_key + ' in variant at line ' +
                             str(rec.line_num + 1) + ' (' + rec.CHROM + ':' + str(rec.POS) +
                             ') - required to test ' + cnf_key)
                else:
                    return True  # PASS

            return op(Filter.filt_cnf[cnf_key], rec.INFO[info_key])

        CnfFilter.__init__(self, cnf_key, check, *args, **kwargs)


class EffectFilter(CnfFilter):
    def __init__(self, cnf_key, *args, **kwargs):
        def check(rec):
            if 'EFF' not in rec.INFO:
                err('Warning: EFF field is missing for variant at line ' + str(rec.line_num))
                return False

            impact_filter_line = Filter.filt_cnf[cnf_key]
            if impact_filter_line:
                important_impacts = [s.upper() for s in impact_filter_line.split('|')]

                return any(eff.impact.upper() in important_impacts
                       for eff in map(Effect, rec.INFO['EFF']))

        CnfFilter.__init__(self, cnf_key, check, *args, **kwargs)


class Filtering:
    def __init__(self, cnf, filt_cnf, vcf_fpath):
        self.cnf = cnf
        self.filt_cnf = filt_cnf
        self.vcf_fpath = vcf_fpath
        self.vardict_mode = False
        Filter.filt_cnf = self.filt_cnf

        self.control_vars = set()
        self.samples = {''}
        self.af_by_varid = defaultdict(list)

        self.round1_filters = [InfoFilter('filt_depth', 'DP', required=False)]
        self.round1_filters.append(InfoFilter('filt_q_mean', 'QUAL', required=False))
        self.round1_filters.append(InfoFilter('filt_p_mean', 'PMEAN', required=False))
        # self.round1_filters.append(Filter('min_q_mean', lambda rec: rec.QUAL >= filt_cnf['filt_q_mean']))

        self.control = self.filt_cnf.get('control')

        self.impact_filter = EffectFilter('impact')

        self.round2_filters = [
            InfoFilter('min_p_mean', 'PMEAN', required=False),
            InfoFilter('min_q_mean', 'QUAL', required=False),
            InfoFilter('min_freq', 'AF', required=False),
            InfoFilter('mean_mq', 'MQ', required=False),
            InfoFilter('signal_noise', 'SN', required=False),
            InfoFilter('mean_vd', 'VD', required=False)]

        self.undet_sample_filter = Filter('UNDET_SAMPLE', lambda rec: rec.var_id() in self.af_by_varid)
        self.multi_filter = Filter('MULTI')
        self.dup_filter = Filter('DUP')
        self.max_rate_filter = CnfFilter('max_ratio')
        self.control_filter = CnfFilter('control', lambda rec: filt_cnf['control'] and rec.var_id() in self.control_vars)
        self.bias_filter = CnfFilter('bias')
        self.nonclnsnp_filter = Filter('NonClnSNP')

    def run_filtering(self):
        step_greetings('Filtering')

        vcf_fpath = self.vcf_fpath

        info('Removing previous FILTER values')
        vcf_fpath = iterate_vcf(self.cnf, vcf_fpath, self.get_proc_line_remove_prev_filter(), suffix='rm_prev')
        info('Saved to ' + vcf_fpath)
        info()

        info('First round')
        vcf_fpath = iterate_vcf(self.cnf, vcf_fpath, self.get_proc_line_1st_round(), suffix='r1')
        info('Saved to ' + vcf_fpath)
        info()

        info('Second round')
        vcf_fpath = iterate_vcf(self.cnf, vcf_fpath, self.get_proc_line_2nd_round(), suffix='r2')
        info('Saved to ' + vcf_fpath)
        info()

        info('Third round')
        vcf_fpath = iterate_vcf(self.cnf, vcf_fpath, self.get_proc_line_3rd_round(), suffix='r3')
        info('Saved to ' + vcf_fpath)
        info()

        return vcf_fpath

    def get_proc_line_remove_prev_filter(self):
        def __f(rec):
            rec.FILTER = ['PASS']
            return rec
        return __f

    # DP, QUAL, PMEAN
    def get_proc_line_1st_round(self):
        def __f(rec):
            [f.apply(rec) for f in self.round1_filters]
            return rec
        return __f

    # counting samples, variants and AF_by_vark
    def get_proc_line_2nd_round(self):
        def __f(rec):
            sample = rec.sample()
            try:
                sample = sample[0]
            except:
                pass

            if sample and self.control and sample == self.control:
                all_passed = all([f.apply(rec, only_check=True) for f in self.round2_filters])

                if all_passed or rec.cls() == 'Novel':
                # So that any novel variants showed up in control won't be filtered:
                    self.control_vars.add(rec.var_id())

            if sample:
                if 'undetermined' not in sample.lower() or self.filt_cnf['count_undetermined']:
                # Undetermined won't count toward samples
                    self.samples.add(sample)
                    self.af_by_varid[rec.var_id()].append(rec.INFO.get('AF', .0))
            return rec
        return __f

    # based on counted samples, variants and AF_by_vark:
    #   var_n    = number of variants for vark       must be >= [sample_cnt]
    #   fraction = var_n / number of total samples   must be > [fraction]
    #   avg_af   = avg AF for this vark              must be < [freq]
    def get_proc_line_3rd_round(self):
        def __f(rec):
            self.impact_filter.apply(rec)

            sample = rec.sample()
            if sample:
                self.undet_sample_filter.apply(rec)
                if rec.is_rejected():
                    return rec

                var_n = len(self.af_by_varid[rec.var_id()])
                fraction = float(var_n) / len(self.samples)
                avg_af = mean(self.af_by_varid[rec.var_id()])
                self.multi_filter.check = lambda _: not (  # novel and present in [max_ratio] samples
                    var_n >= self.filt_cnf['sample_cnt'] and
                    fraction > self.filt_cnf['fraction'] and
                    avg_af < self.filt_cnf['freq'] and
                    rec.ID == '.')  # TODO: check if "." converted to None in the vcf lib
                self.multi_filter.apply(rec)

                pstd = rec.INFO.get('PSTD')
                bias = rec.INFO.get('BIAS')
                # all variants from one position in reads
                if pstd is not None and bias is not None:
                    self.dup_filter.check = lambda: pstd != 0 or bias[-1] in ['0', '1']
                    self.dup_filter.apply(rec)

                max_ratio = self.filt_cnf.get('max_ratio')
                af = rec.INFO.get('AF')
                if af is not None:
                    af = float(af)
                    self.max_rate_filter.check = lambda _: fraction < max_ratio or af < 0.3
                    self.max_rate_filter.apply(rec)

                gmaf = rec.INFO.get('GMAF')
                req_maf = self.filt_cnf['maf']
                # if there's MAF with frequency, it'll be considered
                # dbSNP regardless of COSMIC
                if gmaf:
                    cls = 'dbSNP' if req_maf and gmaf > req_maf else rec.cls()

                    self.control_filter.apply(rec)

                    # Rescue deleterious dbSNP, such as rs80357372 (BRCA1 Q139) that is in dbSNP,
                    # but not in ClnSNP or COSMIC.
                    for eff in map(Effect, rec.FILT.get('EFF', [])):
                        if eff.efftype in ['STOP_GAINED', 'FRAME_SHIFT'] and cls == 'dbSNP':
                            if eff.pos / int(eff.aal) < 0.95:
                                cls = 'dbSNP_del'

                    self.bias_filter.check = lambda _: not (  # Filter novel variants with strand bias.
                        self.filt_cnf['bias'] is True and
                        cls in ['Novel', 'dbSNP'] and
                        bias and bias in ['2;1', '2;0'] and af < 0.3)
                    self.bias_filter.apply(rec)

                    self.nonclnsnp_filter.check = lambda _: rec.check_clnsig() != -1 or cls == 'COSMIC'
                    self.nonclnsnp_filter.apply(rec)

            return rec
        return __f


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