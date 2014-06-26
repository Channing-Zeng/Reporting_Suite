import sys

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

import operator
from collections import defaultdict

from ext_modules import vcf
from ext_modules.vcf.model import _Record

from source.variants import Effect
from source.logger import step_greetings, info, critical
from source.calling_process import call
from source.file_utils import intermediate_fname, convert_file
from source.tools_from_cnf import get_java_tool_cmdline
from source.utils import mean


class Filter:
    filt_cnf = None

    def __init__(self, word, check, required=True):
        self.check = check
        self.required = required
        self.word = word

        self.num_passed = 0
        self.num_rejected = 0

    @staticmethod
    def __remove_pass(rec):
        if rec.FILTER in ['PASS', '.']:
            rec.FILTER = None

    def apply(self, rec):
        if self.check(rec):
            self.num_passed += 1
        else:
            self.num_rejected += 1

            if self.word() not in rec.FILTER:
                self.__remove_pass(rec)
                rec.add_filter(self.word)


class CnfFilter(Filter):
    def __init__(self, key, *args, **kvargs):
        self.key = key
        Filter.__init__(self, key.upper(), *args, **kvargs)

    def apply(self, rec):
        if Filter.filt_cnf.get(self.key) is None:
            return

        Filter.apply(self, rec)


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
                critical('Error: in variant line ' + str(rec.line_num + 1) +
                         ' (' + rec.CHROM + ':' + str(rec.POS) +
                         '), EFF field missing in INFO column')

            req = Filter.filt_cnf[cnf_key]
            if req:
                req_values = [s.upper() for s in req.split('|')]

                for eff in map(Effect, rec.INFO['EFF']):
                    if eff.impact.upper() not in req_values:
                        return False

        CnfFilter.__init__(self, cnf_key, check, *args, **kwargs)


def _proc_vcf(inp_f, out_f, proc_line_fun):
    reader = vcf.Reader(inp_f)
    writer = vcf.Writer(out_f, reader)

    for i, rec in enumerate(reader):
        rec = proc_line_fun(Record(rec, i))
        if rec:
            writer.write_record(rec)


class Record(_Record):
    # noinspection PyMissingConstructor
    def __init__(self, _record, line_num):
        self.__dict__.update(_record.__dict__)
        self.line_num = line_num

    def cls(self):
        cls = 'Novel'
        if 'COSM' in self.ID:
            cls = 'COSMIC'
        elif self.ID.startswith('rs'):
            if self.check_clnsig:
                cls = 'ClnSNP'
            else:
                cls = 'dbSNP'
        return cls

    def is_rejected(self):
        if self.FILTER:
            assert '.' not in self.FILTER
        return self.FILTER and 'PASS' not in self.FILTER

    def check_clnsig(self):
        if not self.INFO.get('CLNSIG'):
            return 0

        for c in self.INFO.get('CLNSIG'):
            if 3 < c < 7 or c == 255:
                return 1

        return -1

    def sample(self):
        return self.INFO.get('SAMPLE')

    def var_id(self):
        return ':'.join(map(str, [self.CHROM, self.POS, self.REF, self.ALT]))


def run_filtering(cnf, filt_cnf, vcf_fpath):
    control_vars = set()
    samples = {''}
    af_by_varid = defaultdict(list)

    Filter.filt_cnf = filt_cnf
    vardict_mode = cnf['vardict_mode']

    first_filters = [InfoFilter('filt_depth', 'DP')]
    if vardict_mode:
        first_filters.append(InfoFilter('filt_q_mean', 'QUAL'))
        first_filters.append(InfoFilter('filt_p_mean', 'PMEAN'))
    else:
        first_filters.append(Filter('min_q_mean', lambda rec, thres: rec.QUAL >= thres))

    def proc_line_1st_round(rec):
        [f.apply(rec) for f in first_filters]
        if rec.is_rejected():
            return rec

    round2_filters = [
        InfoFilter('min_q_mean', 'QUAL'),
        InfoFilter('min_p_mean', 'PMEAN'),
        InfoFilter('min_freq', 'AF'),
        InfoFilter('min_mq', 'MQ'),
        InfoFilter('signal_noise', 'SN'),
        InfoFilter('mean_vd', 'VD', required=False)]

    control = filt_cnf.get('control')

    impact_filter = EffectFilter('impact')

    def proc_line_2nd_round(rec):
        if vardict_mode:
            sample = rec.sample()

            if sample and control and sample == control:
                [f.apply(rec) for f in round2_filters]

                if not rec.is_rejected() or rec.cls() == 'Novel':
                # So that any novel variants showed up in control won't be filtered:
                    control_vars.add(rec)

            if sample:
                if 'undetermined' not in sample.lower() or filt_cnf['count_undetermined']:
                # Undetermined won't count toward samples
                    samples.add(sample)
                    af_by_varid[rec.var_id()].append(rec.INFO.get('AF', .0))

        impact_filter.apply(rec)

        return rec

    undet_sample_filter = Filter('UNDET_SAMPLE', lambda rec: rec.var_id() in af_by_varid)

    multi_filter = Filter('MULTI', None)

    def proc_line_3rd_round(rec):
        if vardict_mode:
            undet_sample_filter.apply(rec)
            if rec.is_rejected():
                return rec

        var_n = len(af_by_varid[rec.var_id()])
        fraction = float(var_n) / len(samples)
        avg_af = mean(af_by_varid[rec.var_id()])

        # novel and present in [max_ratio] samples
        def milti_filter_check(rec):
            return not (fraction > filt_cnf['fraction'] and
                        var_n >= filt_cnf['sample_cnt'] and
                        avg_af < filt_cnf['freq'] and
                        rec.ID == '.')  # TODO: check if "." converted to None in the vcf lib
        multi_filter.check = milti_filter_check

        multi_filter.apply(rec)

        if 'PSTD' in d and d['PSTD'] == 0 and \
           'BIAS' in d and not (d['BIAS'].endswith('0') or d['BIAS'].endswith('1')):
            reject_values.append('DUP')

        if fraction >= filt_cnf['max_ratio'] and 'AF' in d and float(d['AF']) < 0.3:
            reject_values.append('MAXRATE')

        for real_key, test_key in [
             ('QUAL', 'min_q_mean'),
             # ('PMEAN', 'min_p_mean'),
             ('MQ', 'min_mq'),
             # ('SN', 'signal_noise'),
             ('AF', 'min_freq'),
             # ('VD', 'mean_vd')
        ]:
            if less(real_key, test_key):
                reject_values.append(test_key.upper())

        if 'control' in filt_cnf and vark in control_dict:
            reject_values.append('CNTL')

        cls = get_class(d, rec.ID)
        if greater('GMAF', 'maf'):
            # if there's MAF with frequency, it'll be considered
            # dbSNP regardless of COSMIC
            cls = 'dbSNP'

        if 'EFF' not in d:
            critical('Error: in variant line ' + str(line_num + 1) + ', EFF field missing in INFO column')

        # Rescue deleterious dbSNP, such as rs80357372 (BRCA1 Q139* that is
        # in dbSNP, but not in ClnSNP or COSMIC
        for eff in map(Effect, d['EFF']):
            if eff.efftype in ['STOP_GAINED', 'FRAME_SHIFT'] and cls == 'dbSNP':
                if eff.pos / int(eff.aal) < 0.95:
                    cls = 'dbSNP_del'

        if ('bias' in filt_cnf and filt_cnf['bias'] and (cls == 'Novel' or cls == 'dbSNP') and
            'BIAS' in d and (d['BIAS'] == '2;1' or d['BIAS'] == '2;0') and
            'AF' in d and float(d['AF']) < 0.3):
            reject_values.append('BIAS')

        if check_clnsig(d) == -1 and cls != 'COSMIC':
            reject_values.append('NonClnSNP')

        for val in reject_values:
            rec = _add_reject(rec, val)

        return rec

    step_greetings('Filtering')

    proc_vcf_1st_round = lambda inp_f, out_f: _proc_vcf(inp_f, out_f, proc_line_1st_round)

    proc_vcf_2nd_round = lambda inp_f, out_f: _proc_vcf(inp_f, out_f, proc_line_2nd_round)

    proc_vcf_3rd_round = lambda inp_f, out_f: _proc_vcf(inp_f, out_f, proc_line_3rd_round)

    info('First round')
    vcf_fpath = convert_file(cnf, vcf_fpath, proc_vcf_1st_round, suffix='zh1')

    info('Second round')
    vcf_fpath = convert_file(cnf, vcf_fpath, proc_vcf_2nd_round, suffix='zh2')

    info('Third round')
    vcf_fpath = convert_file(cnf, vcf_fpath, proc_vcf_3rd_round, suffix='zh3')

    return vcf_fpath


def remove_previous_filter_rejects(cnf, vcf_fpath):
    step_greetings('Removing previous FILTER values.')

    def proc_line(rec, i, vark, d):
        rec.FILTER = 'PASS'
        return rec

    proc_vcf = lambda inp_f, out_f: _proc_vcf(inp_f, out_f, proc_line)
    vcf_fpath = convert_file(cnf, vcf_fpath, proc_vcf, suffix='rpp')

    info('Done.')

    return vcf_fpath


def run_snpsift(cnf, vcf_cnf, vcf_fpath):
    expression = vcf_cnf.get('expression')
    if not expression:
        return vcf_fpath

    step_greetings('Running SnpSift filter.')

    executable = get_java_tool_cmdline(cnf, 'snpsift')
    cmdline = '{executable} filter -a EXPR -n -p -f ' \
              '{vcf_fpath} "{expression}"'.format(**locals())
    filtered_fpath = intermediate_fname(cnf, vcf_fpath, 'snpsift')
    call(cnf, cmdline, filtered_fpath)

    info('Done.')

    return filtered_fpath