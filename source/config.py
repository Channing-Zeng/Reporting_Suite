import sys
from os import getcwd
from os.path import abspath, expanduser, join, dirname, pardir
from source.file_utils import verify_file, verify_module

from source.logger import info, err, critical

if verify_module('yaml'):
    from yaml import load as load_yaml
    try:
        from yaml import CDumper as Dumper, CLoader as Loader
    except ImportError:
        from yaml import Dumper, Loader
else:
    critical('Cannot import module yaml.')

cur_dirpath = dirname(abspath(__file__))


# noinspection PyClassHasNoInit
class Defaults():
    genome = 'hg19'

    output_dir = getcwd()
    tmp_base_dir = getcwd()

    verbose = True
    threads = 4
    overwrite = None
    reuse_intermediate = True
    keep_intermediate = True

    sys_cnf = join(cur_dirpath, pardir, 'system_info_Waltham.yaml')
    run_cnf = join(cur_dirpath, pardir, 'run_info.yaml')

    bcbio_final_dir = getcwd()
    steps = [
        'VarAnnotate',
        'VarQC',
        'VarFilter',
        'TargetCov',
        'NGScat']
    qualimap = False  # 'True' adds 'QualiMap' to steps
    qsub_runner = 'runner_Waltham.sh'

    coverage_reports = dict(
        depth_thresholds=[1, 5, 10, 25, 50, 100, 500, 1000, 5000, 10000, 50000],
        padding=250,
        report_types='summary,genes'
    )

    ngscat = dict(
        saturation='n',
        depthlist='auto',
        availablefeatures=[
            'percbases', 'saturation', 'specificity', 'coveragefreq',
            'coveragedistr', 'coveragestd', 'gcbias', 'coveragecorr']
    )

    variant_filtering = dict(
        impact='MODERATE|HIGH',
        vardict_mode=False,

        filt_p_mean=0,
        filt_q_mean=0,
        filt_depth=0,

        min_p_mean=5,
        min_q_mean=25,
        min_freq=0.05,
        mean_mq=20,
        mean_vd=2,
        signal_noise=4,

        fraction=0.4,
        max_ratio=1,
        sample_cnt=10,
        freq=0.15,
        count_undetermined=True,
        bias=False,
        maf=0.0025,
    )

    quality_control = dict(
        variant_distribution_scale=1000,
        databases=['dbsnp'],
        novelty=['all', 'known', 'novel'],
        metrics=['nEvalVariants', 'nSNPs', 'nInsertions', 'nDeletions',
           'nVariantsAtComp', 'compRate', 'nConcordant', 'concordantRate',
           'variantRate', 'variantRatePerBp', 'hetHomRatio', 'tiTvRatio'],
    )

    clinical_reporting = False


class Config(object):
    def __init__(self, cmd_line_opts, sys_cnf=None, run_cnf=None, **kwargs):
        object.__setattr__(self, '__d', dict())

        self.level = 0

        if sys_cnf and run_cnf:

            self.genome = None

            self.tmp_base_dir = None
            self.tmp_dir = None
            self.work_dir = None
            self.output_dir = None
            self.log = None
            self.threads = None

            self.overwrite = None
            self.reuse_intermediate = None
            self.keep_intermediate = None

            sys_cnf_fpath, run_cnf_fpath = _check_paths(sys_cnf, run_cnf)
            loaded_dict = _load(sys_cnf_fpath, run_cnf_fpath)
            for k, v in loaded_dict.items():
                self[k] = v

            for k, v in cmd_line_opts.items():
                if v is not None:
                    self[k] = v

            self.sys_cnf = sys_cnf_fpath
            self.run_cnf = run_cnf_fpath

            if self.overwrite is not None:
                self.reuse_intermediate = not self.overwrite
            else:
                self.overwrite = not self.reuse_intermediate

            if not self.tmp_base_dir:
                self.tmp_base_dir = self.work_dir
        else:
            for k, v in cmd_line_opts.items():
                self[k] = v


    def get(self, key, d=None):
        assert 'Please, use [] or . instead'
        return self.__d.get(key)

    def __contains__(self, key):
        return self[key] is not None

    def __getattribute__(self, key):
        d = object.__getattribute__(self, '__d')
        if key == '_Config__d':
            return d
        if key == 'get':
            return self.__d.get
        if key == 'copy':
            return lambda: Config(d)
        if key == 'keys':
            return lambda: [k for k in self.__d.keys() if k != 'level']
        if key == 'values':
            return lambda: [v for k, v in self.__d.items() if k != 'level']
        if key == 'items':
            return lambda: [(k, v) for k, v in self.__d.items() if k != 'level']
        else:
            return d.get(key)

    def __setattr__(self, key, value):
        if isinstance(value, dict) and not isinstance(value, Config):
            value = Config(value)
            object.__setattr__(self, 'level', self.level + 1)
        self.__d[key] = value

    def __delattr__(self, key):
        del self.__d[key]

    def __getitem__(self, key):
        return self.__d.get(key)

    def __setitem__(self, key, value):
        if isinstance(value, dict) and not isinstance(value, Config):
            value = Config(value)
            value.level = self.level + 1
        self.__d[key] = value

    def __delitem__(self, key):
        del self.__d[key]

    def __repr__(self):
        s = ''
        for k, v in self.__d.items():
            if k == 'level':
                continue
            s += '  ' * self.level
            s += str(k) + ': '
            if isinstance(v, Config):
                s += '\n'
            s += repr(v)
            s += '\n'
        return s

    def __len__(self):
        return self.__d.__len__()


def load_yaml_config(fpath):
    if not verify_file(fpath):
        sys.exit(1)
    dic = load_yaml(open(fpath), Loader=Loader)
    return dic


def _load(sys_cnf_fpath, run_cnf_fpath):
    sys_dict = load_yaml(open(sys_cnf_fpath), Loader=Loader)
    run_dict = load_yaml(open(run_cnf_fpath), Loader=Loader)

    loaded_dict = dict(run_dict.items() + sys_dict.items())
    loaded_dict = _fill_dict_from_defaults(loaded_dict, Defaults.__dict__)

    info('Loaded system config ' + sys_cnf_fpath)
    info('Loaded run config ' + run_cnf_fpath)
    info()
    return loaded_dict


def _fill_dict_from_defaults(cur_cnf, defaults_dict):
    for key in defaults_dict:
        if key in cur_cnf:
            if isinstance(cur_cnf[key], dict) and isinstance(defaults_dict[key], dict):
                _fill_dict_from_defaults(cur_cnf[key], defaults_dict[key])
        else:
            cur_cnf[key] = defaults_dict[key]
    return cur_cnf


def _check_paths(sys_cnf, run_cnf):
    to_exit = False

    info('Using ' + sys_cnf + ' as a system configuration file.')
    info('Using ' + run_cnf + ' as a run configuration file.')
    info()

    for fn in [sys_cnf, run_cnf]:
        if not verify_file(fn, 'Config'):
            to_exit = True
    if to_exit:
        sys.exit(1)

    sys_cnf_path = abspath(expanduser(sys_cnf))
    run_cnf_path = abspath(expanduser(run_cnf))

    for fn in [sys_cnf_path, run_cnf_path]:
        if not fn.endswith('.yaml'):
            err(fn + ' does not end with .yaml, maybe incorrect parameter?')
            err()
            to_exit = True
    if to_exit:
        sys.exit(1)

    return sys_cnf_path, run_cnf_path


def join_parent_conf(child_conf, parent_conf):
    bc = parent_conf.copy()
    bc.update(child_conf)
    child_conf.update(bc)
    return child_conf


