import sys
from os import getcwd
from os.path import abspath, expanduser, join, dirname, pardir

from source.logger import info, err, critical
from source.utils import verify_file, verify_module

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
    base_tmp_dir = getcwd()

    verbose = False
    threads = 4
    overwrite = False
    reuse_intermediate = True
    keep_intermediate = True

    sys_cnf = join(cur_dirpath, pardir, 'system_info_Waltham.yaml')
    run_cnf = join(cur_dirpath, pardir, 'run_info.yaml')

    bcbio_final_dir = getcwd()
    steps = ['IndelFilter',
             'VarAnnotate', 'VarQC', 'FilterVariants', 'TargetCoverage', 'NGSCat', 'QualiMap']
    qualimap = False

    coverage_reports = dict(
        depth_thresholds=[1, 5, 10, 25, 50, 100, 500, 1000, 5000, 10000, 50000],
        padding=250,
        report_types='summary,genes')

    ngscat = dict(
        saturation='n',
        depthlist='auto',
        availablefeatures=[
            'percbases', 'saturation', 'specificity', 'coveragefreq',
            'coveragedistr', 'coveragestd', 'gcbias', 'coveragecorr'])

    quality_control = dict(variant_distribution_scale=1000)

    clinical_reporting = False


class Config():
    def __init__(self, cmd_line_opts, sys_cnf, run_cnf):

        sys_cnf_fpath, run_cnf_fpath = _check_paths(sys_cnf, run_cnf)

        self.genome = None

        self.base_tmp_dir = None
        self.tmp_dir = None
        self.work_dir = None
        self.output_dir = None
        self.log = None
        self.threads = None

        self.overwrite = False
        self.reuse_intermediate = None
        self.keep_intermediate = None

        self.__dict__ = _load(sys_cnf_fpath, run_cnf_fpath)

        self.sys_cnf = sys_cnf_fpath
        self.run_cnf = run_cnf_fpath

        for k, v in cmd_line_opts.items():
            if k not in self.__dict__ or v is not None:
                self.__dict__[k] = v

        if self.overwrite:
            self.reuse_intermediate = False

        if not self.base_tmp_dir:
            self.base_tmp_dir = self.work_dir

    def get(self, key):
        return self.__dict__.get(key)

    def __contains__(self, key):
        return key in self.__dict__

    def __getattr__(self, key):
        return self.__dict__.get(key)

    def __setattr__(self, key, value):
        self.__dict__[key] = value

    def __delattr__(self, key):
        del self.__dict__[key]

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __delitem__(self, key):
        del self.__dict__[key]

    def __repr__(self):
        return self.__dict__.__repr__()



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

