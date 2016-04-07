import sys
from contextlib import contextmanager
from os import getcwd
from os.path import abspath, expanduser, join, dirname, pardir
from traceback import format_exc

import bcbio_postproc
from source import logger
from source.file_utils import verify_file, verify_module, adjust_path

from source.logger import info, err, critical, debug

if verify_module('yaml'):
    from yaml import load as load_yaml
    try:
        from yaml import CDumper as Dumper, CLoader as Loader
    except ImportError:
        from yaml import Dumper, Loader
else:
    critical('Error: cannot import module yaml. ')

# this_script_dirpath = dirname(abspath(__file__))
# configs_dirpath = abspath(join(abspath(this_script_dirpath), pardir, 'configs'))


configs_dirpath = abspath(join(abspath(bcbio_postproc.project_dir), 'configs'))
test_dirpath = abspath(join(abspath(bcbio_postproc.project_dir), 'test'))

defaults = dict(
    sys_cnfs = dict(
        us = join(configs_dirpath, 'system_info_Waltham.yaml'),
        uk = join(configs_dirpath, 'system_info_AP.yaml'),
        china = join(configs_dirpath, 'system_info_China.yaml'),
        cloud = join(configs_dirpath, 'system_info_cloud.yaml'),
        local = abspath(join(test_dirpath, 'system_info.yaml')),
    ),
    run_cnf_exome_seq = join(configs_dirpath, 'run_info_ExomeSeq.yaml'),
    run_cnf_wgs = join(configs_dirpath, 'run_info_WGS.yaml'),
    run_cnf_deep_seq = join(configs_dirpath, 'run_info_DeepSeq.yaml'),
    run_cnf_rnaseq = join(configs_dirpath, 'run_info_RNAseq.yaml'),

    load_mongo = False,  # 'True' adds 'LoadMongo' to steps
    qsub_runner = 'runner_Waltham.sh',

    default_min_freq = 0.05,

    threads = 1,
    bcbio_postproc_threads = 30
)
defaults['sys_cnf'] = defaults['sys_cnfs']['us']

defaults_yaml_fpath = join(configs_dirpath, 'RUNINFO_DEFAULTS.yaml')
verify_file(defaults_yaml_fpath, is_critical=True)
run_info_defaults = load_yaml(open(defaults_yaml_fpath), Loader=Loader)
for k, v in run_info_defaults.items():
    defaults[k] = v


class Config(object):
    def __init__(self, d, sys_cnf=None, run_cnf=None, bcbio_genome_build=None, **kwargs):
        object.__setattr__(self, '__d', dict())

        self.level = 0

        if sys_cnf:  # Improtant if! Do no delete because __init__ called recursivly (implicitly)
            sys_cnf_fpath, run_cnf_fpath = _check_paths(sys_cnf, run_cnf)
            loaded_dict = _load(sys_cnf_fpath, run_cnf_fpath)
            for k, v in loaded_dict.items():
                # if k == 'annotation':
                #     pass
                self[k] = v

            if d:
                for k, v in d.items():
                    if v is not None:
                        self[k] = v

            self.sys_cnf = sys_cnf_fpath
            if run_cnf_fpath:
                self.run_cnf = run_cnf_fpath
            self.tmp_base_dir = self.work_dir

            if not self.genome:
                self.genome = bcbio_genome_build
            if self.genomes and self.genome:
                build_name = self.genome
                if build_name not in self.genomes:
                    critical('Genome ' + str(build_name) + ' not in ' + sys_cnf_fpath)
                self.genome = Config(self.genomes[build_name])
                self.genome.name = build_name
        else:
            for k, v in d.items():
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
        if key == '__dict__':
            return dict([(k, v) for k, v in self.__d.items() if k != 'level'])
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

    # def __genome_resources(self):
    #     cnf = self
    #     if cnf.genomes and cnf.genome:
    #         build_name = cnf.genome
    #         cnf.genome = cnf.genomes[build_name]
    #         cnf.genome.name = build_name


class CallCnf:
    def __init__(self, dict_cnf):
        self.work_dir = dict_cnf.get('work_dir')
        self.log = dict_cnf.get('log')
        self.reuse = self.reuse_intermediate = dict_cnf.get('reuse_intermediate')
        self.keep_intermediate = dict_cnf.get('keep_intermediate')
        self.verbose = dict_cnf.get('verbose')
        self.threads = dict_cnf.get('threads')


def load_yaml_config(fpath):
    verify_file(fpath, is_critical=True)
    try:
        dic = load_yaml(open(fpath), Loader=Loader)
    except:
        err(format_exc())
        critical('Could not parse bcbio YAML ' + fpath)
    else:
        return dic


def _load(sys_cnf_fpath=None, run_cnf_fpath=None):
    sys_dict = load_yaml(open(sys_cnf_fpath), Loader=Loader)
    info('Loaded system config ' + sys_cnf_fpath)
    run_dict = dict()
    if run_cnf_fpath:
        run_dict = load_yaml(open(run_cnf_fpath), Loader=Loader)
        info('Loaded run config ' + run_cnf_fpath)

    loaded_dict = dict(run_dict.items() + sys_dict.items())
    loaded_dict = fill_dict_from_defaults(loaded_dict, defaults)
    info()
    return loaded_dict


def fill_dict_from_defaults(cur_cnf, defaults_dict):
    for key in defaults_dict:
        if key in cur_cnf:
            if isinstance(cur_cnf[key], dict) and isinstance(defaults_dict[key], dict):
                fill_dict_from_defaults(cur_cnf[key], defaults_dict[key])
        else:
            cur_cnf[key] = defaults_dict[key]
    return cur_cnf


def _check_paths(sys_cnf=None, run_cnf=None):
    to_exit = False

    debug('System configuration file: ' + str(sys_cnf))
    if run_cnf:
        debug('Run configuration file: ' + str(run_cnf))
    debug()

    sys_cnf = verify_file(sys_cnf, 'System config', is_critical=True)
    if run_cnf:
        run_cnf = verify_file(run_cnf, 'Run config', is_critical=True)

    errors = []
    for fn in [sys_cnf, run_cnf]:
        if fn and not fn.endswith('.yaml'):
            errors.append(fn + ' does not end with .yaml, maybe incorrect parameter?')
    if errors:
        critical(errors)

    return sys_cnf, run_cnf


def join_parent_conf(child_conf, parent_conf):
    bc = parent_conf.copy()
    bc.update(child_conf)
    child_conf.update(bc)
    return child_conf


@contextmanager
def with_cnf(cnf, **kwargs):
    prev_opts = {k: cnf[k] for k in kwargs.keys()}
    try:
        for k, v in kwargs.items():
            cnf[k] = v
        yield cnf
    finally:
        for k, v in prev_opts.items():
            cnf[k] = v
