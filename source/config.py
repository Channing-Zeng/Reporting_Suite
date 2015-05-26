import sys
from os import getcwd
from os.path import abspath, expanduser, join, dirname, pardir
from source.file_utils import verify_file, verify_module, adjust_path

from source.logger import info, err, critical

if verify_module('yaml'):
    from yaml import load as load_yaml
    try:
        from yaml import CDumper as Dumper, CLoader as Loader
    except ImportError:
        from yaml import Dumper, Loader
else:
    critical('Error: cannot import module yaml. ')

this_script_dirpath = dirname(abspath(__file__))


defaults = dict(
    sys_cnfs = dict(
        us=abspath(join(this_script_dirpath, pardir, 'system_info_Waltham.yaml')),
        uk=abspath(join(this_script_dirpath, pardir, 'system_info_AP.yaml')),
        china=abspath(join(this_script_dirpath, pardir, 'system_info_China.yaml')),
        cloud=abspath(join(this_script_dirpath, pardir, 'system_info_cloud.yaml')),
        local=abspath(join(this_script_dirpath, pardir, 'test', 'system_info.yaml')),
    ),
    run_cnf = abspath(join(this_script_dirpath, pardir, 'run_info.yaml')),
    run_cnf_deep_seq = abspath(join(this_script_dirpath, pardir, 'run_info_DeepSeq.yaml')),

    load_mongo = False,  # 'True' adds 'LoadMongo' to steps
    qsub_runner = 'runner_Waltham.sh',

    default_min_freq = 0.05,

    genome='hg19',
)
defaults['sys_cnf'] = defaults['sys_cnfs']['us']


defaults_yaml_fpath = abspath(join(this_script_dirpath, pardir, 'RUNINFO_DEFAULTS.yaml'))
verify_file(defaults_yaml_fpath, is_critical=True)
run_info_defaults = load_yaml(open(defaults_yaml_fpath), Loader=Loader)
for k, v in run_info_defaults.items():
    defaults[k] = v


class Config(object):
    def __init__(self, cmd_line_opts, sys_cnf=None, run_cnf=None, **kwargs):
        object.__setattr__(self, '__d', dict())

        self.level = 0

        if sys_cnf and run_cnf:

            sys_cnf_fpath, run_cnf_fpath = _check_paths(sys_cnf, run_cnf)
            loaded_dict = _load(sys_cnf_fpath, run_cnf_fpath)
            for k, v in loaded_dict.items():
                self[k] = v

            for k, v in cmd_line_opts.items():
                if v is not None:
                    self[k] = v

            self.sys_cnf = sys_cnf_fpath
            self.run_cnf = run_cnf_fpath
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


class CallCnf:
    def __init__(self, dict_cnf):
        self.work_dir = dict_cnf['work_dir']
        self.log = dict_cnf.get('log')
        self.reuse_intermediate = dict_cnf['reuse_intermediate']
        self.keep_intermediate = dict_cnf['keep_intermediate']
        self.verbose = dict_cnf.get('verbose')


def load_yaml_config(fpath):
    verify_file(fpath, is_critical=True)
    dic = load_yaml(open(fpath), Loader=Loader)
    return dic


def _load(sys_cnf_fpath, run_cnf_fpath):
    sys_dict = load_yaml(open(sys_cnf_fpath), Loader=Loader)
    run_dict = load_yaml(open(run_cnf_fpath), Loader=Loader)

    loaded_dict = dict(run_dict.items() + sys_dict.items())
    loaded_dict = fill_dict_from_defaults(loaded_dict, defaults)

    info('Loaded system config ' + sys_cnf_fpath)
    info('Loaded run config ' + run_cnf_fpath)
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


def _check_paths(sys_cnf, run_cnf):
    to_exit = False

    info('System configuration file: ' + sys_cnf)
    info('Run configuration file:    ' + run_cnf)
    info()

    verify_file(sys_cnf, 'System config', is_critical=True)
    verify_file(run_cnf, 'Run config', is_critical=True)

    sys_cnf_path = adjust_path(sys_cnf)
    run_cnf_path = adjust_path(run_cnf)

    errors = []
    for fn in [sys_cnf_path, run_cnf_path]:
        if not fn.endswith('.yaml'):
            errors.append(fn + ' does not end with .yaml, maybe incorrect parameter?')
    if errors:
        critical(errors)

    return sys_cnf_path, run_cnf_path


def join_parent_conf(child_conf, parent_conf):
    bc = parent_conf.copy()
    bc.update(child_conf)
    child_conf.update(bc)
    return child_conf


