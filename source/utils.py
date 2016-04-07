import hashlib
from os import environ
import socket
import re
from collections import OrderedDict
from os.path import join, basename

from source.logger import info, critical, err
from source.file_utils import file_exists, verify_file, file_transaction, adjust_path


class OrderedDefaultDict(OrderedDict):
    def __init__(self, *args, **kwargs):
        if not args:
            self.default_factory = None
        else:
            if not (args[0] is None or callable(args[0])):
                raise TypeError('first argument must be callable or None')
            self.default_factory = args[0]
            args = args[1:]
        super(OrderedDefaultDict, self).__init__(*args, **kwargs)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = default = self.default_factory()
        return default

    def __reduce__(self):  # optional, for pickle support
        args = (self.default_factory,) if self.default_factory else ()
        return self.__class__, args, None, None, self.iteritems()


def _tryint(s):
    try:
        return int(s)
    except ValueError:
        return s


def _alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [_tryint(c) for c in re.split('([0-9]+)', s)]


def human_sorted(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=_alphanum_key)
    return l


def get_chr_lengths_from_seq(seq_fpath):
    chr_lengths = []

    if verify_file(seq_fpath + '.fai', silent=True):
        info('Reading genome index file (.fai) to get chromosome lengths')
        with open(adjust_path(seq_fpath + '.fai'), 'r') as handle:
            for line in handle:
                line = line.strip()
                if line:
                    chrom, length = line.split()[0], line.split()[1]
                    chr_lengths.append((chrom, length))
    elif verify_file(seq_fpath, silent=True):
        info('Reading genome sequence (.fa) to get chromosome lengths')
        with open(adjust_path(seq_fpath), 'r') as handle:
            from Bio import SeqIO
            reference_records = SeqIO.parse(handle, 'fasta')
            for record in reference_records:
                chrom = record.id
                chr_lengths.append((chrom, len(record.seq)))
    else:
        critical('Can\'t find ' + seq_fpath + ' and ' + seq_fpath + '.fai')
    return chr_lengths


def get_chr_len_fpath(cnf):
    chr_len_fpath = join(cnf.work_dir, 'chr_lengths.txt')
    if cnf.reuse_intermediate and file_exists(chr_len_fpath):
        info(chr_len_fpath + ' exists, reusing')
        return chr_len_fpath

    else:
        if not cnf.genome.seq:
            critical('There is no "seq" key in ' + cnf.sys_cnf + ' for "' + cnf.genome.name + '" section')
            return None

        chr_lengths = get_chr_lengths_from_seq(adjust_path(cnf.genome.seq))

        with file_transaction(cnf.work_dir, chr_len_fpath) as tx:
            with open(tx, 'w') as handle:
                for c, l in chr_lengths:
                    handle.write(c + '\t' + str(l) + '\n')
    return chr_len_fpath


def get_ext_tools_dirname(is_common_file=False):
    from sys import platform as _platform
    if not is_common_file and 'darwin' in _platform:
        return join('ext_tools', 'osx')
    else:
        return 'ext_tools'


def get_chr_lengths(cnf):
    chr_len_fpath = get_chr_len_fpath(cnf)
    if not chr_len_fpath:
        return None

    chr_lengths = []
    with open(chr_len_fpath, 'r') as f:
        for line in f:
            if len(line.split()) == 2:
                chr_name = line.split()[0]
                chr_length = int(line.split()[1])
                chr_lengths.append([chr_name, chr_length])
    return chr_lengths


def format_integer(name, value, unit=''):
    value = int(value)
    if value is not None:
        return '{name}: {value:,}{unit}'.format(**locals())
    else:
        return '{name}: -'.format(**locals())


def format_decimal(name, value, unit=''):
    if value is not None:
        return '{name}: {value:.2f}{unit}'.format(**locals())
    else:
        return '{name}: -'.format(**locals())


def mean(values):
    return float(sum(values)) / len(values) if len(values) > 0 else float('nan')


def median(values):
    values = sorted(values)

    if len(values) % 2 == 1:  # odd number of values
        return values[(len(values) - 1) / 2]
    else:  # even number of values - take the avg of central
        return (values[len(values) / 2] + values[len(values) / 2 - 1]) / 2


def get_numeric_value(string_value):
    """ parses string_value and returns only number-like part
    """
    num_chars = ['.', '+', '-']
    number = ''
    for c in string_value:
        if c.isdigit() or c in num_chars:
            number += c
    return number


def is_uk():
    hostname = socket.gethostname()
    return 'ukap' in hostname

def is_china():
    hostname = socket.gethostname()
    return 'cniclhpc' in hostname

def is_local():
    hostname = socket.gethostname()
    return 'local' in hostname or 'Home' in hostname or environ.get('PYTHONUNBUFFERED')

def is_us():
    hostname = socket.gethostname()
    return any(name in hostname for name in ['rask', 'chara', 'blue', 'green', 'espo', 'orr', 'usbod', 'bn0'])

def is_az():
    return is_us() or is_uk() or is_china()

def is_cloud():
    hostname = socket.gethostname()
    return 'starcluster' in hostname


def md5(fpath):
    hash = hashlib.md5()
    with open(fpath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash.update(chunk)
    return hash.hexdigest()