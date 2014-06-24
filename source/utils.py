import sys
import re

from collections import OrderedDict
from os.path import join, basename
from source.calling_process import call_subprocess
from source.tools_from_cnf import get_tool_cmdline
from source.logger import info
from source.utils_from_bcbio import file_exists


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


def remove_quotes(s):
    if s and s[0] == '"':
        s = s[1:]
    if s and s[-1] == '"':
        s = s[:-1]
    return s


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


def index_bam(cnf, bam_fpath):
    samtools = get_tool_cmdline(cnf, 'samtools')
    if not samtools:
        sys.exit(1)

    cmdline = '{samtools} index {bam_fpath}'.format(**locals())
    call_subprocess(cnf, cmdline, None, None)


def bgzip_and_tabix_vcf(cnf, vcf_fpath):
    bgzip = get_tool_cmdline(cnf, 'bgzip', suppress_warn=True)
    tabix = get_tool_cmdline(cnf, 'tabix', suppress_warn=True)

    gzipped_fpath = join(cnf['work_dir'], basename(vcf_fpath) + '.gz')
    tbi_fpath = gzipped_fpath + '.tbi'

    if bgzip and not file_exists(gzipped_fpath):
        info('BGzipping VCF')
        cmdline = '{bgzip} -c {vcf_fpath}'.format(**locals())
        call_subprocess(cnf, cmdline, None, gzipped_fpath, exit_on_error=False)

    if tabix and not file_exists(tbi_fpath):
        info('Tabixing VCF')
        cmdline = '{tabix} -f -p vcf {gzipped_fpath}'.format(**locals())
        call_subprocess(cnf, cmdline, None, tbi_fpath, exit_on_error=False)

    return gzipped_fpath, tbi_fpath


def get_chr_len_fpath(cnf):
    chr_len_fpath = cnf['genome'].get('chr_lengths')
    if chr_len_fpath:
        return chr_len_fpath

    chr_len_fpath = join(cnf['work_dir'], 'chr_lengths.txt')
    if cnf.reuse_intermediate and file_exists(chr_len_fpath):
        info(chr_len_fpath + ' exists, reusing')
        return chr_len_fpath

    chr_lengths = dict()
    with open(cnf['genome']['seq'], 'r') as handle:
        from Bio import SeqIO
        reference_records = SeqIO.parse(handle, 'fasta')
        for record in reference_records:
            chr_lengths[record.id] = len(record.seq)
    with open(chr_len_fpath, 'w') as handle:
        for chr_name in sorted(chr_lengths, key=chr_lengths.get, reverse=True):
            handle.write(chr_name + '\t' + str(chr_lengths[chr_name]) + '\n')
    return chr_len_fpath


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


def mean(ints):
    return float(sum(ints)) / len(ints) if len(ints) > 0 else float('nan')


def median(values):
    if len(values) % 2 == 1:
        return values[(len(values) - 1) / 2]
    else:
        return (values[(len(values) - 1) / 2] + values[(len(values) - 2) / 2]) / 2
