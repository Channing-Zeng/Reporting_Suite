"""Partly taken from bcbio-nextgen.
"""

import shutil
import os
from os.path import isfile, isdir, getsize, exists, basename, join, abspath, splitext, islink, dirname, \
    pardir, realpath
import gzip
import tempfile
import contextlib
import itertools
import functools
import random
import ConfigParser
import collections
import sys
import fnmatch
import time
from ext_modules import yaml

from source.logger import info, err, warn, critical


try:
    from concurrent import futures
except ImportError:
    try:
        import futures
    except ImportError:
        futures = None


@contextlib.contextmanager
def cpmap(cores=1):
    """Configurable parallel map context manager.

    Returns appropriate map compatible function based on configuration:
    - Local single core (the default)
    - Multiple local cores
    """
    if int(cores) == 1:
        yield itertools.imap
    else:
        if futures is None:
            raise ImportError("concurrent.futures not available")
        pool = futures.ProcessPoolExecutor(cores)
        yield pool.map
        pool.shutdown()

def map_wrap(f):
    """Wrap standard function to easily pass into 'map' processing.
    """
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return apply(f, *args, **kwargs)
    return wrapper


def safe_mkdir(dirpath, descriptive_name=''):
    if isdir(dirpath):
        return dirpath

    if not dirpath:
        critical(descriptive_name + ' path is empty.')

    if isfile(dirpath):
        critical(descriptive_name + ' ' + dirpath + ' is a file.')

    num_tries = 0
    max_tries = 10

    while not exists(dirpath):
        # we could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(dirpath)
        except OSError, e:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(2)
    return dirpath


def transform_to(ext):
    """
    Decorator to create an output filename from an output filename with
    the specified extension. Changes the extension, in_file is transformed
    to a new type.

    Takes functions like this to decorate:
    f(in_file, out_dir=None, out_file=None) or,
    f(in_file=in_file, out_dir=None, out_file=None)

    examples:
    @transform(".bam")
    f("the/input/path/file.sam") ->
        f("the/input/path/file.sam", out_file="the/input/path/file.bam")

    @transform(".bam")
    f("the/input/path/file.sam", out_dir="results") ->
        f("the/input/path/file.sam", out_file="results/file.bam")

    """

    def decor(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            out_file = kwargs.get("out_file", None)
            if not out_file:
                in_path = kwargs.get("in_file", args[0])
                out_dir = kwargs.get("out_dir", os.path.dirname(in_path))
                safe_mkdir(out_dir)
                out_name = replace_suffix(os.path.basename(in_path), ext)
                out_file = os.path.join(out_dir, out_name)
            kwargs["out_file"] = out_file
            if not file_exists(out_file):
                out_file = f(*args, **kwargs)
            return out_file
        return wrapper
    return decor


def filter_to(word):
    """
    Decorator to create an output filename from an input filename by
    adding a word onto the stem. in_file is filtered by the function
    and the results are written to out_file. You would want to use
    this over transform_to if you don't know the extension of the file
    going in. This also memoizes the output file.

    Takes functions like this to decorate:
    f(in_file, out_dir=None, out_file=None) or,
    f(in_file=in_file, out_dir=None, out_file=None)

    examples:
    @filter_to(".foo")
    f("the/input/path/file.sam") ->
        f("the/input/path/file.sam", out_file="the/input/path/file.foo.bam")

    @filter_to(".foo")
    f("the/input/path/file.sam", out_dir="results") ->
        f("the/input/path/file.sam", out_file="results/file.foo.bam")

    """

    def decor(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            out_file = kwargs.get("out_file", None)
            if not out_file:
                in_path = kwargs.get("in_file", args[0])
                out_dir = kwargs.get("out_dir", os.path.dirname(in_path))
                safe_mkdir(out_dir)
                out_name = append_stem(os.path.basename(in_path), word)
                out_file = os.path.join(out_dir, out_name)
            kwargs["out_file"] = out_file
            if not file_exists(out_file):
                out_file = f(*args, **kwargs)
            return out_file
        return wrapper
    return decor


def memoize_outfile(ext=None, stem=None):
    """
    Memoization decorator.

    See docstring for transform_to and filter_to for details.
    """
    if ext:
        return transform_to(ext)
    if stem:
        return filter_to(stem)


@contextlib.contextmanager
def chdir(new_dir):
    """Context manager to temporarily change to a new directory.

    http://lucentbeing.com/blog/context-managers-and-the-with-statement-in-python/
    """
    cur_dir = os.getcwd()
    safe_mkdir(new_dir)
    os.chdir(new_dir)
    try:
        yield
    finally:
        os.chdir(cur_dir)

def file_uptodate(fname, cmp_fname):
    """Check if a file exists, is non-empty and is more recent than cmp_fname.
    """
    return (file_exists(fname) and file_exists(cmp_fname) and
            os.path.getmtime(fname) >= os.path.getmtime(cmp_fname))

def create_dirs(config, names=None):
    if names is None:
        names = config["dir"].keys()
    for dname in names:
        d = config["dir"][dname]
        safe_mkdir(d)

def save_diskspace(fname, reason, config):
    """Overwrite a file in place with a short message to save disk.

    This keeps files as a sanity check on processes working, but saves
    disk by replacing them with a short message.
    """
    if config["algorithm"].get("save_diskspace", False):
        with open(fname, "w") as out_handle:
            out_handle.write("File removed to save disk space: %s" % reason)

def read_galaxy_amqp_config(galaxy_config, base_dir):
    """Read connection information on the RabbitMQ server from Galaxy config.
    """
    galaxy_config = add_full_path(galaxy_config, base_dir)
    config = ConfigParser.ConfigParser()
    config.read(galaxy_config)
    amqp_config = {}
    for option in config.options("galaxy_amqp"):
        amqp_config[option] = config.get("galaxy_amqp", option)
    return amqp_config

def add_full_path(dirname, basedir=None):
    if basedir is None:
        basedir = os.getcwd()
    if not dirname.startswith("/"):
        dirname = os.path.join(basedir, dirname)
    return dirname

def symlink_plus(orig, new):
    """Create relative symlinks and handle associated biological index files.
    """
    for ext in ["", ".idx", ".gbi", ".tbi", ".bai"]:
        if os.path.exists(orig + ext) and not os.path.lexists(new + ext):
            with chdir(os.path.dirname(new)):
                os.symlink(os.path.relpath(orig + ext), os.path.basename(new + ext))
    orig_noext = splitext_plus(orig)[0]
    new_noext = splitext_plus(new)[0]
    for sub_ext in [".bai"]:
        if os.path.exists(orig_noext + sub_ext) and not os.path.lexists(new_noext + sub_ext):
            with chdir(os.path.dirname(new_noext)):
                os.symlink(os.path.relpath(orig_noext + sub_ext), os.path.basename(new_noext + sub_ext))

def open_gzipsafe(f, mode='rb'):
    if f.endswith('.gz'):
        try:
            h = gzip.open(f, mode=mode)
        except IOError, e:
            err('Error opening gzip ' + f + ': ' + str(e) + ', opening as plain text')
            return open(f, mode=mode)
        else:
            try:
                h.read(1)
            except IOError, e:
                err('Error opening gzip ' + f + ': ' + str(e) + ', opening as plain text')
                h.close()
                return open(f, mode=mode)
            else:
                h.close()
                h = gzip.open(f, mode=mode)
                return h
    else:
        return open(f, mode=mode)


def append_stem(to_transform, word):
    """
    renames a filename or list of filenames with 'word' appended to the stem
    of each one:
    example: append_stem("/path/to/test.sam", "_filtered") ->
    "/path/to/test_filtered.sam"

    """
    if is_sequence(to_transform):
        return [append_stem(f, word) for f in to_transform]
    elif is_string(to_transform):
        (base, ext) = splitext_plus(to_transform)
        return "".join([base, word, ext])
    else:
        raise ValueError("append_stem takes a single filename as a string or "
                         "a list of filenames to transform.")


def replace_suffix(to_transform, suffix):
    """
    replaces the suffix on a filename or list of filenames
    example: replace_suffix("/path/to/test.sam", ".bam") ->
    "/path/to/test.bam"

    """
    if is_sequence(to_transform):
        transformed = []
        for f in to_transform:
            (base, _) = os.path.splitext(f)
            transformed.append(base + suffix)
        return transformed
    elif is_string(to_transform):
        (base, _) = os.path.splitext(to_transform)
        return base + suffix
    else:
        raise ValueError("replace_suffix takes a single filename as a string or "
                         "a list of filenames to transform.")

# ## Functional programming

def partition_all(n, iterable):
    """Partition a list into equally sized pieces, including last smaller parts
    http://stackoverflow.com/questions/5129102/python-equivalent-to-clojures-partition-all
    """
    it = iter(iterable)
    while True:
        chunk = list(itertools.islice(it, n))
        if not chunk:
            break
        yield chunk

def partition(pred, iterable):
    'Use a predicate to partition entries into false entries and true entries'
    # partition(is_odd, range(10)) --> 0 2 4 6 8   and  1 3 5 7 9
    t1, t2 = itertools.tee(iterable)
    return itertools.ifilterfalse(pred, t1), itertools.ifilter(pred, t2)

# ## Dealing with configuration files

def merge_config_files(fnames):
    """Merge configuration files, preferring definitions in latter files.
    """
    def _load_yaml(fname):
        with open(fname) as in_handle:
            config = yaml.load(in_handle)
        return config
    out = _load_yaml(fnames[0])
    for fname in fnames[1:]:
        cur = _load_yaml(fname)
        for k, v in cur.iteritems():
            if out.has_key(k) and isinstance(out[k], dict):
                out[k].update(v)
            else:
                out[k] = v
    return out


def get_in(d, t, default=None):
    """
    look up if you can get a tuple of values from a nested dictionary,
    each item in the tuple a deeper layer

    example: get_in({1: {2: 3}}, (1, 2)) -> 3
    example: get_in({1: {2: 3}}, (2, 3)) -> {}
    """
    result = reduce(lambda d, t: d.get(t, {}), t, d)
    if not result:
        return default
    else:
        return result


def flatten(l):
    """
    flatten an irregular list of lists
    example: flatten([[[1, 2, 3], [4, 5]], 6]) -> [1, 2, 3, 4, 5, 6]
    lifted from: http://stackoverflow.com/questions/2158395/

    """
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el,
                                                                   basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el


def is_sequence(arg):
    """
    check if 'arg' is a sequence

    example: arg([]) -> True
    example: arg("lol") -> False

    """
    return (not hasattr(arg, "strip") and
            hasattr(arg, "__getitem__") or
            hasattr(arg, "__iter__"))


def is_pair(arg):
    """
    check if 'arg' is a two-item sequence

    """
    return is_sequence(arg) and len(arg) == 2

def is_string(arg):
    return isinstance(arg, basestring)


def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
    supplied root directory.'''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)


def itersubclasses(cls, _seen=None):
    """
    snagged from:  http://code.activestate.com/recipes/576949/
    itersubclasses(cls)

    Generator over all subclasses of a given class, in depth first order.

    >>> list(itersubclasses(int)) == [bool]
    True
    >>> class A(object): pass
    >>> class B(A): pass
    >>> class C(A): pass
    >>> class D(B,C): pass
    >>> class E(D): pass
    >>>
    >>> for cls in itersubclasses(A):
    ...     print(cls.__name__)
    B
    D
    E
    C
    >>> # get ALL (new-style) classes currently defined
    >>> [cls.__name__ for cls in itersubclasses(object)] #doctest: +ELLIPSIS
    ['type', ...'tuple', ...]
    """

    if not isinstance(cls, type):
        raise TypeError('itersubclasses must be called with '
                        'new-style classes, not %.100r' % cls)
    if _seen is None:
        _seen = set()
    try:
        subs = cls.__subclasses__()
    except TypeError:  # fails only when cls is type
        subs = cls.__subclasses__(cls)
    for sub in subs:
        if sub not in _seen:
            _seen.add(sub)
            yield sub
            for sub in itersubclasses(sub, _seen):
                yield sub

def replace_directory(out_files, dest_dir):
    """
    change the output directory to dest_dir
    can take a string (single file) or a list of files
    """
    if is_sequence(out_files):
        filenames = map(os.path.basename, out_files)
        return [os.path.join(dest_dir, x) for x in filenames]
    elif is_string(out_files):
        return os.path.join(dest_dir, os.path.basename(out_files))
    else:
        raise ValueError("in_files must either be a sequence of filenames "
                         "or a string")

def which(program):
    """
    returns the path to an executable or None if it can't be found
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def reservoir_sample(stream, num_items, item_parser=lambda x: x):
    """
    samples num_items from the stream keeping each with equal probability
    """
    kept = []
    for index, item in enumerate(stream):
        if index < num_items:
            kept.append(item_parser(item))
        else:
            r = random.randint(0, index)
            if r < num_items:
                kept[r] = item_parser(item)
    return kept


def compose(f, g):
    return lambda x: f(g(x))

def dictapply(d, fn):
    """
    apply a function to all non-dict values in a dictionary
    """
    for k, v in d.items():
        if isinstance(v, dict):
            v = dictapply(v, fn)
        else:
            d[k] = fn(v)
    return d


def verify_module(name):
    try:
        __import__(name)
        return True
    except:
        return False


def adjust_path(path):
    if path is None: return None

    path = remove_quotes(path)
    if path is None: return None

    path = expanduser(path)
    if path is None: return None

    path = abspath(path)
    if path is None: return None

    path = realpath(path)
    if path is None: return None

    return path


code_base_path = abspath(join(dirname(abspath(__file__)), pardir))

def adjust_system_path(path):
    if path is None: return None

    path = remove_quotes(path)
    if path is None: return None

    path = expanduser(path)
    if path is None: return None

    path = join(code_base_path, path)  # will only join if the tool_path is not absolute:
    if path is None: return None

    path = realpath(path)
    if path is None: return None

    path = abspath(path)
    if path is None: return None

    return path


# Expand paths beginning with '~' or '~user'.
# '~' means $HOME; '~user' means that user's home directory.
# If the path doesn't begin with '~', or if the user or $HOME is unknown,
# the path is returned unchanged (leaving error reporting to whatever
# function is called with the expanded path as argument).
# See also module 'glob' for expansion of *, ? and [...] in pathnames.
# (A function should also be defined to do full *sh-style environment
# variable expansion.)
def expanduser(path):
    """Expand ~ and ~user constructs.

    If user or $HOME is unknown, do nothing."""
    if path[:1] != '~':
        return path
    i, n = 1, len(path)
    while i < n and path[i] not in '/\\':
        i = i + 1

    if 'HOME' in os.environ:
        userhome = os.environ['HOME']
    elif 'USERPROFILE' in os.environ:
        userhome = os.environ['USERPROFILE']
    elif not 'HOMEPATH' in os.environ:
        return path
    else:
        try:
            drive = os.environ['HOMEDRIVE']
        except KeyError:
            drive = ''
        userhome = join(drive, os.environ['HOMEPATH'])

    if i != 1:  # ~user
        userhome = join(dirname(userhome), path[1:i])

    return userhome + path[i:]


def file_exists(fpath):
    """Check if a file exists and is non-empty.
    """
    return fpath and exists(adjust_path(fpath)) and getsize(adjust_path(fpath)) > 0


def _log(msg, silent, is_critical):
    if is_critical:
        critical(msg)
    if not silent:
        warn(msg)

def verify_obj_by_path(path, description='', silent=False, is_critical=False):
    if not path:
        msg = (description + ': f' if description else 'N') + 'ame is empty.'
        _log(msg, silent, is_critical)
        return path

    path = adjust_path(path)
    if not exists(path):
        msg = (description + ': ' if description else '') + path + ' does not exist.'
        _log(msg, silent, is_critical)
        return None

    if isfile(path):
        return verify_file(path, description, silent)
    elif isdir(path):
        return verify_dir(path, description, silent)
    else:
        msg = (description + ': ' if description else '') + path + ' is not a file or a directory.'
        _log(msg, silent, is_critical)
        return None

def verify_file(fpath, description='', silent=False, is_critical=False):
    if not fpath:
        msg = (description + ': f' if description else 'F') + 'ile name is empty.'
        _log(msg, silent, is_critical)
        return fpath

    fpath = adjust_path(fpath)
    if not exists(fpath):
        msg = (description + ': ' if description else '') + fpath + ' does not exist.'
        _log(msg, silent, is_critical)
        return None

    if not isfile(fpath):
        msg = (description + ': ' if description else '') + fpath + ' is not a file.'
        _log(msg, silent, is_critical)
        return None

    if getsize(fpath) <= 0:
        msg = (description + ': ' if description else '') + fpath + ' is empty.'
        _log(msg, silent, is_critical)
        return None

    return fpath

def verify_dir(fpath, description='', silent=False, is_critical=False):
    if not fpath:
        msg = (description + ': d' if description else 'D') + 'ir name is empty.'
        _log(msg, silent, is_critical)
        return None

    fpath = adjust_path(fpath)
    if not exists(fpath):
        msg = (description + ': ' if description else '') + fpath + ' does not exist.'
        _log(msg, silent, is_critical)
        return None

    if not isdir(fpath):
        msg = (description + ': ' if description else '') + fpath + ' is not a directory.'
        _log(msg, silent, is_critical)
        return None

    return fpath


def num_lines(fpath):
    with open(fpath) as f:
        return sum(1 for l in f)


# def make_tmpdir(cnf, prefix='ngs_reporting_tmp', *args, **kwargs):
#     base_dir = cnf.tmp_base_dir or cnf.work_dir or os.getcwd()
#     if not verify_dir(base_dir, 'Base directory for temporary files'):
#         sys.exit(1)
#
#     return tempfile.mkdtemp(dir=base_dir, prefix=prefix)


# @contextlib.contextmanager
# def tmpdir(cnf, *args, **kwargs):
#     prev_tmp_dir = cnf.tmp_dir
#
#     cnf.tmp_dir = make_tmpdir(cnf, *args, **kwargs)
#     try:
#         yield cnf.tmp_dir
#     finally:
#         try:
#             shutil.rmtree(cnf.tmp_dir)
#         except OSError:
#             pass
#         cnf.tmp_dir = prev_tmp_dir


def make_tmpfile(cnf, *args, **kwargs):
    yield tempfile.mkstemp(dir=cnf['work_dir'], *args, **kwargs)


@contextlib.contextmanager
def tmpfile(cnf, *args, **kwargs):
    tmp_file, fpath = make_tmpfile(cnf, *args, **kwargs)
    try:
        yield fpath
    finally:
        try:
            os.remove(fpath)
        except OSError:
            pass


def splitext_plus(fname):
    """Split on file extensions, allowing for zipped extensions.
    """
    base, ext = splitext(fname)
    if ext in [".gz", ".bz2", ".zip"]:
        base, ext2 = splitext(base)
        ext = ext2 + ext
    return base, ext


def add_suffix(fname, suf):
    base, ext = splitext_plus(fname)
    return base + '.' + suf + ext


def intermediate_fname(cnf, fpath, suf):
    output_fname = add_suffix(fpath, suf)
    return join(cnf.work_dir, basename(output_fname))


def remove_quotes(s):
    if s and s[0] in ['"', "'"]:
        s = s[1:]
    if s and s[-1] in ['"', "'"]:
        s = s[:-1]
    return s


def convert_file(cnf, input_fpath, convert_file_fn, suffix=None, check_result=True,
                 overwrite=False, reuse_intermediate=True, ctx=None):

    output_fpath = intermediate_fname(cnf, input_fpath, suf=suffix or 'tmp')
    if output_fpath.endswith('.gz'):
        output_fpath = output_fpath[:-3]

    if islink(output_fpath):
        os.unlink(output_fpath)

    if suffix and cnf.reuse_intermediate and reuse_intermediate and not overwrite and file_exists(output_fpath):
        info(output_fpath + ' exists, reusing')
        return output_fpath
    else:
        info('Writing to ' + output_fpath)

    with file_transaction(cnf.work_dir, output_fpath) as tx_fpath:
        with open_gzipsafe(input_fpath) as inp_f, open_gzipsafe(tx_fpath, 'w') as out_f:
            if ctx:
                convert_file_fn(inp_f, out_f, ctx)
            else:
                convert_file_fn(inp_f, out_f)

    if overwrite or suffix is None:
        shutil.move(output_fpath, input_fpath)
        output_fpath = input_fpath

    if suffix:
        info('Saved to ' + output_fpath)

    verify_file(output_fpath, is_critical=check_result)
    return output_fpath


def iterate_file(cnf, input_fpath, proc_line_fun, check_result=True, *args, **kwargs):
    def _proc_file(inp_f, out_f, ctx=None):
        max_bunch_size = 1000 * 1000
        written_lines = 0
        bunch = []

        for i, line in enumerate(inp_f):
            clean_line = line.strip()
            if clean_line:
                if ctx:
                    new_l = proc_line_fun(clean_line, i, ctx)
                else:
                    new_l = proc_line_fun(clean_line, i)
                if new_l is not None:
                    bunch.append(new_l + '\n')
                    written_lines += 1
            else:
                bunch.append(line)
                written_lines += 1

            if len(bunch) >= max_bunch_size:
                out_f.writelines(bunch)
                info('Written lines: ' + str(written_lines))
                bunch = []

        out_f.writelines(bunch)
        info('Written lines: ' + str(written_lines))

    return convert_file(cnf, input_fpath, _proc_file, check_result=check_result, *args, **kwargs)


def dots_to_empty_cells(config, tsv_fpath):
    """Put dots instead of empty cells in order to view TSV with column -t
    """
    def proc_line(l, i):
        while '\t\t' in l:
            l = l.replace('\t\t', '\t.\t')
        return l
    return iterate_file(config, tsv_fpath, proc_line, 'dots')


def __remove_tmpdirs(fnames):
    if isinstance(fnames, basestring):
        fnames = [fnames]
    for x in fnames:
        xdir = os.path.dirname(os.path.abspath(x))
        if xdir and os.path.exists(xdir):
            shutil.rmtree(xdir, ignore_errors=True)


def __remove_files(fnames):
    if isinstance(fnames, basestring):
        fnames = [fnames]

    for x in fnames:
        if x and os.path.exists(x):
            if os.path.isfile(x):
                os.remove(x)
            elif os.path.isdir(x):
                shutil.rmtree(x, ignore_errors=True)


#################################################
######## Transaction ############################
@contextlib.contextmanager
def file_transaction(work_dir, *rollback_files):
    """Wrap file generation in a transaction, moving to output if finishes.
    """
    exts = {".vcf": ".idx", ".bam": ".bai", "vcf.gz": ".tbi"}
    safe_fpaths, orig_names = _flatten_plus_safe(work_dir, rollback_files)
    __remove_files(safe_fpaths)  # remove any half-finished transactions
    try:
        if len(safe_fpaths) == 1:
            yield safe_fpaths[0]
        else:
            yield tuple(safe_fpaths)
    except:  # failure -- delete any temporary files
        __remove_files(safe_fpaths)
        raise
    else:  # worked -- move the temporary files to permanent location
        for safe, orig in zip(safe_fpaths, orig_names):
            if exists(safe):
                shutil.move(safe, orig)
                for check_ext, check_idx in exts.iteritems():
                    if safe.endswith(check_ext):
                        safe_idx = safe + check_idx
                        if exists(safe_idx):
                            shutil.move(safe_idx, orig + check_idx)


def _flatten_plus_safe(tmp_dir, rollback_files):
    """Flatten names of files and create temporary file names.
    """
    tx_fpaths, orig_files = [], []
    for fnames in rollback_files:
        if isinstance(fnames, basestring):
            fnames = [fnames]
        for fname in fnames:
            tx_file = add_suffix(fname, 'tx')
            tx_fpath = join(tmp_dir, tx_file)
            tx_fpaths.append(tx_fpath)
            orig_files.append(fname)
    return tx_fpaths, orig_files

#################################################