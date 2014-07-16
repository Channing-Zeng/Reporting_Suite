import contextlib
import hashlib
import sys
import tempfile
import os
import shutil
from os.path import isfile, isdir, getsize, exists, expanduser, basename, join, abspath

from source.logger import info, err
from source.transaction import file_transaction
from source.utils_from_bcbio import file_exists, open_gzipsafe, add_suffix


def verify_module(name):
    try:
        __import__(name)
        return True
    except ImportError:
        return False


def verify_file(fpath, description=''):
    if not fpath:
        err((description + ': f' if description else 'F') + 'ile name is empty.')
        return None

    fpath = abspath(expanduser(fpath))
    if not exists(fpath):
        err((description + ': ' if description else '') + fpath + ' does not exist.')
        return None

    if not isfile(fpath):
        err((description + ': ' if description else '') + fpath + ' is not a file.')
        return None

    if getsize(fpath) <= 0:
        err((description + ': ' if description else '') + fpath + ' is empty.')
        return None

    return fpath


def verify_dir(fpath, description=''):
    if not fpath:
        err((description + ': d' if description else 'D') + 'ir name is empty.')
        return None

    fpath = expanduser(fpath)
    if not exists(fpath):
        err((description + ': ' if description else '') + fpath + ' does not exist.')
        return None

    if not isdir(fpath):
        err((description + ': ' if description else '') + fpath + ' is not a directory.')
        return None

    return fpath


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


def intermediate_fname(cnf, fname, suf):
    output_fname = add_suffix(fname, suf)
    return join(cnf['work_dir'], basename(output_fname))


def convert_file(cnf, input_fpath, convert_file_fn, suffix=None,
                 overwrite=False, reuse_intermediate=True):

    output_fpath = intermediate_fname(cnf, input_fpath, suf=suffix or 'tmp')

    if suffix and cnf.reuse_intermediate and reuse_intermediate and not overwrite and file_exists(output_fpath):
        info(output_fpath + ' exists, reusing')
        return output_fpath
    else:
        info('Writing to ' + output_fpath)

    with file_transaction(cnf, output_fpath) as tx_fpath:
        with open(input_fpath) as inp_f, open(tx_fpath, 'w') as out_f:
            convert_file_fn(inp_f, out_f)

    if overwrite:
        shutil.move(output_fpath, input_fpath)
        output_fpath = input_fpath

    if suffix:
        info('Saved to ' + output_fpath)

    return output_fpath


def iterate_file(cnf, input_fpath, proc_line_fun, *args, **kwargs):
    def _proc_file(inp_f, out_f):
        max_bunch_size = 1000 * 1000
        written_lines = 0
        bunch = []

        for i, line in enumerate(inp_f):
            clean_line = line.strip()
            if clean_line:
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

    return convert_file(cnf, input_fpath, _proc_file, *args, **kwargs)


def dots_to_empty_cells(config, tsv_fpath):
    """Put dots instead of empty cells in order to view TSV with column -t
    """
    def proc_line(l, i):
        while '\t\t' in l:
            l = l.replace('\t\t', '\t.\t')
        return l
    return iterate_file(config, tsv_fpath, proc_line, 'dots')
