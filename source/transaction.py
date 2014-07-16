"""Handle file based transactions allowing safe restarts at any point.

To handle interrupts,this defines output files written to temporary
locations during processing and copied to the final location when finished.
This ensures output files will be complete independent of method of
interruption.
"""
import os
import shutil
import contextlib
from os.path import exists, join
from source.utils_from_bcbio import add_suffix


def _remove_tmpdirs(fnames):
    for x in fnames:
        xdir = os.path.dirname(os.path.abspath(x))
        if xdir and os.path.exists(xdir):
            shutil.rmtree(xdir, ignore_errors=True)


def _remove_files(fnames):
    for x in fnames:
        if x and os.path.exists(x):
            if os.path.isfile(x):
                os.remove(x)
            elif os.path.isdir(x):
                shutil.rmtree(x, ignore_errors=True)


@contextlib.contextmanager
def file_transaction(cnf, *rollback_files):
    """Wrap file generation in a transaction, moving to output if finishes.
    """
    exts = {".vcf": ".idx", ".bam": ".bai", "vcf.gz": ".tbi"}
    safe_fpaths, orig_names = _flatten_plus_safe(cnf['work_dir'], rollback_files)
    _remove_files(safe_fpaths)  # remove any half-finished transactions
    try:
        if len(safe_fpaths) == 1:
            yield safe_fpaths[0]
        else:
            yield tuple(safe_fpaths)
    except:  # failure -- delete any temporary files
        _remove_files(safe_fpaths)
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