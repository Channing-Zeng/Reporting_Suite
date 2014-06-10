import sys
from datetime import datetime


log_fpath = None


def timestamp():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S  ")


def step_greetings(name):
    info()
    info('-' * 70)
    info(name)
    info('-' * 70)


def info(msg=''):
    _log(sys.stdout, msg)


def err(msg=''):
    _log(sys.stderr, msg)


def critical(msg=''):
    err(msg)
    sys.exit(1)


def _log(out, msg=''):
    msg = timestamp() + msg

    out.write(msg + '\n')
    out.flush()

    if log_fpath:
        open(log_fpath, 'a').write(msg + '\n')