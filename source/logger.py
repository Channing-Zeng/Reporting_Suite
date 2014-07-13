from os import environ
import sys
from datetime import datetime
from time import sleep


log_fpath = None


def timestamp():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S  ")


def step_greetings(name):
    info()
    info('-' * 70)
    info(name)
    info('-' * 70)


def info(msg='', ending='\n'):
    _log(sys.stdout, msg, ending)


def err(msg='', ending='\n'):
    _log(sys.stderr, msg, ending)


def critical(msg=''):
    err(msg)
    sys.exit(1)


def _log(out, msg='', ending='\n'):
    msg = timestamp() + msg

    out.write(msg + ending)
    sys.stdout.flush()
    sys.stderr.flush()
    if environ.get('PYCHARM'):
        sleep(0.01)

    if log_fpath:
        try:
            open(log_fpath, 'a').write(msg + '\n')
        except IOError:
            sys.stderr.write('Cannot append to ' + log_fpath)
            try:
                open(log_fpath, 'w').write(msg + '\n')
            except IOError:
                sys.stderr.write('Cannot write to ' + log_fpath)