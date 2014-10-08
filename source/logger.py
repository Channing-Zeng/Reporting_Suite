from os import environ
import sys
from datetime import datetime
from time import sleep
import smtplib
from email.mime.text import MIMEText
import traceback


log_fpath = None
project_name = None
proc_name = None
address = 'vladislav.sav@gmail.com'


def timestamp():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S  ")


def step_greetings(name):
    info()
    info('-' * 70)
    info(name)
    info('-' * 70)


def info(msg='', ending='\n', print_date=True):
    _log(sys.stdout, msg, ending, print_date)


def warn(msg='', ending='\n', print_date=True):
    _log(sys.stdout, msg, ending, print_date)


def err(msg='', ending='\n', print_date=True):
    _log(sys.stderr, msg, ending, print_date)
    _send_email(msg)


def _send_email(msg=''):
    if msg:
        msg = MIMEText(msg)
        msg['Subject'] = 'Reporting for ' + (project_name or '') + ' ' +\
                         (proc_name or '')
        msg['From'] = 'AZ reporting'
        msg['To'] = address
        try:
            s = smtplib.SMTP('localhost')
            s.sendmail(msg['From'], [msg['To']], msg.as_string())
            s.quit()
        except:
            traceback.print_exc()
            warn('Could not send email with error. Proceeding.')
            pass


def critical(msg=''):
    if msg:
        err(msg)
    sys.exit(1)


def _log(out, msg='', ending='\n', print_date=True):
    if print_date:
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
            sys.stderr.write('Logging: cannot append to ' + log_fpath + '\n')
            try:
                open(log_fpath, 'w').write(msg + '\n')
            except IOError:
                sys.stderr.write('Logging: cannot write to ' + log_fpath + '\n')