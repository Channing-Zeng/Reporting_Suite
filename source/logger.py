from os import environ
import socket
import sys
from datetime import datetime
from time import sleep
import smtplib
from email.mime.text import MIMEText
import traceback


smtp_host = 'localhost'
log_fpath = None
project_name = None
project_fpath = None
proc_name = None
my_address = 'vladislav.sav@gmail.com'
address = None

import socket
hostname = socket.gethostname()
is_local = 'local' in hostname


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
    _log(sys.stderr, msg, ending, print_date)


def err(msg='', ending='\n', print_date=True):
    warn(msg, ending, print_date)
    if msg:
        send_email(msg)


def send_email(msg='', subj=''):
    if msg:
        addresses = [my_address]
        if address:
            addresses.append(address)

        msg = MIMEText(msg)
        subject = ''
        if project_name:
            subject += project_name
        else:
            subject += 'Reporting'
        if project_fpath:
            subject += ' - ' + project_fpath
        if proc_name:
            subject += ' - ' + proc_name
        if subj:
            subject += ': ' + subj
        msg['Subject'] = subject

        msg['From'] = 'klpf990@rask.usbod.astrazeneca.com'
        msg['To'] = ','.join(addresses)
        try:
            s = smtplib.SMTP(smtp_host)
            s.sendmail(msg['From'], addresses, msg.as_string())
            s.quit()
            info('Mail sent to ' + msg['To'])
        except socket.error:
            warn('Could not send email with exception: ')
            warn('; '.join(traceback.format_exception_only(sys.exc_type, sys.exc_value)))
            for line in msg.as_string().split('\n'):
                print '   | ' + line
            print ''


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
    if is_local:
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