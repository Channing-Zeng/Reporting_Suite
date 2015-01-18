from os import environ
import socket
import sys
from datetime import datetime
from time import sleep
import smtplib
from email.mime.text import MIMEText
import traceback


log_fpath = None
project_name = None
project_fpath = None
proc_name = None
my_address = 'Vlad.Saveliev@astrazeneca.com'
address = None

import socket
hostname = socket.gethostname()
is_local = 'local' in hostname or 'Home' in hostname or environ.get('PYTHONUNBUFFERED')

smtp_host = None  # set up in source/config.py and system_info.yaml


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
    # if msg:
    #     send_email(msg)


def send_email(msg='', subj=''):
    return

    if msg and smtp_host:
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

        def try_send(host):
            s = smtplib.SMTP(host)
            s.sendmail(msg['From'], addresses, msg.as_string())
            s.quit()
            info('Mail sent to ' + msg['To'] + ' using ' + host)

        def print_msg():
            for line in msg.as_string().split('\n'):
                print '   | ' + line
            print ''

        try:
            try_send(smtp_host)
        except socket.error:
            warn('Could not send email using the sever "' + smtp_host + '" with exception: ')
            warn('; '.join(traceback.format_exception_only(sys.exc_type, sys.exc_value)))
            if smtp_host != 'localhost':
                warn('Trying "localhost" as a server...')
                try:
                    try_send('localhost')
                except socket.error:
                    warn('Could not send email using the sever "localhost" with exception: ')
                    warn('; '.join(traceback.format_exception_only(sys.exc_type, sys.exc_value)))
                    print_msg()
            else:
                print_msg()


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