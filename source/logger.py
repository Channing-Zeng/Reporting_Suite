from collections import namedtuple
import getpass
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

smtp_host = None  # set up in source/config.py and system_info.yaml


error_msgs = []
warning_msgs = []
critical_msgs = []


def is_local():
    hostname = socket.gethostname()
    return 'local' in hostname or 'Home' in hostname or environ.get('PYTHONUNBUFFERED')


def timestamp():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S  ")


def step_greetings(name):
    info()
    info('-' * 70)
    info(name)
    info('-' * 70)


def info(msg='', ending='\n', print_date=True, severity='info'):
    _log(sys.stdout, msg, ending, print_date, severity=severity)


def warn(msg='', ending='\n', print_date=True, severity='warning'):
    _log(sys.stderr, msg, ending, print_date, severity=severity)


def silent_err(msg='', ending='\n', print_date=True, severity='silent_err'):
    warn(msg, ending, print_date, severity=severity)


def err(msg='', ending='\n', print_date=True, severity='error'):
    warn(msg, ending, print_date, severity=severity)


email_by_prid = {
    'klpf990': 'Vlad.Saveliev@astrazeneca.com',
    'kjgk478': 'Alexey.Gurevich@astrazeneca.com',
    'knfz728': 'Alla.Bushoy@astrazeneca.com',
    'klrl262': 'Miika.Ahdesmaki@astrazeneca.com',
    'kmtc481': 'Sakina.Saif@astrazeneca.com',
    'kxjn734': 'Justin.Johnson@astrazeneca.com'
}

def send_email(msg='', subj=''):
    if msg and smtp_host:
        addresses = [my_address]
        if address:
            addresses.append(address)
        prid = getpass.getuser()
        if prid in email_by_prid:
            addresses.append(email_by_prid[prid])

        msg += '\n'
        msg += '\n'
        msg += 'Ran by ' + prid + '\n'
        msg += '\n'
        if critical_msgs:
            msg += 'Critical errors during the processing:\n'
            for m in critical_msgs:
                msg += '  ' + m + '\n'
            msg += '\n'

        if error_msgs:
            if critical_msgs:
                msg += 'Other e'
            else:
                msg += 'E'
            msg += 'rrors during the processing:\n'
            for m in error_msgs:
                msg += '  ' + m + '\n'

        msg = MIMEText(msg)
        if not subj:
            if project_name:
                subj += project_name
            else:
                subj += 'Reporting'
            if proc_name:
                subj += ' - ' + proc_name
        msg['Subject'] = subj

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
            warn('  ' + '; '.join(traceback.format_exception_only(sys.exc_type, sys.exc_value)))
            if smtp_host != 'localhost':
                warn('Trying "localhost" as a server...')
                try:
                    try_send('localhost')
                except socket.error:
                    warn('Could not send email using the sever "localhost" with exception: ')
                    warn('  ' + '; '.join(traceback.format_exception_only(sys.exc_type, sys.exc_value)))
                    print_msg()
            else:
                print_msg()


def critical(msg=''):
    if isinstance(msg, basestring):
        err(msg, severity='critical')
    else:
        if not msg:
            return
        for m in msg:
            err(m, severity='critical')
    sys.exit(1)


def _log(out, msg='', ending='\n', print_date=True, severity=None):
    if print_date:
        msg = timestamp() + msg

    out.write(msg + ending)

    if severity == 'critical':
        critical_msgs.append(msg)
    if severity == 'error':
        error_msgs.append(msg)
    if severity == 'warning':
        warning_msgs.append(msg)

    sys.stdout.flush()
    sys.stderr.flush()
    if is_local():
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