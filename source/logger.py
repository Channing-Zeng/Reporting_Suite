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
cnf_address = None

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

def send_email(msg_other='', subj='', only_me=False):
    # for addr in [my_address, ]
    if not msg_other or not smtp_host:
        return

    prid = getpass.getuser()

    me_address = email_by_prid['klpf990']

    other_addresses = []
    if not only_me:
        if cnf_address:
            other_addresses.append(cnf_address)
        if prid in email_by_prid:
            other_addresses.append(email_by_prid[prid])
    other_addresses = [a for a in other_addresses if a !=  me_address]

    msg_other += '\n'
    msg_other += '\n'
    msg_other += 'Ran by ' + prid + '\n'
    msg_other += '\n'
    if critical_msgs:
        msg_other += 'Critical errors during the processing:\n'
        for m in critical_msgs:
            msg_other += '  ' + m + '\n'
        msg_other += '\n'

    if error_msgs:
        if critical_msgs:
            msg_other += 'Other e'
        else:
            msg_other += 'E'
        msg_other += 'rrors during the processing:\n'
        for m in error_msgs:
            msg_other += '  ' + m + '\n'

    msg_me = msg_other[:]
    if warning_msgs:
        msg_me += 'Warnings during the processing:\n'
        for m in warning_msgs:
            msg_me += '  ' + m + '\n'

    if not subj:
        if project_name:
            subj += project_name
        else:
            subj += 'Reporting'
        if proc_name:
            subj += ' - ' + proc_name

    msg_other = MIMEText(msg_other)
    msg_me = MIMEText(msg_me)

    msg_other['Subject'] = msg_me['Subject'] = subj

    msg_other['From'] = msg_me['From'] = 'klpf990@rask.usbod.astrazeneca.com'
    msg_other['To'] = ','.join(other_addresses)
    msg_me['To'] = me_address

    def try_send(host, msg_):
        s = smtplib.SMTP(host)
        s.sendmail(msg_['From'], msg_['To'].split(','), msg_.as_string())
        s.quit()
        info('Mail sent to ' + msg_['To'] + ' using ' + host)

    def print_msg():
        for line in msg_other.as_string().split('\n'):
            print '   | ' + line
        print ''

    for msg in [msg_other, msg_me]:
        try:
            try_send(smtp_host, msg)
        except:
            warn('Could not send email using the sever "' + smtp_host + '" with exception: ')
            warn('  ' + '; '.join(traceback.format_exception_only(sys.exc_type, sys.exc_value)))
            if smtp_host != 'localhost':
                warn('Trying "localhost" as a server...')
                try:
                    try_send('localhost', msg)
                except:
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
    msg = str(msg)
    
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
            open(log_fpath, 'a').write(msg + ending)
        except IOError:
            sys.stderr.write('Logging: cannot append to ' + log_fpath + '\n')
            try:
                open(log_fpath, 'w').write(msg + '\n')
            except IOError:
                sys.stderr.write('Logging: cannot write to ' + log_fpath + '\n')