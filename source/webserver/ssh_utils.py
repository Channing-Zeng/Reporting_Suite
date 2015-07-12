from collections import OrderedDict
import getpass
import os
from os.path import join, isfile, basename, dirname, abspath, isdir, relpath
from traceback import print_exc
from source import verify_file
from source.file_utils import file_transaction
from source.logger import info, critical, err, is_local, warn
from source.tools_from_cnf import get_system_path, get_script_cmdline
from source.utils import is_uk, compatible_with_ngs_webserver


US_CSV = '/ngs/oncology/NGS.Project.csv'
UK_CSV = '/ngs/oncology/reports/NGS.Project.csv'


def sync_with_ngs_server(cnf, jira_case, project_name, sample_names, final_dirpath=None, final_summary_report_fpath=None):
    if not compatible_with_ngs_webserver():
        return None

    info()
    csv_fpath = US_CSV
    country_id = 'US'
    fn = _symlink_report_us
    if is_uk():
        csv_fpath = UK_CSV
        country_id = 'UK'
        fn = _symlink_report_uk

    # Symlink reports
    html_report_url = fn(cnf, cnf.work_dir, final_dirpath, project_name, final_summary_report_fpath)

    # Write to CSV
    if verify_file(csv_fpath, 'Project list'):
        write_to_csv_file(
            work_dir=cnf.work_dir, jira_case=jira_case, project_list_fpath=csv_fpath,
            country_id=country_id, project_name=project_name, samples_num=len(sample_names),
            analysis_dirpath=dirname(final_dirpath) if final_dirpath else None,
            html_report_url=final_summary_report_fpath)

    return html_report_url


def _symlink_report_uk(cnf, work_dir, final_dirpath, project_name, html_report_fpath):
    html_report_url = 'http://ukapdlnx115.ukapd.astrazeneca.net/ngs/reports/' + project_name + '/' + \
        relpath(html_report_fpath, final_dirpath)

    server_path = '/ngs/oncology/reports'
    info('UK, symlinking to ' + server_path)
    link_fpath = join(server_path, project_name)
    cmd = 'rm ' + link_fpath + '; ln -s ' + final_dirpath + ' ' + link_fpath
    info(cmd)
    try:
        os.system(cmd)
    except Exception, e:
        warn('Cannot create symlink')
        warn('  ' + str(e))
        html_report_url = None
    return html_report_url


def _symlink_report_us(cnf, work_dir, final_dirpath, project_name, html_report_fpath):
    server_path = '/opt/lampp/htdocs/reports'

    html_report_url = None
    with connect_to_server() as ssh:
        html_report_url = 'http://ngs.usbod.astrazeneca.net/reports/' + project_name + '/' + \
            relpath(html_report_fpath, final_dirpath)
        final_dirpath_in_ngs = final_dirpath.split('/gpfs')[1]
        link_path = join(server_path, project_name)
        cmd = 'rm ' + link_path + '; ln -s ' + final_dirpath_in_ngs + ' ' + link_path
        ssh.exec_command(cmd)
        info('  ' + cmd)

    info()
    return html_report_url


def symlink_to_ngs(src_fpaths, dst_dirpath):
    if isinstance(src_fpaths, basestring):
        src_fpaths = [src_fpaths]

    dst_fpaths = []

    with connect_to_server() as ssh:
        for src_fpath in src_fpaths:
            dst_fpath = join(dst_dirpath, basename(src_fpath))
            for cmd in ['mkdir ' + dst_dirpath,
                        'rm ' + dst_fpath,
                        'ln -s ' + src_fpath + ' ' + dst_fpath]:
                ssh.exec_command(cmd)
                info('  ' + cmd)

            info('Symlinked ' + src_fpath + ' to ' + dst_fpath)
            dst_fpaths.append(dst_fpath)

    if len(dst_fpaths) == 1:
        return dst_fpaths[0]
    else:
        return dst_fpaths


def write_to_csv_file(work_dir, jira_case, project_list_fpath, country_id, project_name,
                      samples_num=None, analysis_dirpath=None, html_report_url=None):
    info('Reading project list ' + project_list_fpath)
    with open(project_list_fpath) as f:
        lines = f.readlines()
    uncom_lines = [l.strip() for l in lines if not l.strip().startswith('#')]

    header = uncom_lines[0].strip()
    info('header: ' + header)
    header_keys = header.split(',')  # 'Updated By,PID,Name,JIRA URL,HTML report path,Why_IfNoReport,Data Hub,Analyses directory UK,Analyses directory US,Type,Division,Department,Sample Number,Reporter,Assignee,Description,IGV,Notes'
    index_of_pid = header_keys.index('PID')
    if index_of_pid == -1: index_of_pid = 1

    values_by_keys_by_pid = OrderedDict()
    for l in uncom_lines[1:]:
        if l:
            values = l.split(',')
            pid = values[index_of_pid]
            values_by_keys_by_pid[pid] = OrderedDict(zip(header_keys, values))

    pid = project_name
    with file_transaction(work_dir, project_list_fpath) as tx_fpath:
        if pid not in values_by_keys_by_pid:
            info('Adding new record for ' + pid)
            values_by_keys_by_pid[pid] = OrderedDict(zip(header_keys, [''] * len(header_keys)))
        else:
            info('Updating existing record for ' + pid)
        d = values_by_keys_by_pid[pid]

        d['PID'] = pid
        d['Name'] = project_name
        if jira_case:
            d['JIRA URL'] = jira_case.url
            d['Updated By'] = getpass.getuser() if 'Updated By' not in d else d['Updated By']
            if jira_case.description:
                d['Description'] = jira_case.description
            if jira_case.data_hub:
                d['Data Hub'] = jira_case.data_hub
            if jira_case.type:
                d['Type'] = jira_case.type
            if jira_case.department:
                d['Department'] = jira_case.department
            if jira_case.division:
                d['Division'] = jira_case.division
            if jira_case.assignee:
                d['Assignee'] = jira_case.assignee
            if jira_case.reporter:
                d['Reporter'] = jira_case.reporter
        if html_report_url:
            d['HTML report path'] = html_report_url
        if analysis_dirpath:
            d['Analyses directory ' + country_id] = analysis_dirpath
        if samples_num:
            d['Sample Number'] = str(samples_num)

        new_line = ','.join(v or '' for v in d.values())
        info('Adding the new line: ' + new_line)

        with open(tx_fpath, 'w') as f:
            os.umask(0002)
            os.chmod(tx_fpath, 0666)
            for l in lines:
                if not l:
                    pass
                if l.startswith('#'):
                    f.write(l)
                else:
                    if ',' + project_name + ',' in l:
                        info('Old csv line: ' + l)
                        # f.write('#' + l)
                    else:
                        f.write(l)
            f.write(new_line + '\n')


def connect_to_server(server_url='172.18.47.33', username='klpf990', password='123werasd'):
    # html_report_url = 'http://ngs.usbod.astrazeneca.net/reports/' + bcbio_structure.project_name + '/' + \
    #     relpath(html_report_fpath, bcbio_structure.final_dirpath)

    rsa_key_path = get_system_path(None, join(dirname(__file__), 'id_rsa'), is_critical=False)
    if not rsa_key_path:
        err('Could not find key ' + rsa_key_path)

    try:
        from ext_modules.paramiko import SSHClient, RSAKey, AutoAddPolicy
    except ImportError as e:
        print_exc()
        warn()
        err('Cannot improt SSHClient - skipping trasnferring symlinking to the ngs-website')
        warn()
    else:
        ssh = SSHClient()
        ssh.load_system_host_keys()
        # ki = RSAKey.from_private_key_file(filename=rsa_key_path)
        ssh.set_missing_host_key_policy(AutoAddPolicy())
        try:
            key = RSAKey(filename=rsa_key_path, password='%1!6vLaD')
        except Exception, e:
            warn('Cannot read RSAKey from ' + rsa_key_path)
            warn()
            print_exc()
            warn()
        else:
            info('Succesfully read RSAKey from ' + rsa_key_path)
            try:
                ssh.connect(server_url, username=username, password=password, pkey=key)
            except Exception, e:
                warn('Cannot connect to ' + server_url + ':')
                warn()
                print_exc()
                warn()
            else:
                info('Succesfully connected to ' + server_url)
                return ssh

    return None