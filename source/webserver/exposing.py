from collections import OrderedDict
import getpass
import os
from os.path import join, isfile, basename, dirname, abspath, isdir, relpath, realpath, pardir
from traceback import print_exc, format_exc
from source import verify_file
from source.file_utils import file_transaction, safe_mkdir, adjust_path
from source.logger import info, critical, err, is_local, warn
from source.tools_from_cnf import get_system_path, get_script_cmdline
from source.utils import is_uk, is_us, is_az
from source.webserver.ssh_utils import connect_to_server
from source.jira_utils import retrieve_jira_info

ngs_server_url = '172.18.47.33'
ngs_server_username = 'ngsall'
ngs_server_password = 'astra123'

class Location:
    def __init__(self, loc_id, report_url_base, website_url_base, csv_fpath, dirpath, reports_dirpath, proper_path_should_contain):
        self.loc_id = loc_id
        self.report_url_base = report_url_base
        self.website_url_base = website_url_base
        self.csv_fpath = csv_fpath
        self.dirpath = dirpath
        self.reports_dirpath = reports_dirpath
        self.proper_path_should_contain = proper_path_should_contain

us = Location('US',
    report_url_base='http://ngs.usbod.astrazeneca.net/reports/',
    website_url_base='http://ngs.usbod.astrazeneca.net/',
    csv_fpath='/ngs/oncology/NGS.Project.csv',
    dirpath='/opt/lampp/htdocs',
    reports_dirpath='/opt/lampp/htdocs/reports',
    proper_path_should_contain=['/gpfs/ngs/oncology/Analysis/', '/gpfs/ngs/oncology/Datasets/']
)
uk = Location('UK',
    report_url_base='http://ukapdlnx115.ukapd.astrazeneca.net/ngs/reports/',
    website_url_base='http://ngs.usbod.astrazeneca.net/',
    csv_fpath='/ngs/web_content/reports/NGS.Project.csv',
    dirpath='/ngs/web_content',
    reports_dirpath='/ngs/web_content/reports',
    proper_path_should_contain=['/ngs/oncology/analysis/', '/ngs/oncology/datasets/']
)
local = Location('Local',
    report_url_base='http://localhost/reports/',
    website_url_base='http://localhost/ngs_website/',
    csv_fpath='/Users/vlad/Sites/reports/NGS.Project.csv',
    dirpath='/Users/vlad/Sites',
    reports_dirpath='/Users/vlad/Sites/reports',
    proper_path_should_contain=['/Dropbox/az/analysis/', '/Dropbox/az/datasets/']
)
loc_by_id = dict(us=us, uk=uk, local=local)


def sync_with_ngs_server(
        cnf,
        jira_url,
        project_name,
        sample_names,
        summary_report_fpath,
        dataset_dirpath=None,
        bcbio_final_dirpath=None,
        jira_case=None):

    if is_us(): loc = us
    elif is_uk(): loc = uk
    elif is_local(): loc = local
    else:
        return None

    html_report_url = None
    if any(p in realpath((bcbio_final_dirpath or dataset_dirpath)) for p in loc.proper_path_should_contain):
        if jira_case is None and jira_case != 'unreached' and is_az() and jira_url:
            info()
            info('Getting info from JIRA...')
            jira_case = retrieve_jira_info(jira_url)

        proj_dirpath_on_server = _symlink_dirs(
            cnf=cnf,
            loc=loc,
            project_name=project_name,
            final_dirpath=bcbio_final_dirpath,
            dataset_dirpath=dataset_dirpath)
            # html_report_fpath=summary_report_fpath,
            # html_report_url=html_report_url)

        if bcbio_final_dirpath:
            html_report_url = join(loc.report_url_base,                                          # http://ngs.usbod.astrazeneca.net/reports/
                                   relpath(proj_dirpath_on_server, loc.reports_dirpath),         # project_name/dataset/project_name
                                   relpath(summary_report_fpath, dirname(bcbio_final_dirpath)))  # final/2015_01_01_project/project.html
        elif dataset_dirpath:
            html_report_url = join(loc.report_url_base,
                                   relpath(proj_dirpath_on_server, loc.reports_dirpath),
                                   relpath(summary_report_fpath, dataset_dirpath))

        html_report_full_url = join(loc.website_url_base, 'samples.php?project_name=' + project_name + '&file=' + html_report_url)
        info('HTML url: ' + html_report_full_url)

        if verify_file(loc.csv_fpath, 'Project list'):
            write_to_csv_file(
                work_dir=cnf.work_dir,
                jira_case=jira_case,
                project_list_fpath=loc.csv_fpath,
                country_id=loc.loc_id,
                project_name=project_name,
                samples_num=len(sample_names),
                analysis_dirpath=dirname(bcbio_final_dirpath) if bcbio_final_dirpath else None,
                html_report_url=html_report_url)
    return html_report_url


def _symlink_dirs(cnf, loc, project_name, final_dirpath, dataset_dirpath):  #, html_report_fpath, html_report_url):
    info(loc.loc_id + ', symlinking to ' + loc.reports_dirpath)
    dst_project_dirpath = None

    if dataset_dirpath:
        dst_project_dirpath = join(loc.reports_dirpath, project_name, 'dataset', project_name)
        (symlink_to_ngs if is_us() else local_symlink)(dataset_dirpath, dst_project_dirpath)

    if final_dirpath:
        dst_project_dirpath = join(loc.reports_dirpath, project_name, 'bcbio', project_name)  # before "final"!
        (symlink_to_ngs if is_us() else local_symlink)(dirname(final_dirpath), dst_project_dirpath)

    return dst_project_dirpath


def local_symlink(src, dst):
    if os.path.exists(dst):
        try:
            os.unlink(dst)
        except Exception, e:
            err('Cannot remove link ' + dst + ': ' + str(e))
            return None
    if not os.path.exists(dst):
        safe_mkdir(dirname(dst))
        try:
            os.symlink(src, dst)
        except Exception, e:
            err('Cannot create link ' + dst + ': ' + str(e))


# def symlink_uk(cnf, final_dirpath, project_name, dataset_dirpath, html_report_fpath):
#     html_report_url = UK_URL + project_name + '/' + relpath(html_report_fpath, final_dirpath)
#
#     info('UK, symlinking to ' + UK_SERVER_PATH)
#     link_fpath = join(UK_SERVER_PATH, project_name)
#     cmd = 'rm ' + link_fpath + '; ln -s ' + final_dirpath + ' ' + link_fpath
#     info(cmd)
#     try:
#         os.system(cmd)
#     except Exception, e:
#         warn('Cannot create symlink')
#         warn('  ' + str(e))
#         html_report_url = None
#     return html_report_url


# def symlink_us(cnf, final_dirpath, project_name, dataset_dirpath, html_report_fpath):
#     html_report_url = None
#     ssh = connect_to_server()
#     if ssh is None:
#         return None
#     html_report_url = US_URL + project_name + '/' +   relpath(html_report_fpath, final_dirpath)
#     final_dirpath_in_ngs = realpath(final_dirpath).split('/gpfs')[1]
#     link_path = join(US_SERVER_PATH, project_name)
#     cmd = 'rm ' + link_path + '; ln -s ' + final_dirpath_in_ngs + ' ' + link_path
#     ssh.exec_command(cmd)
#     info('  ' + cmd)
#     ssh.close()
#
#     info()
#     return html_report_url


def symlink_to_ngs(src_path, dst_fpath):
    ssh = connect_to_server(ngs_server_url, ngs_server_username, ngs_server_password)
    if ssh is None:
        return None

    src_path = src_path.replace('/gpfs/', '/')
    for cmd in ['mkdir ' + dirname(dirname(dst_fpath)),
                'mkdir ' + dirname(dst_fpath),
                'rm ' + dst_fpath,
                'ln -s ' + src_path + ' ' + dst_fpath,
                'chmod -R g+w ' + dirname(dirname(dst_fpath))]:
        info('Executing on the server:  ' + cmd)
        try:
            ssh.exec_command(cmd)
        except Exception, e:
            err('Cannot execute command: ' + str(e))
        continue
    info('Symlinked ' + src_path + ' to ' + dst_fpath)
    ssh.close()

    return dst_fpath


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
            values = map(__unquote, l.split(','))
            pid = values[index_of_pid]
            values_by_keys_by_pid[pid] = OrderedDict(zip(header_keys, values))

    pid = project_name
    with file_transaction(work_dir, project_list_fpath) as tx_fpath:
        if pid not in values_by_keys_by_pid.keys():
            # info(pid + ' not in ' + str(values_by_keys_by_pid.keys()))
            info('Adding new record for ' + pid)
            values_by_keys_by_pid[pid] = OrderedDict(zip(header_keys, [''] * len(header_keys)))
        else:
            info('Updating existing record for ' + pid)
        d = values_by_keys_by_pid[pid]
        for k in header_keys:
            if k not in d:
                err('Error: ' + k + ' not in ' + project_list_fpath + ' for ' + pid)

        d['PID'] = pid

        if analysis_dirpath:
            d['Analyses directory ' + (country_id if not is_local() else 'US')] = analysis_dirpath
        if project_name and (analysis_dirpath or not __unquote(d['Name'])):  # update only if running after bcbio, or no value there at all
            d['Name'] = project_name
        if html_report_url and (analysis_dirpath or not __unquote(d['HTML report path'])):  # update only if running after bcbio, or no value there at all
            d['HTML report path'] = html_report_url

        if jira_case:
            d['JIRA URL'] = jira_case
            # if 'Updated By' in d and __unquote(d['Updated By']):
            d['Updated By'] = getpass.getuser()
            if jira_case.description:
                d['Description'] = jira_case.summary
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
        if samples_num:
            d['Sample Number'] = str(samples_num)

        new_line = ','.join(__requote(d.get(k, '').replace(',', ';').replace('\n', ' | ')) or '' for k in header_keys)

        with open(tx_fpath, 'w') as f:
            os.umask(0002)
            try:
                os.chmod(tx_fpath, 0666)
            except OSError:
                err(format_exc())
            for l in lines:
                if not l:
                    pass
                if l.startswith('#'):
                    f.write(l)
                else:
                    l = unicode(l, 'utf-8')
                    l_ascii = l.encode('ascii', 'ignore')
                    if ',' + project_name + ',' in l_ascii or ',"' + project_name + '",' in l_ascii:
                        info('Old csv line: ' + l_ascii)
                        # f.write('#' + l)
                    else:
                        f.write(l)
            f.write(new_line + '\n')
        info()
        info('New line: ' + new_line)
        info()


def __unquote(s):
    if s.startswith('"') or s.startswith("'"):
        s = s[1:]
    if s.endswith('"') or s.endswith("'"):
        s = s[:-1]
    return s


def __requote(s):
    if s.startswith('"') or s.startswith("'"):
        s = s[1:]
    if s.endswith('"') or s.endswith("'"):
        s = s[:-1]
    s = s.replace('"', "'")
    return '"' + s + '"'


def convert_path_to_url(path):
    if is_us():
        return adjust_path(path).replace('/gpfs/ngs/', 'http://blue.usbod.astrazeneca.net/~klpf990/ngs/')
    else:
        return ''

