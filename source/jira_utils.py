from traceback import format_exc
from source.logger import err, info, warn


JIRA_SERVER = 'https://jira.rd.astrazeneca.net'


class JiraCase:
    def __init__(self, case_id, url, assignee=None, reporter=None, type_=None, department=None, division=None,
                 data_hub=None, analysis_path=None, project_name=None, project_id=None, summary=None, description=None):
        self.case_id = case_id
        self.url = url
        self.assignee = assignee
        self.reporter = reporter
        self.type = type_
        self.department = department
        self.division = division
        self.data_hub = data_hub
        self.analysis_path = analysis_path
        self.project_name = project_name
        self.project_id = project_id
        self.summary = summary
        self.description = description


def retrieve_jira_info(url):
    try:
        from jira import JIRA
    except ImportError, e:
        warn('Cannot import JIRA: ' + str(e))
        return None

    """
    :param jira_url:  https://jira.rd.astrazeneca.net/i#browse/NGSG-38
                      https://jira.rd.astrazeneca.net/browse/NGSG-38
                      https://jira.rd.astrazeneca.net/i#browse/NGSG-38?filter=-1
    :return: instance of JiraCase
    """
    try:
        info('Creating JIRA instance...')
        jira_inst = JIRA(server=JIRA_SERVER, basic_auth=('NGSG_user', 'todngs'), options={'verify': False})
    except:
        warn('Cannot create JIRA obj:')
        warn(format_exc())
        return None

    info('Adding JIRA proxy...')
    jira_inst._session.proxies = {
        'http': 'http://hpcproxy.usbod.astrazeneca.net:3128',
        'https': 'http://hpcproxy.usbod.astrazeneca.net:3128'
    }

    info('Parsing JIRA url ' + str(url))
    case_id = __parse_id(url)
    if case_id:
        info('Parsing the JIRA case ' + case_id)
    else:
        err('Could not parse JIRA case from ' + str(url) + ', skipping connecting to JIRA.')
        return None

    info('Getting issue ' + case_id)
    try:
        issue = jira_inst.issue('NGSG-' + case_id)
    except:
        return None

    case = JiraCase(case_id=case_id, url=url)
    # print issue.fields.project.key             # 'JRA'
    case.reporter = issue.fields.reporter.displayName    # 'Greenawalt, Danielle'
    case.assignee = issue.fields.assignee.displayName    # 'Saif, Sakina'
    case.summary = issue.fields.summary              # Bio_029 - M2Gen - 100 Post Treatment Exomes
    case.type = issue.fields.customfield_12711.value if issue.fields.customfield_12711 else None  # Exome/Panel/etc
    case.department = issue.fields.customfield_12701.value if issue.fields.customfield_12701 else None  # BIO/EXT/etc
    # case.division = None  # always 'ONC' in NGS.Project.csv, no such field in Jira
    case.data_hub = issue.fields.customfield_12704  # projects/ProcessedDataHub/Patients/BRCA/Bio_029_M2Gen_RR
    case.analysis_path = issue.fields.customfield_12714  # /analysis/bioscience/Bio_031_M2Gen_RR
    case.project_name = issue.fields.customfield_12707  # M2Gen_RR
    case.project_id = issue.fields.customfield_12709  # Bio_031_M2Gen_RR
    case.description = issue.fields.description  # 100 post treatment exomes sourced...

    if case.reporter:
        info('reporter: ' + case.reporter)
    if case.assignee:
        info('assignee: ' + case.assignee)
    if case.summary:
        info('summary: ' + case.summary)
    if case.type:
        info('type: ' + case.type)
    if case.department:
        info('department: ' + case.department)
    if case.data_hub:
        info('data_hub: ' + case.data_hub)
    if case.analysis_path:
        info('analysis_path: ' + case.analysis_path)
    if case.project_name:
        info('project_name: ' + case.project_name)
    if case.project_id:
        info('project_id: ' + case.project_id)
    info()

    return case


def __parse_id(url):
    t = url.split('NGSG-')
    if len(t) == 1:
        err('Incorrect JIRA URL ' + url)
        return None
    case_id = t[1].split('?')[0]
    return case_id
