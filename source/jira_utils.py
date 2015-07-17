from traceback import format_exc
from source.logger import err


JIRA_SERVER = 'https://jira.rd.astrazeneca.net'


class JiraCase:
    def __init__(self, case_id, url, assignee=None, reporter=None, type_=None, department=None, division=None, data_hub=None):
        self.case_id = case_id
        self.url = url
        self.assignee = assignee
        self.reporter = reporter
        self.type = type_
        self.department = department
        self.division = division
        self.data_hub = data_hub


def retrieve_jira_info(url):
    from ext_modules.jira import JIRA

    """
    :param jira_url:  https://jira.rd.astrazeneca.net/i#browse/NGSG-38
                      https://jira.rd.astrazeneca.net/browse/NGSG-196
                      https://jira.rd.astrazeneca.net/i#browse/NGSG-38?filter=-1
    :return: instance of JiraCase
    """
    try:
        jira_inst = JIRA(server=JIRA_SERVER,
                    basic_auth=('klpf990', '123qweasd'),
                    options={'verify': False})
    except:
        err(format_exc())
        return None

    # retrieve everything
    case_id = __parse_id(url)
    issue = jira_inst.issue('NGSG-' + case_id)
    case = JiraCase(case_id=case_id, url=url)
    # print issue.fields.project.key             # 'JRA'
    case.reporter = issue.fields.reporter.displayName    # 'Mike Cannon-Brookes [Atlassian]'
    case.assignee = issue.fields.assignee.displayName    # 'Mike Cannon-Brookes [Atlassian]'
    case.description = issue.fields.summary              # HiSeq4000_2x75 Whole genome sequencing of 5 AURA plasma
    # case.type = issue.fields.project.type
    # case.department = issue.fields.project.group
    # case.data_hub = issue.fields.data_hub_location
    # case.division = None
    return case


def __parse_id(url):
    t = url.split('NGSG-')
    if len(t) == 1:
        err('Incorrect JIRA URL ' + url)
        return None
    case_id = t[1].split('?')[0]
    return case_id