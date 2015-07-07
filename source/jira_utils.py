from ext_modules.jira import JIRA


class JiraCase:
    def __init__(self, url, assignee=None, reporter=None, type_=None, department=None, data_hub=None):
        self.url = url
        self.assignee = assignee
        self.reporter = reporter
        self.type = type_
        self.department = department
        self.data_hub = data_hub


def retrieve_jira_info(jira_url):
    jira = JIRA(jira_url)  # https://jira.rd.astrazeneca.net/i#browse/NGSG-38
                           # https://jira.rd.astrazeneca.net/browse/NGSG-196
                           # https://jira.rd.astrazeneca.net/i#browse/NGSG-38?filter=-1
    t = jira_url.split('NGSG-')
    if len(t) == 1:
        return None
    case_id = t[1].split('?')[0]
    issue = jira.issue('NGSG-' + case_id)
    # retrieve everything
    case = JiraCase(jira_url)
    # print issue.fields.project.key             # 'JRA'
    # print issue.fields.issuetype.name          # 'New Feature'
    case.reporter = issue.fields.reporter.displayName    # 'Mike Cannon-Brookes [Atlassian]'
    case.type = issue.fields.project.type
    case.department = issue.fields.project.group
    case.data_hub = issue.fields.data_hub_location
    return case