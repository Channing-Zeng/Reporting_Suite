from os.path import isfile
from source.calling_process import call
from source.tools_from_cnf import get_tool_cmdline


def index_bam(cnf, bam_fpath):
    if not isfile(bam_fpath + '.bai'):
        samtools = get_tool_cmdline(cnf, 'samtools')
        cmdline = '{samtools} index {bam}'.format(**locals())
        call(cnf, cmdline)