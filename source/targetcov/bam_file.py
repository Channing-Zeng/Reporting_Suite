from os.path import isfile
import sys
from source.calling_process import call
from source.logger import info
from source.tools_from_cnf import get_tool_cmdline


def index_bam(cnf, bam_fpath):
    indexed_bam = bam_fpath + '.bai'
    if not isfile(bam_fpath + '.bai'):
        info('Indexing to ' + indexed_bam + '...')
        samtools = get_tool_cmdline(cnf, 'samtools')
        if not samtools:
            sys.exit(1)
        cmdline = '{samtools} index {bam_fpath}'.format(**locals())
        call(cnf, cmdline)
    info('Index: ' + indexed_bam)