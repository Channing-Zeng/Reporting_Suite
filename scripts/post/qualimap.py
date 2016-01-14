#!/usr/bin/python
# noinspection PyUnresolvedReferences
import bcbio_postproc


from os.path import isfile, join

from source.bcbio.bcbio_structure import BCBioStructure
from source.calling_process import call
from source.targetcov.bam_and_bed_utils import index_bam, number_of_dup_reads
from source.tools_from_cnf import get_system_path
from source.logger import info, critical
from source.main import read_opts_and_cnfs


def main():
    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--bam'], dict(
                dest='bam',
                help='path to the BAM file')
             ),
            (['--bed', '--capture', '--amplicons'], dict(
                dest='bed',
                help='capture panel/amplicons')
             ),
            (['--pcr'], dict(
                dest='pcr',
                action='store_true',
                help='deduplication was not perfomed, thus do not try to dedup')
             ),
        ],
        required_keys=['bam'],
        file_keys=['bam', 'bed'],
        key_for_sample_name='bam',
        proc_name=BCBioStructure.qualimap_name)

    if not isfile(cnf.bam + '.bai'):
        info('Indexing bam ' + cnf.bam)
        index_bam(cnf, cnf.bam)
    info('Using alignment ' + cnf.bam)

    bed = ''
    if cnf.bed:
        bed = ' -gff ' + cnf.bed + ' '
        info('Using amplicons/capture panel ' + cnf.bed)

    qualimap = get_system_path(cnf, 'qualimap', is_critical=True)
    if not qualimap:
        critical('Cannot find qualimap')

    info()
    cmdline = ('{qualimap} bamqc --skip-duplicated --skip-dup-mode 1 -nt ' + str(cnf.threads) + ' --java-mem-size=24G -nr 5000 '
        '-bam {cnf.bam} -outdir {cnf.output_dir} {bed} -c -gd HUMAN').format(**locals())
    report_fpath = join(cnf.output_dir, 'qualimapReport.html')

    call(cnf, cmdline, output_fpath=report_fpath, stdout_to_outputfile=False, env_vars=dict(DISPLAY=None))

    info('Qualimap report: ' + str(report_fpath))


if __name__ == '__main__':
    main()
