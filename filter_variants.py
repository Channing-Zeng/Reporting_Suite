#!/usr/bin/env python
import subprocess

import sys
from source.config import Defaults
from source.logger import err

if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from source.main import read_opts_and_cnfs, check_system_resources, load_genome_resources, check_inputs
from source.runner import run_one
from source.utils import info, get_java_tool_cmdline, intermediate_fname, call, iterate_file


def main(args):
    cnf = read_opts_and_cnfs(
        extra_opts=[
            (['--expression'], dict(
                dest='expression',
                help='filtering line for snpsift'
            )),
            (['-f'], 'DOUBLE', dict(
                dest='freq',
                type='double',
                help='When the ave allele frequency is also below the [freq], '
                     'the variant is considered likely false positive. '
                     'Default %d. Used with -r and -n' % Defaults.variant_filtering.freq,
            )),
        ],
        required_keys=['vcf'],
        file_keys=['vcf'],
        key_for_sample_name='vcf')

    check_system_resources(cnf,
                           required=['java', 'snpsift'],
                           optional=[])

    if 'expression' in cnf:
        cnf['variant_filtering'] = {'expression': cnf['expression']}

    run_one(cnf, process_one, finalize_one)


def process_one(cnf):
    anno_vcf_fpath = filter_variants(cnf, cnf['vcf'])
    return [anno_vcf_fpath]


def finalize_one(cnf, filtered_vcf_fpath):
    if filtered_vcf_fpath:
        err('Saved filtered VCF to ' + filtered_vcf_fpath)


def filter_variants(cnf, vcf_fpath):
    executable = get_java_tool_cmdline(cnf, 'snpsift')
    expression = cnf['variant_filtering']['expression'] or ''

    err('')
    err('*' * 70)
    vcf_fpath = iterate_file(cnf, vcf_fpath, lambda l: l.replace('\tPASS\t', '\t\t'))

    # awk_sort = subprocess.Popen( ["-c", "awk -f script.awk | sort > outfile.txt" ],
    # stdin= subprocess.PIPE, shell=True )
    # awk_sort.communicate( "input data\n" )
    # awk_sort.wait()

    cmdline = '{executable} filter -i PASS -f {vcf_fpath} "{expression}"'.format(**locals())
    filtered_fpath = intermediate_fname(cnf, vcf_fpath, 'filtered')
    subprocess.call(cmdline, stdin=sys.stdin, stdout=)
    return filtered_fpath


if __name__ == '__main__':
    main(sys.argv[1:])