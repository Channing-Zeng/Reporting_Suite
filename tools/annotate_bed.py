#!/usr/bin/env python

from os.path import abspath, dirname, realpath, join
from site import addsitedir
project_dir = abspath(dirname(dirname(realpath(__file__))))
addsitedir(join(project_dir))
addsitedir(join(project_dir, 'ext_modules'))
import sub_scripts.__check_python_version  # do not remove it: checking for python version and adding site dirs inside

from collections import defaultdict
import sys
from source.logger import err, is_local


def main():
    """
    Usage:
        annotate_bed.py [Reference_BED_file] < Input_BED_file > Annotated_BED_file
    """

    ref_bed_fpath = '/ngs/reference_data/genomes/Hsapiens/hg19/bed/Exons/Exons.bed'
    if len(sys.argv) > 1:
        ref_bed_fpath = sys.argv[1]

        bedtools = get_system_path(cnf, 'bedtools')
    cmdline = 'cut -f1,2,3 {amplicons_bed} ' \
              '| {bedtools} intersect -a - -b {exons_bed} -loj'
    call(cnf, cmdline, output_fpath)


    with sys.stdin as inp, sys.stdout as out:
        for l in inp:
            if l.startswith('#'):
                out.write(l)
            else:


            not_approved_gene_names = _proc_ensembl(inp, out, approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym)
        else:
            not_approved_gene_names = _proc_ucsc(inp, out, approved_gene_by_name, approved_gnames_by_prev_gname, approved_gnames_by_synonym)

    sys.stderr.write('\n')
    sys.stderr.write('Not approved ' + str(len(not_approved_gene_names.keys())) + ' genes.\n')
    if not_approved_fpath:
        with open(not_approved_fpath, 'w') as f:
            f.write('#Searched as\tStatus\n')
            f.writelines((gn + '\t' + st + '\n' for gn, st in not_approved_gene_names.items()))
        sys.stderr.write('Saved not approved to ' + not_approved_fpath + '\n')


if __name__ == '__main__':
    main()