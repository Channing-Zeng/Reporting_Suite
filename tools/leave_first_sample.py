#!/usr/bin/env python
import sys
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from os.path import join, dirname, abspath, realpath
from site import addsitedir
source_dir = abspath(dirname(realpath(__file__)))
addsitedir(join(source_dir, 'ext_modules'))

import os
from source.logger import info
from source.bcbio_structure import BCBioStructure
from source.summary import summary_script_proc_params


def main():
    info(' '.join(sys.argv))
    info()

    cnf, bcbio_structure = summary_script_proc_params(
        BCBioStructure.varfilter_name,
        dir=BCBioStructure.varfilter_dir)

    info('*' * 70)
    info()

    proc_all(cnf, bcbio_structure)


def proc_all(cnf, bcbio_structure):
    info('Starting.')
    info('-' * 70)
    for _, caller in bcbio_structure.variant_callers.items():
        proc_var_caller(caller, cnf, bcbio_structure)


def proc_var_caller(caller, cnf, bcbio_structure):
    info('Running for ' + caller.name)

    filt_vcf_fpath_by_sample = caller.find_pass_filt_vcf_by_samples()

    for sample, vcf_fpath in filt_vcf_fpath_by_sample.items():
        info(sample.name)
        info('Processing ' + vcf_fpath)
        new_vcf_fpath = vcf_fpath + '.tx'

        with open(vcf_fpath) as f, open(new_vcf_fpath, 'w') as out:
            for l in f:
                if l.startswith('#CHROM'):
                    ts = l.strip().split('\t')
                    l2 = '\t'.join(ts[:10])
                    info('  ' + l.strip() + '    ->    ' + l2)
                    l = l2 + '\n'
                out.write(l)

        # new_vcf_fpath = leave_first_sample(cnf, vcf_fpath)

        info('Saved to ' + new_vcf_fpath)
        if new_vcf_fpath != vcf_fpath:
            os.remove(vcf_fpath)
            os.rename(new_vcf_fpath, vcf_fpath)
        info('Renamed into ' + vcf_fpath)
        info()

    info('-' * 70)
    info()

    return caller

if __name__ == '__main__':
    main()

















