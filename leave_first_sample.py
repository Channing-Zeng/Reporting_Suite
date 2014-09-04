#!/usr/bin/env python
from genericpath import isfile
import sys
from source.file_utils import safe_mkdir
from source.main import load_genome_resources
from source.variants.vcf_processing import leave_first_sample
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from os.path import join, pardir, basename, dirname, abspath, realpath, islink, isdir
from site import addsitedir
source_dir = abspath(dirname(realpath(__file__)))
addsitedir(join(source_dir, 'ext_modules'))

import os
from source.variants.filtering import filter_for_variant_caller
from source.config import Defaults
from source.logger import info
from source.bcbio_structure import BCBioStructure
from source.summary import summary_script_proc_params


def main():
    info(' '.join(sys.argv))
    info()

    description = ''

    defaults = Defaults.variant_filtering

    cnf, bcbio_structure = summary_script_proc_params(
        BCBioStructure.varfilter_name,
        dir=BCBioStructure.varfilter_dir,
        description=description)

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

    filt_vcf_fpath_by_sample = caller.get_filt_vcf_by_samples()

    for sample, vcf_fpath in filt_vcf_fpath_by_sample.items():
        info(sample.name)
        info('Processing ' + vcf_fpath)
        new_vcf_fpath = leave_first_sample(cnf, vcf_fpath)
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

















