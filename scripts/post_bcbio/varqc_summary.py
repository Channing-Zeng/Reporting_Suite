#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc


import sys
from collections import defaultdict

from source.variants.summarize_qc import make_summary_reports
from source.bcbio.bcbio_structure import BCBioStructure, bcbio_summary_script_proc_params
from source.logger import info


def main():
    info(' '.join(sys.argv))
    info()

    cnf, bcbio_structure = bcbio_summary_script_proc_params(BCBioStructure.varqc_name, BCBioStructure.varqc_dir)

    jsons_by_sample_by_caller = defaultdict(dict)
    htmls_by_sample_by_caller = defaultdict(dict)
    for vc in bcbio_structure.variant_callers.values():
        jsons_by_sample_by_caller[vc.name] = vc.find_fpaths_by_sample(cnf.proc_dir_name, cnf.proc_name, 'json', bcbio_structure.final_dirpath)
        htmls_by_sample_by_caller[vc.name] = vc.find_fpaths_by_sample(cnf.proc_dir_name, cnf.proc_name, 'html', bcbio_structure.final_dirpath)

    make_summary_reports(cnf, 1, cnf.output_dir, bcbio_structure.variant_callers.values(),
         bcbio_structure.samples, jsons_by_sample_by_caller, htmls_by_sample_by_caller,
         varqc_name=BCBioStructure.varqc_name, caption='Variant QC')


if __name__ == '__main__':
    main()

