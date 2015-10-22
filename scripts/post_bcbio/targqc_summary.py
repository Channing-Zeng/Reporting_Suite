#!/usr/bin/env python
import bcbio_postproc

import sys
from source.bcbio.bcbio_structure import BCBioStructure, summary_script_proc_params
from source.logger import info
from source.targetcov.summarize_targetcov import summarize_targqc, get_bed_targqc_inputs


def main():
    info(' '.join(sys.argv))
    info()

    cnf, bcbio_structure = summary_script_proc_params(
        BCBioStructure.targqc_name,
        BCBioStructure.targqc_summary_dir,
        extra_opts=[
           (['--bed', '--capture', '--amplicons'], dict(
                dest='bed',
                help='BED file to run targetSeq and Seq2C analysis on.')
            ),
           (['--exons', '--exome'], dict(
                dest='exons',
                help='Exons BED file to make targetSeq exon/amplicon regions reports.')
            )
        ])

    bed_fpath, exons_bed_fpath, genes_fpath = get_bed_targqc_inputs(cnf, bcbio_structure.bed)

    summarize_targqc(cnf, cnf.threads or len(bcbio_structure.samples),
        cnf.output_dir, bcbio_structure.samples, bed_fpath, exons_bed_fpath, genes_fpath)


if __name__ == '__main__':
    main()