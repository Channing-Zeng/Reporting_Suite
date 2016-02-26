#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc


import sys
from source.bcbio.bcbio_structure import BCBioStructure, bcbio_summary_script_proc_params
from source.file_utils import adjust_path
from source.logger import info
from source.targetcov.summarize_targetcov import summarize_targqc, get_bed_targqc_inputs


def main():
    info(' '.join(sys.argv))
    info()

    cnf, bcbio_structure = bcbio_summary_script_proc_params(
        BCBioStructure.targqc_name,
        BCBioStructure.targqc_summary_dir,
        extra_opts=[
           (['--bed', '--capture', '--amplicons'], dict(
                dest='bed',
                help='BED file to run targetSeq and Seq2C analysis on.')
            ),
           (['--exons', '--exome', '--features'], dict(
                dest='features',
                help='Annotated CDS/Exons/Gene/Transcript BED file to make targetSeq exon/amplicon regions reports.')
            )
        ])

    bed_fpath, features_bed_fpath = adjust_path(cnf.bed), adjust_path(cnf.features)

    summarize_targqc(cnf, cnf.threads or len(bcbio_structure.samples),
        cnf.output_dir, bcbio_structure.samples, bed_fpath=bed_fpath, features_fpath=features_bed_fpath)


if __name__ == '__main__':
    main()