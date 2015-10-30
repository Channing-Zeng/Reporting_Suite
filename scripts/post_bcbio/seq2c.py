#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc


import sys
from os.path import join

from source.bcbio.bcbio_structure import BCBioStructure, summary_script_proc_params
from source.calling_process import call
from source.copy_number import cnv_reports
from source.logger import info
from source.tools_from_cnf import get_system_path


def main():
    info(' '.join(sys.argv))
    info()

    cnf, bcbio_structure = summary_script_proc_params(
        BCBioStructure.seq2c_name,
        BCBioStructure.cnv_summary_dir,
        extra_opts=[
           (['--controls', '-c'], dict(
                dest='controls',
                help='Optional control sample names for Seq2C. For multiple controls, separate them using :',
                default=''
           )),
           (['--seq2c_opts'], dict(
                dest='seq2c_opts',
                help='Options for the final lr2gene.pl script.',
                default=''
           )),
           (['--bed', '--capture', '--amplicons'], dict(
                dest='bed',
                help='BED file to run targetSeq and Seq2C analysis on.')
            ),
           (['--reannotate'], dict(
                dest='reannotate',
                help='re-annotate BED file with gene names',
                action='store_true',
                default=False)
            ),
           # (['--dedup'], dict(
           #      dest='dedup',
           #      help='Remove duplicates from the input bedfile.')
           #  ),
        ],
    )

    cnv_reports(cnf, bcbio_structure)

    # exons_bed_fpath = cnf.exons if cnf.exons else cnf.genome.exons
    # _, _, target_bed, seq2c_bed = prepare_beds(cnf, exons_bed=exons_bed_fpath, target_bed=bcbio_structure.sv_bed)
    #
    # bcbio_structure.sv_bed = target_bed or cnf.genome.refseq or seq2c_bed
    # standalone_cnv(cnf, bcbio_structure)
    # TODO: replace with /group/ngs/src/az.reporting/Seq2C/seq2c.sh sample2bam.tsv seq2c_regions.bed "" "" "-q batch.q"


def standalone_cnv(cnf, bcbio_structure):
    sample2bam_tsv_fpath = join(cnf.work_dir, 'seq2c_sample2bam.tsv')
    with open(sample2bam_tsv_fpath, 'w') as f:
        for s in bcbio_structure.samples:
            f.write(s.name + '\t' + s.bam + '\n')

    seq2c_sh = get_system_path(cnf, join('Seq2C', 'seq2c.sh'))
    samtools = get_system_path(cnf, 'samtools')
    cmdl = '{seq2c_sh} {sample2bam_tsv_fpath} {bcbio_structure.bed} {cnf.controls} {cnf.seq2c_opts} "-q {cnf.queue}" {samtools}'.format(**locals())
    seq2c_report_fpath = join(cnf.output_dir, BCBioStructure.seq2c_name + '.tsv')
    return call(cnf, cmdl, output_fpath=seq2c_report_fpath)


if __name__ == '__main__':
    main()

