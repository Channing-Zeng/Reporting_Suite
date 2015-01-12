from os.path import join
from collections import OrderedDict

varfilter_name           = 'varFilter'
varannotate_name         = 'varAnnotate'
targetseq_name           = 'targetSeq'
cnv_dir                  = 'cnv'

varqc_name               = 'varQC'
varqc_after_name         = 'varQC_postVarFilter'
ngscat_name              = 'ngscat'
qualimap_name            = 'qualimap'
targqc_name              = 'targQC'
fastqc_name              = 'fastqc'

varqc_repr               = 'Var QC'
varqc_after_repr         = 'Var QC after filtering'
ngscat_repr              = 'Ngscat'
qualimap_repr            = 'Qualimap'
targqc_repr              = 'Target QC'
fastqc_repr              = 'FastQC'


class BaseSample:
    def __init__(self, name, dirpath=None, bam=None, bed=None, vcf=None, genome=None):
        self.name = name
        self.bam = bam
        self.bed = bed
        self.qualimap_bed = None
        self.vcf_by_callername = OrderedDict()  # string -> vcf_fpath
        self.vcf = vcf
        self.dirpath = dirpath
        self.phenotype = None
        self.genome = None
        self.var_dirpath = None
        self.normal_match = None
        self.min_af = None

    def make_fpath(self, path_template, **kwargs):
        keys = dict(dict(dirpath=self.dirpath, sample=self.name).items() + kwargs.items())
        return join(*path_template.split('/')).format(**keys)

