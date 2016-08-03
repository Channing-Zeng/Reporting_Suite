import re

from source import verify_file
from source.targetcov.Region import SortableByChrom
from source.targetcov.bam_and_bed_utils import get_gene_keys


class SVEvent(SortableByChrom):
    min_sv_depth = 15

    class Annotation:
        KNOWN = 1
        ON_PRIORITY_LIST = 2
        NOT_PRIORITISED = 3

        def __init__(self):
            self.type = None
            self.effect = None
            self.genes = []
            self.transcript = None
            self.priority = None
            self.prioritisation = None

            self.known = False
            self.event = None

        def get_key(self):
            if self.event:
                return self.event.get_key(), self.type, self.effect, self.transcript

        @staticmethod
        def parse_annotation(string):
            def filter_digits(s):
                return ''.join(c for c in s if c.isdigit())

            fs = string.split('|')
            a = SVEvent.Annotation()
            a.type = fs[0]
            a.effect = fs[1]
            genes_val = fs[2]
            if a.type == 'BND':
                if '/' in genes_val:
                    a.genes = sorted(genes_val.split('/'))
                elif '&' in genes_val:
                    a.genes = sorted(genes_val.split('&'))
                elif genes_val.count('-') == 1:
                    a.genes = sorted(genes_val.split('-'))
                else:
                    return None
            elif genes_val:
                if '&' in genes_val:
                    a.genes = sorted(genes_val.split('&'))
                else:
                    a.genes = [genes_val]
            a.transcript = fs[3]
            a.exon_info = fs[4]
            if len(fs) > 5:
                a.exon_info = filter_digits(a.exon_info)
                a.prioritisation = fs[4]
                a.priority = int(fs[5])
            return a

    @staticmethod
    def parse_sv_event(chr_order, key_gene_by_name_chrom, **kwargs):  # caller  sample  chrom  start  end  svtype  lof  annotation  split_read_support  paired_support_PE  paired_support_PR
        e = SVEvent(chrom=kwargs.get('chrom'), chrom_ref_order=chr_order.get(kwargs.get('chrom')))
        e.caller = kwargs.get('caller')
        e.start = int(kwargs.get('start'))
        e.sample = kwargs.get('sample')
        e.end = int(kwargs.get('end')) if kwargs.get('end') else None

        e.type = None
        e.id = None
        e.mate_id = None
        svt = kwargs.get('svtype')
        if svt:  # BND:MantaBND:12:0:1:0:0:0:1:MantaBND:12:0:1:0:0:0:0 or BND:71_2:71_1 or
            e.type = svt.split(':', 1)[0]
            if e.type == 'BND' and ':' in svt:
                if 'MantaBND' in svt:
                    m = re.match(r'(?P<id1>MantaBND[:0-9]+):(?P<id2>MantaBND[:0-9]+)', svt.split(':', 1)[1])
                else:
                    m = re.match(r'(?P<id1>.+):(?P<id2>.+)', svt.split(':', 1)[1])
                e.id = m.group('id1')
                e.mate_id = m.group('id2')

        # e.known_gene_val = kwargs.get('known')
        # e.known_gene = e.known_gene_val.split('-with-')[1] if '-with-' in e.known_gene_val else e.known_gene_val
        if kwargs.get('end_gene'):
            e.end_genes = kwargs.get('end_gene').split(',')
        else:
            e.end_genes = []
        e.lof = kwargs.get('lof')
        e.annotations = []
        if kwargs.get('annotation'):
            for s in kwargs.get('annotation').split(','):
                a = SVEvent.Annotation.parse_annotation(s)
                if a:
                    assert a.type == e.type, 'Annotation type and event type does not match: ' + str(e.type) + ', ' + str(a.type)
                    e.annotations.append(a)
                    known_genes = [g for g in a.genes if (g, e.chrom) in key_gene_by_name_chrom]
                    if known_genes:
                        e.known_gene = known_genes[0]
        paired_end_manta_header = 'paired_support_PR' if 'paired_support_PR' in kwargs else 'paired_end_support'
        paired_end_lumpy_header = 'paired_support_PE' if 'paired_support_PE' in kwargs else 'paired_end_support'
        if e.caller == 'manta':  # Manta has comma separated REF/ALT depths in third last and last column
            e.split_read_support = int(kwargs.get('split_read_support').split(',')[1]) if kwargs.get('split_read_support') else None
            e.paired_end_support = int(kwargs.get(paired_end_manta_header).split(',')[1]) if kwargs.get(paired_end_manta_header) else None
        else:
            e.split_read_support = int(kwargs.get('split_read_support').split(',')[0]) if kwargs.get('split_read_support') else None
            e.paired_end_support = int(kwargs.get(paired_end_lumpy_header).split(',')[0]) if kwargs.get(paired_end_lumpy_header) else None
        return e

        # lof_genes = []
        # if kwargs.get('end_gene'):
        #     for a in kwargs.get('end_gene').split(','):
        #         lof_genes.append(a[1:-1].split('|')[0])
        # with open(adjust_path('~/t.tsv'), 'a') as f:
        #     f.write(str(self.type) + '\t' + kwargs.get('known') + '\t' + kwargs.get('end_gene') + '\t' + ', '.join(lof_genes) + '\n')

    def __init__(self, chrom, chrom_ref_order):
        SortableByChrom.__init__(self, chrom, chrom_ref_order)
        self.caller = None
        self.start = None
        self.sample = None
        self.end = None
        self.chrom2 = None
        self.type = None
        self.id = None
        self.mate_id = None
        self.known_gene_val = None
        self.known_gene = ''
        self.end_genes = []
        self.lof = None
        self.annotations = []
        self.key_annotations = set()
        self.split_read_support = None
        self.paired_end_support = None

    def is_known_fusion(self, annotation):
        return annotation.prioritisation and 'KNOWN_FUSION' in annotation.prioritisation

    def is_fusion(self):
        return self.type == 'BND'

    def get_possible_fusion_pairs(self):
        if self.type == 'BND':
            return [a.genes for a in self.annotations]

    def is_deletion(self):
        return self.type == 'DEL'

    def is_insertion(self):
        return self.type == 'INS'

    def is_duplication(self):
        return self.type == 'DUP'

    def __str__(self):
        return str(self.chrom) + ':' + str(self.start) + ' ' + str(self.annotations)

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash((self.caller, self.chrom, self.start, self.type, self.id, self.mate_id))

    def get_key(self):
        return self.chrom, self.type  #, tuple(tuple(sorted(a.genes)) for a in self.annotations)

    def get_chrom_key(self):
        return SortableByChrom.get_key(self)


# class FusionEvent(SVEvent):
#     def __init__(self, **kwargs):
#         SVEvent.__init__(self, **kwargs)
#         self.is_know_sv = False
#         self.fusion_pair = False


def get_key_or_target_bed_genes(bed_fpath, key_genes_fpath):
    use_custom_panel = False
    key_gene_names_chroms = None
    if bed_fpath:
        key_gene_names_chroms, gene_names_list = get_gene_keys(bed_fpath)
        if key_gene_names_chroms and len(key_gene_names_chroms) < 2000:
            use_custom_panel = True
    if not use_custom_panel:
        key_gene_names = get_key_genes(key_genes_fpath)
        key_gene_names_chroms = [(gn, None) for gn in key_gene_names]
    key_gene_names_chroms = [(gn, c) for (gn, c) in key_gene_names_chroms if (gn and gn != '.')]
    return key_gene_names_chroms, use_custom_panel


def get_key_genes(key_genes_fpath):
    key_genes_fpath = verify_file(key_genes_fpath, is_critical=True, description='820 AZ key genes')
    with open(key_genes_fpath) as f:
        key_gene_names = set([l.strip() for l in f.readlines() if l.strip() != ''])

    return key_gene_names

