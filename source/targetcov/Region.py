from collections import defaultdict, OrderedDict
import copy
import math
import os
from os.path import isfile, join
import sys
from source.logger import info, err
from source.ngscat.bed_file import verify_bed

class Region:
    def __init__(self, sample_name=None, gene_name=None, exon_num=None, strand=None, biotype=None,
                 feature=None, extra_fields=list(),
                 chrom=None, start=None, end=None, size=None, min_depth=None,
                 avg_depth=None, std_dev=None, rate_within_normal=None, bases_by_depth=None):

        self.sample_name = sample_name
        self.gene_name = gene_name
        self.feature = feature
        self.extra_fields = extra_fields  # for exons, extra_fields is [Gene, Exon number, Strand]
        self.exon_num = exon_num
        self.strand = strand
        self.biotype = biotype

        self.chrom = chrom
        self.start = start  # int
        self.end = end      # int
        self.size = size    # int
        self.min_depth = min_depth
        self.bases_by_depth = bases_by_depth or defaultdict(int)  # filled in from the "bedcoverage hist" output

        # Calculated once on "sum_up()", when all self.bases_by_depth are there:
        self.avg_depth = avg_depth
        self.std_dev = std_dev
        self.rate_within_normal = rate_within_normal
        self.bases_within_threshs = None    # OrderedDict((depth, 0) for depth in depth_thresholds)
        self.percent_within_threshs = None  # defaultdict(float)

        self.missed_by_db = dict()
        self.var_num = None

    # def get_chrom_num(self):
    #     digits = [c for c in self.chrom if c.isdigit()]
    #     if digits:
    #         return int(''.join(digits))
    #     if 'M' in self.chrom:
    #         return 0
    #     if 'X' in self.chrom:
    #         return 23
    #     if 'Y' in self.chrom:
    #         return 24
    #     else:
    #         return 25
    #
    # @staticmethod
    # def get_order_key(r):
    #     return r.get_chrom_num(), r.get_start(), r.get_end(), r.gene_name

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

    def get_avg_depth(self):
        return self.avg_depth

    def get_std_dev(self):
        return self.std_dev

    def get_bases_by_depth(self):
        return self.bases_by_depth

    def get_size(self):
        if self.size is not None:
            return self.size
        if self.start and self.end:
            return abs(self.end - self.start)
        return None

    def add_bases_for_depth(self, depth, bases):
        assert depth not in self.bases_by_depth, \
            'Duplicated depth ' + str(depth) + ' for the region ' + str(self)
        self.bases_by_depth[depth] += bases

    def __str__(self):
        ts = [self.sample_name, self.chrom, self.start, self.end,
              self.gene_name, self.feature]
        return '"' + '\t'.join(map(str, ts)) + '"'

    def __repr__(self):
        return self.__str__()

    def calc_bases_within_threshs(self, depth_thresholds):
        if self.bases_within_threshs is not None:
            return self.bases_within_threshs

        if self.bases_by_depth is None:
            err('Error: self.bases_by_depth is None for ' + str(self))

        self.bases_within_threshs = OrderedDict((depth, 0) for depth in depth_thresholds)
        for depth, bases in self.bases_by_depth.iteritems():
            for depth_thres in depth_thresholds:
                if depth >= depth_thres:
                    self.bases_within_threshs[depth_thres] += bases
        return self.bases_within_threshs

    def calc_avg_depth(self):
        if self.avg_depth is not None:
            return self.avg_depth

        if self.bases_by_depth:
            depth_sum = sum(
                depth * bases
                for depth, bases
                in self.bases_by_depth.items())
            self.avg_depth = float(depth_sum) / self.get_size() if self.get_size() else None
            return self.avg_depth

    def calc_std_dev(self, avg_depth):
        if avg_depth is None:
            return None

        if self.std_dev is not None:
            return self.std_dev

        if self.bases_by_depth:
            sum_of_sq_var = sum(
                (depth - avg_depth) ** 2 * bases
                for depth, bases
                in self.bases_by_depth.items())
            self.std_dev = math.sqrt(float(sum_of_sq_var) / self.get_size())
            return self.std_dev

    def calc_rate_within_normal(self, avg_depth):
        if avg_depth is None:
            return None

        if self.rate_within_normal is not None:
            return self.rate_within_normal

        if self.bases_by_depth:
            bases_within_normal = sum(
                bases
                for depth, bases
                in self.bases_by_depth.items()
                if math.fabs(avg_depth - depth) < 0.2 * avg_depth)

            self.rate_within_normal = 1.0 * bases_within_normal / self.get_size() \
                if self.get_size() else None
            return self.rate_within_normal

    def sum_up(self, depth_thresholds):
        self.calc_avg_depth()
        self.calc_std_dev(self.avg_depth)
        self.calc_bases_within_threshs(depth_thresholds)
        self.calc_rate_within_normal(self.avg_depth)
        return self.bases_within_threshs,\
               self.bases_by_depth, \
               self.avg_depth, \
               self.std_dev, \
               self.rate_within_normal

    def intersect(self, reg2):
        return self.chrom == reg2.chrom and \
               (reg2.start < self.start < reg2.end or
                self.start < reg2.start < self.end)


class GeneInfo(Region):
    """ Collects assisiated exons and overlapping amlicons.
        - Knows its sample, gene name and chromosome.

        - Stores amplicons and exons ("Region" instances).

        - Supports extending with exons in sorted by starting position order;
          when adding a new exon, recalculates start, end, size and based_by_depth.

    """
    def __init__(self, sample_name, gene_name, chrom=None, strand=None, feature='Whole-Gene', exon_num=None):
        Region.__init__(self, sample_name=sample_name, gene_name=gene_name, exon_num=exon_num, strand=strand,
                              feature=feature, chrom=chrom)
        self.exons = []
        self.amplicons = []
        self.non_overlapping_exons = []
        self.size = 0
        self.min_depth = None

    def get_exons(self):
        return self.exons  # self.subregions_by_feature['Exon']['regions']

    def get_amplicons(self):
        return self.amplicons  # self.subregions_by_feature['Capture']['regions']

    def add_exon(self, exon):  # exons come sorted by start
        if self.exons == []:
            if self.start is None:
                self.start = exon.start
            if self.end is None:
                self.end = exon.end
        else:
            if self.end is None:
                if exon.end > self.end:
                    self.end = exon.end
        self.size += exon.get_size()
        self.exons.append(exon)
        for depth, bases in exon.bases_by_depth.items():
            self.bases_by_depth[depth] += bases

        self.min_depth = min(self.min_depth, exon.min_depth) if self.min_depth else exon.min_depth

    def add_amplicon(self, amplicon):
        # amplicon = copy.copy(amplicon)
        amplicon.gene_name = amplicon.gene_name or self.gene_name
        self.amplicons.append(amplicon)


def proc_regions(regions, fn, *args, **kwargs):
    i = 0
    for region in regions:
        i += 1

        fn(region, *args, **kwargs)

        if i % 10000 == 0:
            info('Processed {0:,} regions.'.format(i))

    if not i % 10000:
        info('Processed {0:,} regions.'.format(i))


def save_regions_to_bed(cnf, regions, f_basename, save_original_fields=False):
    bed_fpath = join(cnf.work_dir, f_basename + '.bed')
    info('Writing regions to ' + bed_fpath)

    if isfile(bed_fpath):
        if cnf.reuse_intermediate:
            verify_bed(bed_fpath, is_critical=True)
            return bed_fpath
        else:
            os.remove(bed_fpath)

    with open(bed_fpath, 'w') as f:
        for r in regions:
            ts = [r.chrom, str(r.get_start()), str(r.get_end()), r.gene_name or '.']
            if save_original_fields:
                ts.extend(r.extra_fields)
            else:
                ts.append(r.feature or '.')
            f.write('\t'.join(ts) + '\n')

                # r.size, r.avg_depth, r.std_dev, r.percent_within_normal
            # f.write('\t'.join([
            #     region.chrom,
            #     str(region.start), str(region.end), str(region.avg_depth)]) + '\n')

    info('Saved to ' + bed_fpath)
    return bed_fpath