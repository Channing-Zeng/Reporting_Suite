from collections import defaultdict, OrderedDict
import copy
import math
import os
from os.path import isfile, join
import sys
from source.logger import info
from source.ngscat.bed_file import verify_bed


class Region:
    def __init__(self, sample_name=None, chrom=None, start=None, end=None,
                 gene_name=None, exon_num=None, strand=None, feature=None, size=None, avg_depth=None,
                 std_dev=None, percent_within_normal=None, extra_fields=list(), bases_by_depth=None):
        self.sample_name = sample_name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.size = size
        self.extra_fields = extra_fields  # for exons, extra_fields is [Gene, Exon number, Strand]
        self.bases_by_depth = bases_by_depth or defaultdict(int)

        self.gene_name = gene_name
        self.exon_num = exon_num
        self.strand = strand

        self.feature = feature
        self.avg_depth = avg_depth
        self.std_dev = std_dev
        self.percent_within_normal = percent_within_normal

        self.missed_by_db = dict()
        self.var_num = None
        self.bases_within_threshs = OrderedDict()
        self.percent_within_threshs = defaultdict(float)

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
        self.bases_by_depth[depth] += bases

    def key(self):
        return hash((self.sample_name, self.chrom, self.start, self.end, self.extra_fields))

    def __str__(self, depth_thresholds=None):
        ts = [self.sample_name, self.chrom, self.start, self.end, self.feature] + self.extra_fields
        return '"' + '\t'.join(map(str, ts)) + '"'

    def calc_bases_within_threshs(self, depth_thresholds):
        if self.bases_by_depth and not self.bases_within_threshs:
            self.bases_within_threshs = OrderedDict((depth, 0) for depth in depth_thresholds)
            for depth, bases in self.bases_by_depth.iteritems():
                for depth_thres in depth_thresholds:
                    if depth >= depth_thres:
                        self.bases_within_threshs[depth_thres] += bases

        return self.bases_within_threshs

    def calc_avg_depth(self):
        if self.bases_by_depth:
            depth_sum = sum(
                depth * bases
                for depth, bases
                in self.bases_by_depth.items())
            self.avg_depth = float(depth_sum) / self.get_size()
        return self.avg_depth

    def calc_std_dev(self, avg_depth):
        if self.bases_by_depth:
            sum_of_sq_var = sum(
                (depth - avg_depth) ** 2 * bases
                for depth, bases
                in self.bases_by_depth.items())
            self.std_dev = math.sqrt(float(sum_of_sq_var) / self.get_size())
        return self.std_dev

    def calc_percent_within_normal(self, avg_depth):
        if self.bases_by_depth:
            bases_within_normal = sum(
                bases
                for depth, bases
                in self.bases_by_depth.items()
                if math.fabs(avg_depth - depth) < 0.2 * avg_depth)

            self.percent_within_normal = 100.0 * bases_within_normal / self.get_size() \
                if self.get_size() else None

        return self.percent_within_normal

    def sum_up(self, depth_thresholds):
        self.calc_bases_within_threshs(depth_thresholds)
        self.calc_avg_depth()
        self.calc_std_dev(self.avg_depth)
        self.calc_percent_within_normal(self.avg_depth)
        return self.bases_within_threshs,\
               self.bases_by_depth, \
               self.avg_depth, \
               self.std_dev, \
               self.percent_within_normal

    def intersect(self, reg2):
        return self.chrom == reg2.chrom and \
               (reg2.start < self.start < reg2.end or
                self.start < reg2.start < self.end)

    def get_order_key(r):
        chr_n = int(''.join(c for c in r.chrom if c.isdigit()))
        return chr_n, r.start, r.end


class GeneInfo:
    """ Collects subregions.
        - Knows it's sample, gene and chromosome.
        - Stores feature-specific subregions, start, end, size and based_by_depth which get recalculated
          on each new subregion.
        - Not a Region, so does not support sum_up, but returns Region object for a given feature.
    """
    def __init__(self, sample_name, gene_name, chrom, feature):
        self.sample_name = sample_name
        self.gene_name = gene_name
        self.feature = feature
        self.chrom = chrom

        self.subregions_by_feature = dict((f,
             dict(
                 regions=[],
                 start=None,
                 end=None,
                 size=0,
                 bases_by_depth=defaultdict(int)))
             for f in ['Exon', 'Amplicon'])

    def get_exons(self):
        return self.subregions_by_feature['Exon']['regions']

    def get_amplicons(self):
        return self.subregions_by_feature['Amplicon']['regions']

    def get_start(self):
        return self._get_exon_value('start')

    def get_end(self):
        return self._get_exon_value('end')

    def get_size(self):
        return self._get_exon_value('size')

    def get_bases_by_depth(self):
        return self._get_exon_value('bases_by_depth')

    def _get_exon_value(self, key):
        return self.subregions_by_feature['Exon'][key]

    def _add_subregion(self, region, feature):
        self.subregions_by_feature[feature]['regions'].append(region)

        self.subregions_by_feature[feature]['start'] = \
            min(self.subregions_by_feature[feature].get('start') or region.start, region.start)

        self.subregions_by_feature[feature]['end'] = \
            max(self.subregions_by_feature[feature].get('end') or region.end, region.end)

        self.subregions_by_feature[feature]['size'] = \
            self.subregions_by_feature[feature]['size'] + region.size

        for d, bs in region.bases_by_depth.items():
            self.subregions_by_feature[feature]['bases_by_depth'][d] += bs

    def add_exon(self, exon):
        self._add_subregion(exon, 'Exon')

    def add_amplicon(self, amplicon):
        amplicon = copy.copy(amplicon)
        amplicon.gene_name = amplicon.gene_name or self.gene_name
        self._add_subregion(amplicon, 'Amplicon')

    def get_summary_region(self, regions, feature, depth_thresholds):
        region = Region(
            sample_name=self.sample_name,
            chrom=self.chrom,
            feature='Gene-' + feature,
            gene_name=self.gene_name,
            start=self.subregions_by_feature[feature]['start'],
            end=self.subregions_by_feature[feature]['end'],
            size=self.subregions_by_feature[feature]['size'],
            bases_by_depth=self.subregions_by_feature[feature]['bases_by_depth'])
        region.sum_up(depth_thresholds)
        return region

    # def add_subregion(self, subregion):
    #     self.size += subregion.get_size()
    #     self.end = max(subregion.end, self.end)
    #     self.subregions.append(subregion)
        # for depth, bases in subregion.bases_by_depth.items():
        #     self.add_bases_for_depth(depth, bases)


def _proc_regions(regions, fn, *args, **kwargs):
    i = 0
    for region in regions:
        i += 1

        fn(region, *args, **kwargs)

        if i % 10000 == 0:
            info('Processed {0:,} regions.'.format(i))

    if not i % 10000:
        info('Processed {0:,} regions.'.format(i))


def save_regions_to_bed(cnf, regions, f_basename, save_feature=True):
    bed_fpath = join(cnf.work_dir, f_basename + '.bed')
    info()
    info('Saving regions to ' + bed_fpath)

    if isfile(bed_fpath):
        if cnf.reuse_intermediate:
            if not verify_bed(bed_fpath):
                sys.exit(1)
            return bed_fpath
        else:
            os.remove(bed_fpath)

    with open(bed_fpath, 'w') as f:
        for r in regions:
            ts = [r.chrom, str(r.get_start()), str(r.get_end()), r.gene_name or '.']
            if save_feature:
                ts.append(r.feature or '.')
            f.write('\t'.join(ts) + '\n')

                # r.size, r.avg_depth, r.std_dev, r.percent_within_normal
            # f.write('\t'.join([
            #     region.chrom,
            #     str(region.start), str(region.end), str(region.avg_depth)]) + '\n')

    info('Saved to ' + bed_fpath)
    return bed_fpath