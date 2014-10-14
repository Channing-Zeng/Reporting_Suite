from collections import defaultdict, OrderedDict
import math
import os
from os.path import isfile, join
import sys
from source.logger import info
from source.ngscat.bed_file import verify_bed


class Region:
    def __init__(self, sample_name=None, chrom=None, start=None, end=None,
                 gene_name=None, feature=None, size=None, avg_depth=None, std_dev=None,
                 percent_within_normal=None, extra_fields=list()):
        self.sample_name = sample_name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.size = size
        self.gene_name = gene_name
        self.feature = feature
        self.avg_depth = avg_depth
        self.std_dev = std_dev
        self.percent_within_normal = percent_within_normal
        self.bases_within_threshs = OrderedDict()

        self.extra_fields = extra_fields
        self.bases_by_depth = defaultdict(int)
        self.subregions = []
        self.missed_by_db = dict()
        self.var_num = None
        self.percent_within_threshs = defaultdict(float)

    def add_subregion(self, subregion):
        self.subregions.append(subregion)
        self.size += subregion.get_size()
        self.end = max(subregion.end, self.end)
        for depth, bases in subregion.bases_by_depth.items():
            self.add_bases_for_depth(depth, bases)

    def get_size(self):
        if self.size is not None:
            return self.size
        if self.start and self.end:
            return abs(self.end - self.start)
        return None

    def add_bases_for_depth(self, depth, bases):
        self.bases_by_depth[depth] += bases

    def key(self):
        return hash((self.sample_name, self.chrom, self.start, self.end))

    def __str__(self, depth_thresholds=None):
        ts = [self.sample_name, self.chrom, self.start, self.end, self.feature] + self.extra_fields
        return '"' + '\t'.join(map(str, ts)) + '"'

    def calc_bases_within_threshs(self, depth_thresholds):
        self.bases_within_threshs = OrderedDict((depth, 0) for depth in depth_thresholds)
        for depth, bases in self.bases_by_depth.iteritems():
            for depth_thres in depth_thresholds:
                if depth >= depth_thres:
                    self.bases_within_threshs[depth_thres] += bases

        return self.bases_within_threshs

    def calc_avg_depth(self):
        depth_sum = sum(
            depth * bases
            for depth, bases
            in self.bases_by_depth.items())
        self.avg_depth = float(depth_sum) / self.get_size()
        return self.avg_depth

    def calc_std_dev(self, avg_depth):
        sum_of_sq_var = sum(
            (depth - avg_depth) ** 2 * bases
            for depth, bases
            in self.bases_by_depth.items())
        self.std_dev = math.sqrt(float(sum_of_sq_var) / self.get_size())
        return self.std_dev

    def calc_percent_within_normal(self, avg_depth):
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


def _proc_regions(regions, fn, *args, **kwargs):
    i = 0
    for region in regions:
        i += 1

        fn(region, *args, **kwargs)

        if i % 10000 == 0:
            info('Processed {0:,} regions.'.format(i))

    if not i % 10000:
        info('Processed {0:,} regions.'.format(i))


def _save_regions_to_bed(cnf, regions, f_basename):
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
            f.write('\t'.join(map(str, [
                r.chrom, r.start, r.end, r.gene_name, r.feature])) + '\n')

                # r.size, r.avg_depth, r.std_dev, r.percent_within_normal
            # f.write('\t'.join([
            #     region.chrom,
            #     str(region.start), str(region.end), str(region.avg_depth)]) + '\n')


    info('Saved to ' + bed_fpath)
    return bed_fpath