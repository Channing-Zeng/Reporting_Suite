from collections import defaultdict, OrderedDict
import math


class Region:
    def __init__(self, sample=None, chrom=None, start=None, end=None, size=None,
                 gene_name=None, feature=None, extra_fields=list()):
        self.sample = sample
        self.chrom = chrom
        self.start = start
        self.end = end
        self.size = size
        self.gene_name = gene_name
        self.feature = feature
        self.extra_fields = extra_fields
        self.bases_by_depth = defaultdict(int)
        self.subregions = []

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
        return hash((self.sample, self.chrom, self.start, self.end))

    def __str__(self, depth_thresholds=None):
        ts = [self.sample, self.chrom, self.start, self.end, self.feature] + self.extra_fields
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