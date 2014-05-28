from collections import defaultdict, OrderedDict
import math


class Region():
    def __init__(self, sample=None, chrom=None,
                 start=None, end=None, size=None, extra_fields=list(),
                 feature=None):
        self.sample = sample
        self.chrom = chrom
        self.start = start
        self.end = end
        self.size = size
        self.feature = feature
        self.extra_fields = extra_fields
        self.bases_by_depth = defaultdict(int)

    def get_size(self):
        if self.size:
            return self.size
        if self.start is None or self.end is None:
            return None
        return abs(self.end - self.start)

    def add_bases_for_depth(self, depth, bases):
        self.bases_by_depth[depth] += bases

    def key(self):
        return hash((self.sample, self.chrom, self.start, self.end))

    def __str__(self, depth_thresholds):
        ts = [self.sample, self.chrom, self.start, self.end, self.feature] + self.extra_fields
        return '\t'.join(map(str, ts))

    def bases_within_threshs(self, depth_thresholds):
        bases_within_threshs = OrderedDict((depth, 0) for depth in depth_thresholds)
        for depth, bases in self.bases_by_depth.iteritems():
            for depth_thres in depth_thresholds:
                if depth >= depth_thres:
                    bases_within_threshs[depth_thres] += bases

        return bases_within_threshs

    def avg_depth(self):
        depth_sum = sum(
            depth * bases
            for depth, bases
            in self.bases_by_depth.items())
        avg_depth = float(depth_sum) / self.get_size()
        return avg_depth

    def std_dev(self, avg_depth):
        sum_of_sq_var = sum(
            (depth - avg_depth) ** 2 * bases
            for depth, bases
            in self.bases_by_depth.items())
        std_dev = math.sqrt(float(sum_of_sq_var) / self.get_size())
        return std_dev

    def percent_within_normal(self, avg_depth):
        bases_within_normal = sum(
            bases
            for depth, bases
            in self.bases_by_depth.items()
            if math.fabs(avg_depth - depth) < 0.2 * avg_depth)
        percent_within_normal = 100.0 * bases_within_normal / self.get_size() \
            if self.get_size() else None

        return percent_within_normal

    def sum_up(self, depth_thresholds):
        avg_depth = self.avg_depth()
        return self.bases_within_threshs(depth_thresholds), avg_depth, \
               self.std_dev(avg_depth), self.percent_within_normal(avg_depth)

