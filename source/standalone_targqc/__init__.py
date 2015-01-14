import source
from source.logger import info, err, warn, critical
from source.file_utils import verify_file
from source import BaseSample
from source.targetcov import detail_gene_report_baseending


class StandaloneSample(source.BaseSample):
    def __init__(self, name, dirpath, *args, **kwargs):
        BaseSample.__init__(self, name, dirpath, '{dirpath}/{sample}_{name}/', *args, **kwargs)


# def make_fpaths_per_sample(samples, output_dir, dir_name, f_template):
#     fpaths_dict = dict()
#     for s in samples:
#         fpath = join(output_dir, s.name + '_' + dir_name, f_template.format(**dict(sample=s.name, dir_name=dir_name)))
#         fpaths_dict[s.name] = fpath
#
#     return fpaths_dict


# def find_fpaths_per_sample(samples, output_dir, dir_name, f_template):
#     fpaths_dict = make_fpaths_per_sample(samples, output_dir, dir_name, f_template)
#     return {sn: fpath for sn, fpath in fpaths_dict.items() if verify_file(fpath)}
