import source
from source.logger import info, err, warn, critical
from source.file_utils import verify_file
from source import BaseSample
from source.targetcov import detail_gene_report_baseending


class Sample(source.BaseSample):
    def __init__(self, name, output_dir, *args, **kwargs):
        BaseSample.__init__(self, name, *args, **kwargs)
        self.output_dir = output_dir

        self.targetcov_html_fpath = self._mk_fpath('{output_dir}/{sample}_{name}/{sample}.{name}.html', source.targetseq_name)
        self.targetcov_json_fpath = self._mk_fpath('{output_dir}/{sample}_{name}/{sample}.{name}.json', source.targetseq_name)
        self.targetcov_detailed_tsv = self._mk_fpath('{output_dir}/{sample}_{name}/{sample}.{name}{ending}.txt', source.targetseq_name, ending=detail_gene_report_baseending)
        self.ngscat_html_fpath = self._mk_fpath('{output_dir}/{sample}_{name}/captureQC.html', source.ngscat_name)
        self.qualimap_html_fpath = self._mk_fpath('{output_dir}/{sample}_{name}/qualimapReport.html', source.qualimap_name)

    def _mk_fpath(self, path_template, name, **kwargs):
        return self.make_fpath(path_template, name=name, output_dir=self.output_dir, **kwargs)

    def targetcov_done(self):
        if verify_file(self.targetcov_json_fpath) \
           and verify_file(self.targetcov_html_fpath) \
           and verify_file(self.targetcov_detailed_tsv):
            info(self.targetcov_json_fpath + ', ' + self.targetcov_html_fpath +
                 ', and ' + self.targetcov_detailed_tsv + ' exist.')
            return True
        return False

    def ngscat_done(self):
        if verify_file(self.ngscat_html_fpath):
            info(self.ngscat_html_fpath + ' exists.')
            return True
        return False

    def qualimap_done(self):
        if verify_file(self.qualimap_html_fpath):
            info(self.ngscat_html_fpath + ' exists.')
            return True
        return False


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
