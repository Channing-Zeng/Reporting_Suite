from source import info
from source.clinical_reporting.clinical_parser import clinical_sample_info_from_bcbio_structure


def run_clinical_target2wgs(cnf, wgs_bs, target_bs, shared_sample_names, output_dirpath):
    info('Running clinical reporting comparison')

    for sname in shared_sample_names:
        target_sample = next(s for s in target_bs.samples if s.name == sname)
        wgs_sample = next(s for s in wgs_bs.samples if s.name == sname)

        clinsample_target_info = clinical_sample_info_from_bcbio_structure(cnf, target_bs, target_sample)
        clinsample_wgs_info = clinical_sample_info_from_bcbio_structure(cnf, wgs_bs, wgs_sample)

        run_sample_clinreport_target2wgs(clinsample_target_info, clinsample_wgs_info, output_dirpath)


def run_sample_clinreport_target2wgs(clinsample_target_info, clinsample_wgs_info, output_dirpath):
    pass
