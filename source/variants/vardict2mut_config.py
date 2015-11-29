from os.path import join

min_depth = 0
min_num_reads = 0
min_allele_freq = 0.05
max_rate = 1.0  # if a variant is present in > [double] fraction of samples, it's deemed not a mutation.
min_allele_freq_hotspot = None

is_output_fm = False
reg_exp_sample = None
report_reason = False
suppressors_fpath = None
oncogenes_fpath = None

rule_dirpath = '/ngs/reference_data/genomes/Hsapiens/hg19/variation/rules'
annotation_dirpath = '/ngs/reference_data/genomes/Hsapiens/hg19/variation/cancer_informatics'

filter_common_snp_fpath = join(annotation_dirpath, 'filter_common_snp.txt')
snpeffect_export_polymorphic_fpath = join(annotation_dirpath, 'snpeffect_export_Polymorphic.txt')
filter_common_artifacts_fpath = join(annotation_dirpath, 'filter_common_artifacts.txt')
actionable_hotspot_txt_fpath = join(annotation_dirpath, 'actionable_hotspot.txt')
actionable_txt_fpath = join(annotation_dirpath, 'actionable.txt')
compendia_ms7_hotspot_fpath = join(annotation_dirpath, 'Compendia.MS7.Hotspot.txt')