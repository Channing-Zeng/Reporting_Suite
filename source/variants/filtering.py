from collections import OrderedDict
from genericpath import isdir, exists
import os
import shutil
import pickle
import operator

from os.path import basename, join, isfile, dirname, splitext, islink, pardir, abspath
from joblib import Parallel, delayed
##from memory_profiler import profile
import sys

from source.bcbio_structure import BCBioStructure
from source.calling_process import call
from source.config import defaults
from source.tools_from_cnf import get_script_cmdline, get_system_path
from source.variants.Effect import Effect
from source.logger import step_greetings, info, critical, err, warn
from source.variants.anno import _snpsift_annotate
from source.variants.vcf_processing import iterate_vcf, vcf_one_per_line, \
    get_sample_column_index, bgzip_and_tabix, vcf_merge, leave_main_sample
from source.utils import mean
from source.file_utils import safe_mkdir, add_suffix, verify_file, open_gzipsafe
from source.variants.tsv import make_tsv
from source.variants.vcf_processing import remove_rejected, vcf_is_empty, igvtools_index
from source.logger import info


def prep_vcf(vcf_fpath, sample_name, caller_name):
    global glob_cnf
    cnf = glob_cnf

    main_sample_index = get_sample_column_index(vcf_fpath, sample_name)

    # vcf_fpath = leave_main_sample(cnf, vcf_fpath, sample_name)

    def fix_fields(rec):
        if rec.FILTER and rec.FILTER != 'PASS':
            return None

        # ads = rec.get_sample_val('AD', main_sample_index)
        # if ads is None:
        #     aos = rec.get_sample_val('AO', main_sample_index)
        #     ro = rec.get_info_val('RO', main_sample_index)
        #     missing = [n for f, n in [(aos, 'AO'), (ro, 'RO')] if f is None]
        #     if missing:
        #         err('No AD or ' + ','.join(missing) + ' for ' + rec.get_variant() + ', ' + vcf_fpath)
        #         return None
        #     try:
        #         rec.INFO['AD'] = [ro] + aos
        #     except TypeError:
        #         rec.INFO['AD'] = [ro] + [aos]

            # rec.FORMAT += 'AO'
            # call_data = rec.genotype(sample_name).data._asdict()
            # try:
            #     call_data.AD = [ro] + aos
            # except TypeError:
            #     call_data.AD = [ro] + [aos]

            # af = float(t_alt_count) / dp
            # rec.INFO['AF'] = af

        # for f in ['QUAL', 'PMEAN', 'MQ', 'SN', 'VD']:
        #     val = rec.get_info_val(f)
        #     if val is None:
        #         val = rec.get_sample_val(f, main_sample_index)
        #         if val is None:
        #             rec.INFO[f] = 999999999

        # if rec.get_val('PSTD', main_sample_index) is None:
        #     rec.INFO['PSTD'] = '1.0'

        return rec

    vcf_fpath = iterate_vcf(cnf, vcf_fpath, fix_fields, 'vcf2txt')
    return vcf_fpath


def filter_with_vcf2txt(cnf, bcbio_structure, vcf_fpaths, sample_by_name, caller, sample_min_freq):
    # info()
    # info('Keeping only first sample info in VCFs...')
    # vcf_fpaths = Parallel(n_jobs=n_threads)(delayed(leave_main_sample)(cnf, fpath, s) for fpath, s in zip(vcf_fpaths, sample_names))

    # info()
    # info('Indexing with tabix ')
    # vcf_fpaths = Parallel(n_jobs=n_threads)(delayed(tabix_vcf)(cnf, vcf_fpath) for vcf_fpath in vcf_fpaths)

    safe_mkdir(bcbio_structure.var_dirpath)

    # combined_vcf_fpath = join(bcbio_structure.var_dirpath, caller.name + '.vardict.vcf')
    #
    # info()
    # info('Merging VCFs...')
    # vcf_merge(cnf, vcf_fpaths, combined_vcf_fpath)

    global glob_cnf
    glob_cnf = cnf

    info()
    info('Preparing VCFs for vcf2txt')
    prep_vcf_fpaths = Parallel(
        n_jobs=cnf.threads)(delayed(prep_vcf)(fpath, s, caller.name)
        for fpath, s in zip(vcf_fpaths, sample_by_name.keys()))

    caller.vcf2txt_res_fpath = join(bcbio_structure.var_dirpath, caller.name + '.txt')
    res = run_vcf2txt(cnf, prep_vcf_fpaths, sample_by_name, caller.vcf2txt_res_fpath, sample_min_freq)
    if not res:
        err('vcf2txt run returned non-0')
        return None

    res = run_pickline(cnf, caller, caller.vcf2txt_res_fpath, sample_by_name)
    if not res:
        err('pickLine run returned non-0')
        return None

    write_vcfs(cnf, sample_by_name.keys(), vcf_fpaths, caller, caller.vcf2txt_res_fpath, caller.pickline_res_fpath)

    return res


def run_pickline(cnf, caller, vcf2txt_res_fpath, sample_by_name):
    pick_line = get_script_cmdline(cnf, 'perl', join('VarDict', 'pickLine'))
    if not pick_line:
        sys.exit(1)
        return None

    with open(vcf2txt_res_fpath) as f:
        try:
            pass_col_num = f.readline().split().index('PASS') + 1
        except ValueError:
            critical('No PASS column in the vcf2txt result ' + vcf2txt_res_fpath)

    caller.pickline_res_fpath = add_suffix(vcf2txt_res_fpath, 'PASS')

    cmdline = '{pick_line} -l PASS:TRUE -c {pass_col_num} {vcf2txt_res_fpath} | grep -vw dbSNP | ' \
              'grep -v UTR_ | grep -vw SILENT | grep -v intron_variant | grep -v upstream_gene_variant | ' \
              'grep -v downstream_gene_variant | grep -v intergenic_region | grep -v intragenic_variant | ' \
              'grep -v NON_CODING'

    polymorphic_variants = cnf.genomes[sample_by_name.values()[0].genome].polymorphic_variants
    if polymorphic_variants:
        poly_vars = abspath(polymorphic_variants)
        cmdline += ' | {pick_line} -v -i 12:3 -c 14:11 {poly_vars}'
    cmdline = cmdline.format(**locals())
    res = call(cnf, cmdline, caller.pickline_res_fpath, exit_on_error=False)
    if not res:
        return None
    else:
        return res


glob_cnf = None


def postprocess_vcf(sample, caller_name, anno_vcf_fpath, variants, passed_variants, vcf2txt_res_fpath):
    global glob_cnf
    cnf = glob_cnf

    filt_vcf_fpath = sample.get_filt_vcf_fpath_by_callername(caller_name, gz=False)
    pass_filt_vcf_fpath = sample.get_pass_filt_vcf_fpath_by_callername(caller_name, gz=False)
    safe_mkdir(join(sample.dirpath, BCBioStructure.varfilter_dir))

    with open_gzipsafe(anno_vcf_fpath) as vcf_f, \
         open(filt_vcf_fpath, 'w') as filt_f, \
         open(pass_filt_vcf_fpath, 'w') as pass_f:
        for l in vcf_f:
            if l.startswith('#'):
                filt_f.write(l)
                pass_f.write(l)
            else:
                ts = l.split('\t')
                chrom, pos, alt = ts[0], ts[1], ts[4]
                if (sample.name, chrom, pos, alt) in passed_variants:
                    ts[6] = 'PASS'
                    filt_f.write('\t'.join(ts))
                    pass_f.write('\t'.join(ts))
                else:
                    if ts[6] in ['', '.', 'PASS']:
                        ts[6] = ''
                        filter_value = variants.get((sample.name, chrom, pos, alt))
                        if filter_value is None:
                            warn(chrom + ':' + str(pos) + ' ' + str(alt) + ' for ' + anno_vcf_fpath + ' is not at ' + vcf2txt_res_fpath)
                            ts[6] += 'vcf2txt'
                        elif filter_value == 'TRUE':
                            ts[6] += 'EFFECT'
                        else:
                            ts[6] += filter_value
                    filt_f.write('\t'.join(ts))

    # Indexing
    info()
    info('Indexing')
    igvtools_index(cnf, pass_filt_vcf_fpath)
    igvtools_index(cnf, filt_vcf_fpath)
    bgzip_and_tabix(cnf, filt_vcf_fpath)


def write_vcfs(cnf, sample_names, anno_vcf_fpaths, caller, vcf2txt_res_fpath, pickline_res_fpath):
    variants = dict()
    passed_variants = set()

    with open(pickline_res_fpath) as puckline_res_f:
        for l in puckline_res_f:
            ts = l.split('\t')
            s_name, chrom, pos, alt = ts[0], ts[1], ts[2], ts[5]
            passed_variants.add((s_name, chrom, pos, alt))

    with open(vcf2txt_res_fpath) as vcf2txt_f:
        pass_col = None
        for l in vcf2txt_f:
            if l.startswith('Sample'):
                pass_col = l.split('\t').index('PASS')
            else:
                ts = l.split('\t')
                s_name, chrom, pos, alt = ts[0], ts[1], ts[2], ts[5]
                filt = ts[pass_col]
                variants[(s_name, chrom, pos, alt)] = filt

    Parallel(n_jobs=cnf.threads) \
        (delayed(postprocess_vcf) \
            (next(s for s in caller.samples if s.name == s_name), caller.name, anno_vcf_fpath, variants, passed_variants, vcf2txt_res_fpath)
            for s_name, anno_vcf_fpath in zip(sample_names, anno_vcf_fpaths))


def filter_for_variant_caller(caller, cnf, bcbio_structure):
    info('Running for ' + caller.name)

    vcf_by_sample = caller.find_anno_vcf_by_sample()
    # only those samples that were found

    if len(vcf_by_sample) == 0:
        err('No vcfs for ' + caller.name + '. Skipping.')
        return caller

    info('*' * 70)
    info('Samples (total ' + str(len(vcf_by_sample)) + '):')

    # Set up cnfs for each sample
    # global cnfs_for_sample_names
    # for sample in caller.samples:
    #     info(sample.name)
    #     cnf_copy = cnf.copy()
    #     cnf_copy['name'] = sample.name
    #     if cnf_copy['variant_filtering'].get('min_freq') is not None:
    #         sample.min_af = cnf_copy['variant_filtering']['min_freq']
    #     elif sample.min_af is None:
    #         sample.min_af = 0
    #     cnfs_for_sample_names[sample.name] = cnf_copy

    info('Number of threads for filtering: ' + str(cnf.threads))

    vcf_fpaths = vcf_by_sample.values()
    sample_names = vcf_by_sample.keys()
    sample_by_name = OrderedDict((sn, next(s for s in caller.samples if s.name == sn)) for sn in sample_names)

    # global filtering
    # filtering = Filtering(cnf, bcbio_structure, caller)

    # info('Filtering by impact')
    # vcf_fpaths = Parallel(n_jobs=n_threads)(delayed(impact_round)(vcf_fpath, s) for vcf_fpath, s in zip(vcf_fpaths, sample_names))

    info()
    info('-' * 70)
    info('Filtering using vcf2txt...')
    res = filter_with_vcf2txt(cnf, bcbio_structure, vcf_fpaths, sample_by_name, caller, caller.samples[0].min_af)
    if not res:
        return None

    # symlinking
    pass_txt_basefname = basename(caller.pickline_res_fpath)
    pass_txt_fpath_symlink = join(bcbio_structure.date_dirpath, pass_txt_basefname)

    if islink(pass_txt_fpath_symlink):
        try:
            os.unlink(pass_txt_fpath_symlink)
        except OSError:
            pass
    if isfile(pass_txt_fpath_symlink):
        try:
            os.remove(pass_txt_fpath_symlink)
        except OSError:
            pass
    try:
        os.symlink(caller.pickline_res_fpath, pass_txt_fpath_symlink)
    except OSError:
        err('Cannot symlink ' + caller.pickline_res_fpath + ' -> ' + pass_txt_fpath_symlink)

    for sample in caller.samples:
        filt_vcf = sample.find_filt_vcf_by_callername(caller.name)
        if filt_vcf:
            for link in [join(sample.dirpath, basename(filt_vcf)),
                         join(bcbio_structure.var_dirpath, basename(filt_vcf))]:
                if islink(link):
                    try:
                        os.unlink(link)
                    except OSError:
                        pass
                if isfile(link):
                    try:
                        os.remove(link)
                    except OSError:
                        pass
                try:
                    os.symlink(filt_vcf, link)
                except OSError:
                    err('Cannot symlink ' + filt_vcf + ' -> ' + link)

        BCBioStructure.move_vcfs_to_var(sample)

    info('-' * 70)
    info()

    return caller


def combine_mafs(cnf, maf_fpaths, output_basename):
    output_fpath = output_basename + '.maf'
    output_pass_fpath = output_basename + '.pass.maf'

    if isfile(output_fpath): os.remove(output_fpath)
    if isfile(output_pass_fpath): os.remove(output_pass_fpath)

    if not maf_fpaths:
        warn('No MAFs - no combined MAF will be made.')
        return None, None

    if not isdir(dirname(output_fpath)): safe_mkdir(dirname(output_fpath))

    with open(output_fpath, 'w') as out, \
         open(output_pass_fpath, 'w') as out_pass:

        for i, fpath in enumerate(maf_fpaths):
            with open(fpath) as inp:
                for j, line in enumerate(inp):
                    if i > 0 and j in [0, 1]:
                        continue
                    out.write(line)
                    if '\tInvalid\t' not in line:
                        out_pass.write(line)
    return output_fpath, output_pass_fpath


def run_vcf2txt(cnf, vcf_fpaths, sample_by_name, final_maf_fpath, sample_min_freq=None):
    info()
    info('Running VarDict vcf2txt...')

    vcf2txt = get_script_cmdline(cnf, 'perl', join('VarDict', 'vcf2txt.pl'))

    if not vcf2txt:
        sys.exit(1)

    c = cnf.variant_filtering
    min_freq = cnf.min_freq or c.min_freq or sample_min_freq or defaults.default_min_freq

    cmdline = '{vcf2txt} ' \
        '-f {min_freq} -n {c.sample_cnt} -F {c.ave_freq} -p {c.min_p_mean} -q {c.min_q_mean} ' \
        '-r {c.fraction} -R {c.max_ratio} -P {c.filt_p_mean} -Q {c.filt_q_mean} -D {c.filt_depth} ' \
        '-M {c.min_mq} -V {c.min_vd} -G {c.maf} -o {c.signal_noise} '.format(**locals())

    if c.bias:
        cmdline += ' -b '

    if c.count_undetermined is False:
        cmdline += ' -u '

    dbsnp_multi_mafs = cnf.genomes[sample_by_name.values()[0].genome].dbsnp_multi_mafs
    if dbsnp_multi_mafs and verify_file(dbsnp_multi_mafs):
        cmdline += ' -A ' + dbsnp_multi_mafs

    cmdline += ' ' + ' '.join(vcf_fpaths)

    res = call(cnf, cmdline, final_maf_fpath, exit_on_error=False)
    return res

