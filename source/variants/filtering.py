from collections import OrderedDict
import os
from os.path import basename, join, isfile, dirname, islink, abspath, isdir

from joblib import Parallel, delayed

import source

from source.bcbio_structure import BCBioStructure
from source.calling_process import call
from source.config import defaults
from source.tools_from_cnf import get_script_cmdline
from source.logger import critical, err, warn
from source.variants.vcf_processing import iterate_vcf, get_sample_column_index
from source.file_utils import safe_mkdir, add_suffix, verify_file, open_gzipsafe, symlink_plus
from source.logger import info


def prep_vcf(vcf_fpath, sample_name, caller_name):
    global glob_cnf
    cnf = glob_cnf

    main_sample_index = get_sample_column_index(vcf_fpath, sample_name)

    def fix_fields(rec):
        if rec.FILTER and rec.FILTER != 'PASS':
            return None

        rec.FILTER = 'PASS'
        return rec

    vcf_fpath = iterate_vcf(cnf, vcf_fpath, fix_fields, 'vcf2txt')
    return vcf_fpath


def filter_with_vcf2txt(cnf, bcbio_structure, vcf_fpaths, vcf2txt_out_fpath, sample_by_name, caller_name,
                        sample_min_freq, threads_num):
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
    prep_vcf_fpaths = Parallel(n_jobs=threads_num) \
       (delayed(prep_vcf)(fpath, s, caller_name)
        for fpath, s in zip(vcf_fpaths, sample_by_name.keys()))

    res = run_vcf2txt(cnf, prep_vcf_fpaths, sample_by_name, vcf2txt_out_fpath, sample_min_freq)
    if not res:
        err('vcf2txt run returned non-0')
        return None

    res = run_vardict2mut(cnf, vcf2txt_out_fpath, sample_by_name)
    if not res:
        err('vardict2mut.pl run returned non-0')
        return None

    mut_fpath = res

    # symlinking
    pass_txt_basefname = basename(mut_fpath)
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
        symlink_plus(mut_fpath, pass_txt_fpath_symlink)
    except OSError:
        err('Cannot symlink ' + mut_fpath + ' -> ' + pass_txt_fpath_symlink)

    write_vcfs(cnf, sample_by_name.keys(), bcbio_structure.samples, vcf_fpaths, caller_name, vcf2txt_out_fpath, res, threads_num)

    return res


def run_vardict2mut(cnf, vcf2txt_res_fpath, sample_by_name):
    pick_line = get_script_cmdline(cnf, 'perl', join('VarDict', 'pickLine'), is_critical=True)

    with open(vcf2txt_res_fpath) as f:
        try:
            pass_col_num = f.readline().split().index('PASS') + 1
        except ValueError:
            critical('No PASS column in the vcf2txt result ' + vcf2txt_res_fpath)

    vardict2mut_res_fpath = add_suffix(vcf2txt_res_fpath, source.mut_pass_suffix)

    cmdline = '{pick_line} -l PASS:TRUE -c {pass_col_num} {vcf2txt_res_fpath} | grep -vw dbSNP | ' \
              'grep -v UTR_ | grep -vw SILENT | grep -v intron_variant | grep -v upstream_gene_variant | ' \
              'grep -v downstream_gene_variant | grep -v intergenic_region | grep -v intragenic_variant | ' \
              'grep -v NON_CODING'

    polymorphic_variants = cnf.genomes[sample_by_name.values()[0].genome].polymorphic_variants
    if polymorphic_variants:
        poly_vars = abspath(polymorphic_variants)
        cmdline += ' | {pick_line} -v -i 12:3 -c 14:11 {poly_vars}'

    cmdline = cmdline.format(**locals())
    res = call(cnf, cmdline, vardict2mut_res_fpath, exit_on_error=False)
    if not res:
        return None
    else:
        return res


glob_cnf = None


def postprocess_vcf(sample, caller_name, anno_vcf_fpath, variants, passed_variants, vcf2txt_res_fpath):
    global glob_cnf
    cnf = glob_cnf

    info(sample.name + ', ' + caller_name + ': writing filtered VCFs')

    filt_vcf_fpath = sample.get_filt_vcf_fpath_by_callername(caller_name, gz=False)
    pass_filt_vcf_fpath = sample.get_pass_filt_vcf_fpath_by_callername(caller_name, gz=False)
    safe_mkdir(join(sample.dirpath, BCBioStructure.varfilter_dir))

    with open_gzipsafe(anno_vcf_fpath) as vcf_f, \
         open(filt_vcf_fpath, 'w') as filt_f, \
         open(pass_filt_vcf_fpath, 'w') as pass_f:

        info(sample.name + ', ' + caller_name + ': opened ' + anno_vcf_fpath + ', writing to ' + filt_vcf_fpath + ' and ' + pass_filt_vcf_fpath)

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
                            # warn(chrom + ':' + str(pos) + ' ' + str(alt) + ' for ' + anno_vcf_fpath + ' is not at ' + vcf2txt_res_fpath)
                            ts[6] += 'vcf2txt'
                        elif filter_value == 'TRUE':
                            ts[6] += 'EFFECT'
                        else:
                            ts[6] += filter_value
                    filt_f.write('\t'.join(ts))

    info(sample.name + ', ' + caller_name + ': saved filtered VCFs to ' + filt_vcf_fpath + ' and ' + pass_filt_vcf_fpath)


def write_vcfs(cnf, sample_names, samples, anno_vcf_fpaths, caller_name, vcf2txt_res_fpath, mut_res_fpath, threads_num):
    info('')
    info('-' * 70)
    info('Writing VCFs')

    variants = dict()
    passed_variants = set()

    info('Collecting passed variants...')
    with open(mut_res_fpath) as puckline_res_f:
        for l in puckline_res_f:
            ts = l.split('\t')
            s_name, chrom, pos, alt = ts[0], ts[1], ts[2], ts[5]
            passed_variants.add((s_name, chrom, pos, alt))

    info('Collecting all vcf2txt variants...')
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

    info()
    info('Writing filtered VCFs')
    Parallel(n_jobs=threads_num) \
        (delayed(postprocess_vcf) \
            (next(s for s in samples if s.name == s_name), caller_name, anno_vcf_fpath, variants, passed_variants, vcf2txt_res_fpath)
            for s_name, anno_vcf_fpath in zip(sample_names, anno_vcf_fpaths))


def filter_for_variant_caller(caller, cnf, bcbio_structure):
    info('Running for ' + caller.name)

    all_vcf_by_sample = caller.find_anno_vcf_by_sample()
    if len(all_vcf_by_sample) == 0:
        err('No vcfs for ' + caller.name + '. Skipping.')
        return caller

    def fill_in(batches):
        vcf_by_sample = OrderedDict()
        for b in batches:
            info('Batch ' + b.name)
            if b.normal:
                info('  normal sample: ' + b.normal.name)
                vcf_fpath = all_vcf_by_sample.get(b.normal.name)
                if vcf_fpath:
                    vcf_by_sample[b.normal.name] = vcf_fpath
                    info('  normal VCF: ' + vcf_fpath)
            if len(b.tumor) > 1:
                err('  ERROR: ' + caller.name + ': ' + str(len(b.tumor)) + ' tumor samples (' + ', '.join(t.name for t in b.tumor) + ') for batch ' + b.name)
            if len(b.tumor) > 0:
                info('  tumor sample: ' + b.tumor[0].name)
                vcf_fpath = all_vcf_by_sample.get(b.tumor[0].name)
                if vcf_fpath:
                    vcf_by_sample[b.tumor[0].name] = vcf_fpath
                    info('  tumor VCF: ' + vcf_fpath)
        return vcf_by_sample

    single_vcf_by_sample = fill_in([b for b in bcbio_structure.batches.values() if not b.paired])
    paired_vcf_by_sample = fill_in([b for b in bcbio_structure.batches.values() if b.paired])

    if single_vcf_by_sample:
        info('*' * 70)
        info('Single samples (total ' + str(len(paired_vcf_by_sample)) + '):')

        vcf2txt_fname = source.mut_fname_template.format(caller_name=caller.name)
        if paired_vcf_by_sample:
            vcf2txt_fname = add_suffix(vcf2txt_fname, source.mut_single_suffix)
        vcf2txt_fpath = join(bcbio_structure.var_dirpath, vcf2txt_fname)

        vcf2txt_fpath, mut_fpath = __filter_for_vcfs(cnf, bcbio_structure, caller.name, single_vcf_by_sample, vcf2txt_fpath)
        caller.single_vcf2txt_res_fpath = vcf2txt_fpath
        caller.single_mut_res_fpath = mut_fpath

    if paired_vcf_by_sample:
        info('*' * 70)
        info('Paired samples (total ' + str(len(paired_vcf_by_sample)) + '):')

        vcf2txt_fname = source.mut_fname_template.format(caller_name=caller.name)
        if single_vcf_by_sample:
            vcf2txt_fname = add_suffix(vcf2txt_fname, source.mut_paired_suffix)
        vcf2txt_fpath = join(bcbio_structure.var_dirpath, vcf2txt_fname)

        vcf2txt_fpath, mut_fpath = __filter_for_vcfs(cnf, bcbio_structure, caller.name, paired_vcf_by_sample, vcf2txt_fpath)
        caller.paired_vcf2txt_res_fpath = vcf2txt_fpath
        caller.paired_mut_res_fpath = mut_fpath

    info('-' * 70)
    info()

    return caller


def __filter_for_vcfs(cnf, bcbio_structure, caller_name, vcf_fpaths, vcf2txt_res_fpath):
    threads_num = min(len(vcf_fpaths), cnf.threads)
    info('Number of threads for filtering: ' + str(threads_num))

    sample_names = vcf_fpaths.keys()
    sample_by_name = OrderedDict((sn, next(s for s in bcbio_structure.samples if s.name == sn)) for sn in sample_names)

    info()
    info('-' * 70)
    info('Filtering using vcf2txt...')

    mut_fpath = filter_with_vcf2txt(cnf, bcbio_structure, vcf_fpaths.values(), vcf2txt_res_fpath, sample_by_name, caller_name,
         bcbio_structure.samples[0].min_af, threads_num)
    if not mut_fpath:
        return None, None

    return vcf2txt_res_fpath, mut_fpath


# def combine_mafs(cnf, maf_fpaths, output_basename):
#     output_fpath = output_basename + '.maf'
#     output_pass_fpath = output_basename + '.pass.maf'
#
#     if isfile(output_fpath): os.remove(output_fpath)
#     if isfile(output_pass_fpath): os.remove(output_pass_fpath)
#
#     if not maf_fpaths:
#         warn('No MAFs - no combined MAF will be made.')
#         return None, None
#
#     if not isdir(dirname(output_fpath)): safe_mkdir(dirname(output_fpath))
#
#     with open(output_fpath, 'w') as out, \
#          open(output_pass_fpath, 'w') as out_pass:
#
#         for i, fpath in enumerate(maf_fpaths):
#             with open(fpath) as inp:
#                 for j, line in enumerate(inp):
#                     if i > 0 and j in [0, 1]:
#                         continue
#                     out.write(line)
#                     if '\tInvalid\t' not in line:
#                         out_pass.write(line)
#     return output_fpath, output_pass_fpath


def run_vcf2txt(cnf, vcf_fpaths, sample_by_name, final_maf_fpath, sample_min_freq=None):
    info()
    info('Running VarDict vcf2txt...')

    vcf2txt = get_script_cmdline(cnf, 'perl', join('VarDict', 'vcf2txt.pl'), is_critical=True)

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

