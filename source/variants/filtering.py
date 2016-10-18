from collections import OrderedDict, defaultdict
import os
from os.path import basename, join, isfile, islink, splitext, isdir, dirname
from random import random
from time import sleep
import traceback
import shutil

from ext_modules.joblib import Parallel, delayed
import source
from source.calling_process import call, call_check_output
from source.config import defaults
from source.profiling import fn_timer
from source.tools_from_cnf import get_script_cmdline, get_java_tool_cmdline
from source.logger import err, warn, send_email, critical
from source.utils import is_us, OrderedDefaultDict
from source.variants.tsv import make_tsv
from source.variants.vcf_processing import bgzip_and_tabix, verify_vcf
from source.file_utils import safe_mkdir, add_suffix, verify_file, open_gzipsafe, \
    symlink_plus, file_transaction, num_lines
from source.logger import info
from source.webserver.exposing import convert_gpfs_path_to_url


def combine_vcfs(cnf, vcf_fpath_by_sname, combined_vcf_fpath, additional_parameters=''):
    gatk = get_java_tool_cmdline(cnf, 'gatk')
    if not gatk:
        info('GATK is not found, skipping merging VCFs')
        return None

    cmdl = '{gatk} -T CombineVariants -R {cnf.genome.seq} {additional_parameters}'.format(**locals())
    for s_name, vcf_fpath in vcf_fpath_by_sname.items():
        if vcf_fpath:
            cmdl += ' --variant:' + s_name + ' ' + vcf_fpath
    if ' --variant:' not in cmdl:
        err('No VCFs to combine')
        return None

    if cnf.reuse_intermediate and isfile(combined_vcf_fpath + '.gz') and verify_vcf(combined_vcf_fpath + '.gz'):
        info(combined_vcf_fpath + '.gz exists, reusing')
        return combined_vcf_fpath + '.gz'

    cmdl += ' -o ' + combined_vcf_fpath
    res = call(cnf, cmdl, output_fpath=combined_vcf_fpath, stdout_to_outputfile=False, exit_on_error=False)
    if res:
        info('Joined VCFs, saved into ' + combined_vcf_fpath)
        if isfile(combined_vcf_fpath + '.tx.idx'):
            try:
                os.remove(combined_vcf_fpath + '.tx.idx')
            except OSError:
                err(traceback.format_exc())
                info()
        return bgzip_and_tabix(cnf, combined_vcf_fpath)
    else:
        warn('Could not join VCFs')
        return None


def index_vcf(cnf, sample_name, filt_vcf_fpath, caller_name=None):
    if cnf is None:
        global glob_cnf
        cnf = glob_cnf

    info()
    info(sample_name + ((', ' + caller_name) if caller_name else '') + ': indexing')

    # for fpath in [pass_vcf_fpath, filt_vcf_fpath]:
    #     if not cnf.reuse_intermediate and not verify_file(fpath, silent=True):
    #         err(fpath + ' does not exist - cannot IGV index')
    #     else:
    #         if cnf.reuse_intermediate and verify_file(fpath + '.idx', silent=True):
    #             info('Reusing existing ' + fpath + '.idx')
    #         else:
    #             igvtools_index(cnf, fpath)

    if not cnf.reuse_intermediate and not verify_file(filt_vcf_fpath, silent=True):
        err(filt_vcf_fpath + ' does not exist - cannot gzip and tabix')
    else:
        if cnf.reuse_intermediate and verify_file(filt_vcf_fpath + '.gz', silent=True) \
                and verify_file(filt_vcf_fpath + '.gz.tbi', silent=True):
            info(filt_vcf_fpath + '.gz and .gz.tbi exist; reusing')
        else:
            bgzip_and_tabix(cnf, filt_vcf_fpath)


def run_vcf2txt_vardict2mut_for_samples(
        cnf, var_samples, output_dirpath, vcf2txt_out_fpath,
        caller_name=None, threads_num=1):

    threads_num = min(len(var_samples), cnf.threads)
    info('Number of threads for filtering: ' + str(threads_num))

    safe_mkdir(output_dirpath)

    vcf_fpath_by_sample = {s.name: s.anno_vcf_fpath for s in var_samples}
    res = run_vcf2txt(cnf, vcf_fpath_by_sample, vcf2txt_out_fpath)
    if not res:
        err('vcf2txt run returned non-0')
        return None

    # vardict2mut_py = get_script_cmdline(cnf, 'python', join('scripts', 'post', 'vardict2mut.py'))
    # if not vardict2mut_py:
    #     critical('vardict2mut_py not found')

    info('Running vardict2mut')
    res = run_vardict2mut(cnf, vcf2txt_out_fpath, add_suffix(vcf2txt_out_fpath, source.mut_pass_suffix))
    if not res:
        critical('vardict2mut.py run returned non-0')
    mut_fpath = res
    mut_fpath = convert_gpfs_path_to_url(mut_fpath)
    info()

    info('Done filtering with vcf2txt/vardict2mut, saved to ' + str(mut_fpath))
    return mut_fpath


def check_filtering_results(fpath):
    if not isfile(fpath):
        return False

    with open(fpath) as f:
        l = next(f, None)
        if not l:
            return False
        if 'CLN_GENE' not in l.split('\t'):
            info(fpath + ' exists, but CLN_GENE columns not found. Removing ' + fpath)
            os.remove(fpath)
            return False
        else:
            return True


@fn_timer
def run_vardict2mut(cnf, vcf2txt_res_fpath, vardict2mut_res_fpath=None, vardict2mut_executable=None):
    cmdline = None
    if vardict2mut_res_fpath is None:
        vardict2mut_res_fpath = add_suffix(vcf2txt_res_fpath, source.mut_pass_suffix)
    vardict2mut_reject_fpath = add_suffix(vcf2txt_res_fpath, source.mut_reject_suffix)

    check_filtering_results(vardict2mut_res_fpath)

    if not vardict2mut_executable:
        # vardict2mut_executable = get_script_cmdline(cnf, 'python', join('scripts', 'post', 'vardict2mut.py'))
        vardict2mut_executable = 'vardict2mut'

    c = cnf.variant_filtering

    cmdline = '{vardict2mut_executable} {vcf2txt_res_fpath} '
    if vardict2mut_executable.endswith('.pl'):
        cmdline += ' --report_reason '
        if c.min_hotspot_freq is not None and c.min_hotspot_freq != 'default':
            cmdline += ' -F ' + str(c.min_hotspot_freq)
        if c.max_ratio_vardict2mut is not None:
            cmdline += ' -R ' + str(c.max_ratio_vardict2mut)
        if cnf.genome.filter_common_snp: cmdline += ' --filter_common_snp {cnf.genome.filter_common_snp} '
        if cnf.genome.filter_common_artifacts: cmdline += ' --filter_common_artifacts {cnf.genome.filter_common_artifacts} '
        if cnf.genome.actionable: cmdline += ' --actionable {cnf.genome.actionable} '
        if cnf.genome.compendia_ms7_hotspot: cmdline += ' --compendia_ms7_hotspot {cnf.genome.compendia_ms7_hotspot} '
        if cnf.snpeffect_export_polymorphic: cmdline += ' --snpeffect_export_polymorphic {cnf.snpeffect_export_polymorphic} '
        if cnf.actionable_hotspot: cmdline += ' --actionable_hotspot {cnf.actionable_hotspot} '
        if cnf.ruledir: cmdline += ' --ruledir {cnf.ruledir} '
        cmdline = cmdline.format(**locals())
        res = call(cnf, cmdline, vardict2mut_res_fpath, exit_on_error=False)

    else:
        filt_yaml_fpath = join(cnf.work_dir, 'filt_cnf.yaml')
        info('Writing filtering yaml into ' + filt_yaml_fpath)
        with file_transaction(cnf.work_dir, filt_yaml_fpath) as tx, open(filt_yaml_fpath, 'w') as out:
            with open(cnf.run_cnf) as run_cnf:
                lines = []
                met_variant_filtering = False
                for l in run_cnf:
                    if l.startswith('variant_filtering:'):
                        met_variant_filtering = True
                        continue
                    if met_variant_filtering:
                        if l.startswith(' '):
                            out.write(l.lstrip())
                        else:
                            break

        cmdline += ' --filt-cnf ' + filt_yaml_fpath
        cmdline += ' --work-dir ' + cnf.work_dir
        cmdline += (' --debug ' if cnf.debug else '')
        cmdline += ' --genome ' + cnf.genome.name
        cmdline += ' -o ' + vardict2mut_res_fpath
        cmdline += ' --o-reject ' + vardict2mut_reject_fpath

        if cnf.cohort_freqs_fpath:
            cmdline += ' --cohort-freqs ' + cnf.cohort_freqs_fpath

        cmdline = cmdline.format(**locals())
        res = call(cnf, cmdline, output_fpath=vardict2mut_res_fpath, stdout_to_outputfile=False)

    if not res:
        return None
    else:
        return res


glob_cnf = None


def postprocess_vcf(cnf,
        work_dir, var_sample, caller_name,
        variants, mutations, vcf2txt_res_fpath):
    if cnf is None:
        global glob_cnf
        cnf = glob_cnf

    info(var_sample.name + ((', ' + caller_name) if caller_name else '') + ': writing filtered VCFs')

    filter_values = set(variants.values())

    # Saving .anno.filt.vcf.gz and .anno.filt.pass.vcf
    ungz, gz = None, None
    if var_sample.filt_vcf_fpath.endswith('.gz'):
        ungz = splitext(var_sample.filt_vcf_fpath)[0]
        gz = var_sample.filt_vcf_fpath
    else:
        ungz = var_sample.filt_vcf_fpath
        gz = var_sample.filt_vcf_fpath + '.gz'
    if not var_sample.filt_tsv_fpath:
        var_sample.filt_tsv_fpath = splitext(ungz)[0] + '.tsv'

    if cnf.reuse_intermediate \
            and verify_file(var_sample.filt_vcf_fpath, silent=True) \
            and verify_file(var_sample.pass_filt_vcf_fpath, silent=True) \
            and verify_file(var_sample.filt_tsv_fpath, silent=True):
        info(var_sample.filt_vcf_fpath + ' and ' + var_sample.pass_filt_vcf_fpath + ' exist; reusing.')

    else:
        safe_mkdir(dirname(var_sample.filt_vcf_fpath))
        safe_mkdir(dirname(var_sample.pass_filt_vcf_fpath))

        with open_gzipsafe(var_sample.anno_vcf_fpath) as vcf_f, \
             file_transaction(work_dir, ungz) as filt_tx, \
             file_transaction(work_dir, var_sample.pass_filt_vcf_fpath) as pass_tx:
            with open(filt_tx, 'w') as filt_f, open(pass_tx, 'w') as pass_f:
                info(var_sample.name + ((', ' + caller_name) if caller_name else '') + ': opened ' +
                     var_sample.anno_vcf_fpath + ', writing to ' +
                     ungz + ' and ' + var_sample.pass_filt_vcf_fpath)

                for l in vcf_f:
                    if l.startswith('#'):
                        if l.startswith('#CHROM'):
                            filt_f.write('##FILTER=<ID=vcf2txt,Description="Hard-filtered by vcf2txt.pl">\n')
                            filt_f.write('##FILTER=<ID=vardict2mut,Description="Hard-filtered by vardict2mut.pl">\n')
                            for filt_val in filter_values:
                                if filt_val != 'PASS':
                                    filt_f.write('##FILTER=<ID=' + filt_val + ',Description="">\n')
                        filt_f.write(l)
                        pass_f.write(l)
                    else:
                        ts = l.split('\t')
                        chrom, pos, alt = ts[0], ts[1], ts[4]
                        if (chrom, pos, alt) in mutations:
                            ts[6] = 'PASS'
                            filt_f.write('\t'.join(ts))
                            pass_f.write('\t'.join(ts))
                        else:
                            if ts[6] in ['', '.', 'PASS']:
                                ts[6] = ''
                                filter_value = variants.get((chrom, pos, alt))
                                if filter_value is None:
                                    ts[6] += 'vcf2txt'
                                elif filter_value == 'TRUE':
                                    ts[6] += 'vardict2mut'
                                else:
                                    ts[6] += filter_value
                            filt_f.write('\t'.join(ts))

        info(var_sample.name + ((', ' + caller_name) if caller_name else '') + ': saved filtered VCFs to ' +
             ungz + ' and ' + var_sample.pass_filt_vcf_fpath)

        if False:
            info()
            info(var_sample.name + ((', ' + caller_name) if caller_name else '') + ': writing filtered TSVs')
            # Converting to TSV - saving .anno.filt.tsv
            if 'tsv_fields' in cnf.annotation and cnf.tsv:
                tmp_tsv_fpath = make_tsv(cnf, ungz, var_sample.name)
                if not tmp_tsv_fpath:
                    err('TSV convertion didn\'t work')
                else:
                    if isfile(var_sample.filt_tsv_fpath):
                        os.remove(var_sample.filt_tsv_fpath)
                    shutil.copy(tmp_tsv_fpath, var_sample.filt_tsv_fpath)

                info(var_sample.name + ((', ' + caller_name) if caller_name else '') +
                     ': saved filtered TSV to ' + var_sample.filt_tsv_fpath)

    info('Done postprocessing filtered VCF.')
    return ungz


def write_vcf(cnf, sample, output_dirpath, caller_name, vcf2txt_res_fpath, mut_res_fpath):
    info('')
    info('-' * 70)
    info('Writing VCF')

    variants_dict = dict()
    mutations = set()

    info('Collecting passed variants...')
    with open(mut_res_fpath) as fh:
        for l in fh:
            if l.strip():
                ts = l.split('\t')
                s_name, chrom, pos, alt = ts[0], ts[1], ts[2], ts[5]
                mutations.add((chrom, pos, alt))

    info('Collecting all vcf2txt variants...')
    with open(vcf2txt_res_fpath) as vcf2txt_f:
        pass_col = None
        for l in vcf2txt_f:
            if l.strip():
                if l.startswith('Sample'):
                    pass_col = l.split('\t').index('PASS')
                else:
                    ts = l.split('\t')
                    s_name, chrom, pos, alt = ts[0], ts[1], ts[2], ts[5]
                    filt = ts[pass_col]
                    variants_dict[(chrom, pos, alt)] = filt

    info()
    info('Writing filtered VCFs')
    return postprocess_vcf(cnf, cnf.work_dir, sample, caller_name, variants_dict, mutations, vcf2txt_res_fpath)


def write_vcfs(cnf, var_samples, output_dirpath,
               caller_name, vcf2txt_res_fpath, mut_res_fpath, threads_num):
    info('')
    info('-' * 70)
    info('Writing VCFs')

    variants_by_sample = defaultdict(dict)
    mutations_by_sample = defaultdict(set)

    info('Collecting passed variants...')
    with open(mut_res_fpath) as fh:
        for l in fh:
            ts = l.split('\t')
            s_name, chrom, pos, alt = ts[0], ts[1], ts[2], ts[5]
            mutations_by_sample[s_name].add((chrom, pos, alt))

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
                variants_by_sample[s_name][(chrom, pos, alt)] = filt

    info()

    info('Writing filtered VCFs in ' + str(threads_num) + ' threads')
    try:
        Parallel(n_jobs=threads_num) \
            (delayed(postprocess_vcf) \
                (None,
                 cnf.work_dir,
                 var_sample,
                 caller_name,
                 variants_by_sample[var_sample.name],
                 mutations_by_sample[var_sample.name],
                 vcf2txt_res_fpath)
                 for var_sample in var_samples)
        info('Done postprocessing all filtered VCFs.')

    except OSError:
        err(traceback.format_exc())
        warn('Running sequencially instead in ' + str(threads_num) + ' threads')
        try:
            Parallel(n_jobs=1) \
                (delayed(postprocess_vcf) \
                    (None,
                     cnf.work_dir,
                     var_sample,
                     caller_name,
                     variants_by_sample[var_sample.name],
                     mutations_by_sample[var_sample.name],
                     vcf2txt_res_fpath)
                     for var_sample in var_samples)
            info('Done postprocessing all filtered VCFs.')

        except OSError:
            err(traceback.format_exc())
            err('Cannot postprocess VCF - skipping')
            err()

    info('Filtered VCFs are written.')


def make_vcf2txt_cmdl_params(cnf, vcf_fpath_by_sample):
    c = cnf.variant_filtering
    min_freq = c.act_min_freq

    cmdline = \
        '-r 1.0 -R 1.0 -P {c.filt_p_mean} -Q {c.filt_q_mean} -D {c.filt_depth} -V {c.min_vd} ' \
        '-f {min_freq} -p {c.min_p_mean} -q {c.min_q_mean} ' \
        '-M {c.min_mq} -o {c.signal_noise} -L'.format(**locals())

    if c.bias:
        cmdline += ' -b '

    dbsnp_multi_mafs = cnf.genome.dbsnp_multi_mafs
    if dbsnp_multi_mafs and verify_file(dbsnp_multi_mafs):
        cmdline += ' -A ' + dbsnp_multi_mafs
    else:
        cmdline += ' -A ""'

    if c.amplicon_based:
        cmdline += ' -a '

    # corr_vcf_fpath_by_sample = dict()
    # for sn, vcf_fpath in vcf_fpath_by_sample.items():
    #     ungz = vcf_fpath
    #     if vcf_fpath.endswith('.gz'):
    #         ungz = splitext(vcf_fpath)[0]
    #         call(cnf, 'gunzip ' + vcf_fpath, output_fpath=ungz)
    #     corr_vcf_fpath_by_sample[sn] = ungz

    cmdline += ' ' + ' '.join(vcf_fpath_by_sample.values())
    return cmdline


def run_vcf2txt(cnf, vcf_fpath_by_sample, vcf2txt_out_fpath):
    info()
    info('Running VarDict vcf2txt...')

    vcf2txt = get_script_cmdline(cnf, 'perl', 'vcf2txt', is_critical=True)

    cmdline = vcf2txt + ' ' + make_vcf2txt_cmdl_params(cnf, vcf_fpath_by_sample)

    check_filtering_results(vcf2txt_out_fpath)

    res = run_vcf2txt_with_retries(cnf, cmdline, vcf2txt_out_fpath)
    return res


def join_vcf2txt_results(cnf, vcf_fpath_by_sample, vcf2txt_out_fpath):
    info('WGS; running vcftxt separately for each sample to save memory.')
    vcf2txt_outputs_by_vcf_fpath = OrderedDict()
    for vcf_fpath in vcf_fpath_by_sample.values():
        sample_output_fpath = add_suffix(vcf2txt_out_fpath, splitext(basename(vcf_fpath))[0])
        vcf2txt_outputs_by_vcf_fpath[vcf_fpath] = sample_output_fpath
        info()

    info('Joining vcf2txt ouputs... (' + str(len(vcf2txt_outputs_by_vcf_fpath)) +
         ' out of ' + str(len(vcf_fpath_by_sample)) + ' successful), ' +
         'writing to ' + vcf2txt_out_fpath)
    with file_transaction(cnf.work_dir, vcf2txt_out_fpath) as tx:
        with open(tx, 'w') as out:
            for i, (vcf_fpath, sample_output_fpath) in enumerate(vcf2txt_outputs_by_vcf_fpath.items()):
                info('   Reading ' + sample_output_fpath)
                with open(sample_output_fpath) as inp:
                    for j, l in enumerate(inp):
                        if j == 0 and i != 0:
                            continue
                        out.write(l)
    if verify_file(vcf2txt_out_fpath):
        info('Saved ' + vcf2txt_out_fpath)
        return vcf2txt_out_fpath
    else:
        return None


def run_vcf2txt_with_retries(cnf, cmdline, output_fpath):
    res = None
    tries = 0
    MAX_TRIES = 1
    WAIT_MINUTES = int(random() * 60) + 30
    err_fpath = join(cnf.work_dir, 'varfilter_' + splitext(basename(output_fpath))[0] + '.err')
    while True:
        stderr_dump = []
        output_didnt_exist = not verify_file(output_fpath, silent=True)
        res = call(cnf, cmdline, output_fpath, stderr_dump=stderr_dump, exit_on_error=False)
        if res is not None:
            return res
        else:
            tries += 1
            msg = 'vcf2txt.pl crashed:\n' + cmdline + ' > ' + output_fpath + '\n' + \
                  (''.join(['\t' + l for l in stderr_dump]) if stderr_dump else '')
            if tries < MAX_TRIES:
                msg += '\n\nRerunning in ' + str(WAIT_MINUTES) + ' minutes (tries ' + str(tries) + '/' + str(MAX_TRIES) + ')'

            send_email(cnf, msg_other=msg,
                       subj='vcf2txt.pl crashed [' + str(cnf.project_name) + ']',
                       only_me=True)
            err(msg)
            if tries == MAX_TRIES:
                break
            sleep(WAIT_MINUTES * 60)
            info()
            if output_didnt_exist and verify_file(output_fpath, silent=True):
                info('Output was created while sleeping: ' + output_fpath)
                break
    return res


# def count_cohort_freqs(cnf, samples, cohort_freqs_fpath, max_ratio):
#     info('Calculating frequences of varaints in the cohort')
#     if cnf.reuse_intermediate and verify_file(cohort_freqs_fpath, silent=True):
#         info(cohort_freqs_fpath + ' exists, reusing')
#         return cohort_freqs_fpath
#
#     freq_in_sample_by_vark = defaultdict(int)
#     for varks in Parallel(n_jobs=len(samples))(delayed(get_counts)(s.anno_vcf_fpath) for s in samples):
#         for vark in varks:
#             freq_in_sample_by_vark[vark] += 1
#
#     info('Counted ' + str(len(freq_in_sample_by_vark)) + ' variants in ' + str(len(samples)) + ' samples, '
#        'writing to ' + cohort_freqs_fpath)
#
#     with file_transaction(cnf.work_dir, cohort_freqs_fpath) as tx:
#         with open(tx, 'w') as out:
#             lines_written = 0
#             for vark, count in freq_in_sample_by_vark.items():
#                 freq = float(count) / len(samples)
#                 if freq > max_ratio:
#                     chrom, pos, ref, alt = vark
#                     out.write('chr' + chrom + ':' + str(pos) + ':' + ref + ':' + alt +
#                               '\t' + str(freq) + '\n')
#                     lines_written += 1
#     freq_in_sample_by_vark = None
#     info('Done, written ' + str(lines_written) + ' varks with freq > '
#          + str(max_ratio) + ' to ' + cohort_freqs_fpath)
#     return verify_file(cohort_freqs_fpath)
#
#
# def get_counts(vcf_fpath):
#     varks = []
#
#     with open_gzipsafe(vcf_fpath) as f:
#         for l in f:
#             if l.startswith('#'):
#                 continue
#             fs = l.split()
#             if len(fs) < 7:
#                 continue
#             pass_ = fs[6]
#             if pass_ == 'PASS':
#                 vark = (fs[0].replace('chr', ''), int(fs[1]), fs[3], fs[4])
#                 varks.append(vark)
#     info('Counted ' + str(len(varks)) + ' variants in ' + vcf_fpath)
#     return varks


def combine_results(cnf, samples, vcf2txt_fpaths, variants_fpath, pass_variants_fpath=None, reject_variants_fpath=None):
    info('Combining vcf2txt variants')
    not_existing_snames = []
    if cnf.reuse_intermediate and isfile(variants_fpath) and verify_file(variants_fpath):
        info('Combined filtered results ' + variants_fpath + ' exist, reusing.')
    else:
        for sample_i, (sample, vcf2txt_fpath) in enumerate(zip(samples, vcf2txt_fpaths)):
            if not verify_file(vcf2txt_fpath, description='variants file'):
                not_existing_snames.append(sample.name)
        if not_existing_snames:
            critical('For some samples do not exist, variants file was not found: ' + ', '.join(not_existing_snames))
        with file_transaction(cnf.work_dir, variants_fpath) as tx:
            with open(tx, 'w') as out:
                for sample_i, (sample, vcf2txt_fpath) in enumerate(zip(samples, vcf2txt_fpaths)):
                    with open(vcf2txt_fpath) as f:
                        for line_i, l in enumerate(f):
                            if line_i == 0 and sample_i == 0:
                                out.write(l)
                            if line_i > 0:
                                out.write(l)
        verify_file(variants_fpath, is_critical=True, description='combined mutation calls')
        info('Saved vcf2txt variants to ' + variants_fpath)

    info()
    info('Combining PASSed mutations')
    pass_variants_fpath = pass_variants_fpath or add_suffix(variants_fpath, source.mut_pass_suffix)
    reject_variants_fpath = reject_variants_fpath or add_suffix(variants_fpath, source.mut_reject_suffix)
    not_existing_pass_snames = []
    if cnf.reuse_intermediate and isfile(pass_variants_fpath) and verify_file(pass_variants_fpath)\
            and isfile(reject_variants_fpath) and verify_file(reject_variants_fpath):
        info('Combined PASSed filtered results ' + pass_variants_fpath + ' exist, reusing.')
    else:
        for sample_i, (sample, vcf2txt_fpath) in enumerate(zip(samples, vcf2txt_fpaths)):
            if not verify_file(add_suffix(vcf2txt_fpath, source.mut_pass_suffix), description='PASS variants file'):
                not_existing_pass_snames.append(sample.name)
        if not_existing_pass_snames:
            critical('For some samples do not exist, PASS variants file was not found: ' + ', '.join(not_existing_pass_snames))
        info('*' * 70)
        if cnf.variant_filtering.max_ratio < 1.0:
            info('Max ratio set to ' + str(cnf.variant_filtering.max_ratio))
        else:
            info('Max ratio set to ' + str(cnf.variant_filtering.max_ratio) + ', i.e. no filter')

        info('Calculating frequences of variants in the cohort')
        info('*' * 70)
        freq_in_cohort_by_vark, count_in_cohort_by_vark = count_mutations_freq(cnf, samples, vcf2txt_fpaths)
        reject_freq_in_cohort_by_vark, reject_count_in_cohort_by_vark = count_mutations_freq(
                cnf, samples, vcf2txt_fpaths, suffix=source.mut_reject_suffix)
        info()

        if cnf.variant_filtering.max_ratio < 1.0:
            info('Saving passing threshold if cohort freq < ' + str(cnf.variant_filtering.max_ratio) +
                 ' to ' + pass_variants_fpath)

        artefacts_samples, artefacts_data, variants_count, written_lines_count = write_combined_results(
            cnf, pass_variants_fpath, samples, vcf2txt_fpaths, freq_in_cohort_by_vark, count_in_cohort_by_vark,
            suffix=source.mut_pass_suffix, do_cohort_filtering=True)

        _, _, _, reject_written_lines_count = write_combined_results(cnf, reject_variants_fpath, samples, vcf2txt_fpaths,
            reject_freq_in_cohort_by_vark, reject_count_in_cohort_by_vark,
            suffix=source.mut_reject_suffix, do_cohort_filtering=False)

        if len(artefacts_samples.keys()) > 0:
            reason = 'cohort freq > ' + str(cnf.variant_filtering.max_ratio)
            with open(reject_variants_fpath) as f:
                line = f.readline().split()
                reason_col = line.index('Reason') if 'Reason' in line else None
            with open(reject_variants_fpath, 'a') as f:
                for vark, samples in artefacts_samples.items():
                    fs = artefacts_data[vark]
                    if reason_col:
                        fs[reason_col] = reason
                    else:
                        fs.append(reason)
                    f.write('\t'.join(fs) + '\n')

            info('Skipped artefacts with cohort freq > ' + str(cnf.variant_filtering.max_ratio) +
                 ' and sample count > ' + str(cnf.variant_filtering.max_sample_cnt) + ': ' + str(len(artefacts_samples.keys())))
            info('Added artefacts into ' + reject_variants_fpath)

        info('All variants not under filtering: ' + str(variants_count['not_filtered']))
        if len(artefacts_samples.keys()) > 0:
            info('Variants not under filtering with freq > ' + str(cnf.variant_filtering.max_ratio) + ': ' + str(variants_count['good_freq']))

        verify_file(pass_variants_fpath, 'PASS variants file', is_critical=True)
        info('Written ' + str(written_lines_count) + ' records to ' + pass_variants_fpath)
        info('Written ' + str(reject_written_lines_count + len(artefacts_samples.keys())) + ' rejected records to ' + reject_variants_fpath)

        variants_fpath = verify_file(variants_fpath, is_critical=True)
        pass_variants_fpath = verify_file(pass_variants_fpath, is_critical=True)

        if not_existing_snames or not_existing_pass_snames:
            return None, None

    return variants_fpath, pass_variants_fpath


def count_mutations_freq(cnf, samples, vcf2txt_fpaths, suffix=source.mut_pass_suffix):
    count_in_cohort_by_vark = defaultdict(int)
    total_varks = 0
    total_duplicated_count = 0
    total_records_count = 0
    for sample_i, (sample, vcf2txt_fpath) in enumerate(zip(samples, vcf2txt_fpaths)):
        met_in_this_sample = set()
        processed_fpath = add_suffix(vcf2txt_fpath, suffix)
        if not isfile(processed_fpath):
            critical(processed_fpath + ' does not exist; please, rerun VarFilter.')
        with open(processed_fpath) as f:
            for line_i, l in enumerate(f):
                if line_i > 0:
                    fs = l.replace('\n', '').split()
                    if not fs:
                        continue
                    chrom, pos, db_id, ref, alt = fs[1:6]
                    vark = ':'.join([chrom, pos, ref, alt])
                    if vark in met_in_this_sample:
                        if suffix == source.mut_pass_suffix:
                            total_duplicated_count += 1
                    else:
                        count_in_cohort_by_vark[vark] += 1
                        if suffix == source.mut_pass_suffix:
                            met_in_this_sample.add(vark)
                            total_varks += 1
                    total_records_count += 1

    if suffix == source.mut_pass_suffix:
        info('Counted ' + str(len(count_in_cohort_by_vark)) + ' different variants ' +
             'in ' + str(len(samples)) + ' samples with total ' + str(total_varks) + ' records')
        info('Duplicated varks for this sample: ' + str(total_duplicated_count) + ' out of total ' +
             str(total_records_count) + ' records. Duplicated were not counted into cohort frequencies.')

    freq_in_cohort_by_vark = dict()
    max_freq = 0
    for vark, count in count_in_cohort_by_vark.items():
        f = float(count) / len(samples)
        freq_in_cohort_by_vark[vark] = f
        if f > max_freq:
            max_freq = f

    if suffix == source.mut_pass_suffix:
        info('Maximum frequency in cohort is ' + str(max_freq))
    return freq_in_cohort_by_vark, count_in_cohort_by_vark


def write_combined_results(cnf, variants_fpath, samples, vcf2txt_fpaths, freq_in_cohort_by_vark, count_in_cohort_by_vark,
                           suffix=source.mut_pass_suffix, do_cohort_filtering=True):
    artefacts_samples = OrderedDefaultDict(list)
    artefacts_data = OrderedDict()

    variants_count = defaultdict(int)
    written_lines_count = 0
    status_col, reason_col, n_samples_col, n_var_col, pcnt_sample_col, ave_af_col, incidentalome_col \
        = None, None, None, None, None, None, None

    with file_transaction(cnf.work_dir, variants_fpath) as tx:
        with open(tx, 'w') as out:
            for sample_i, (sample, vcf2txt_fpath) in enumerate(zip(samples, vcf2txt_fpaths)):
                mut_fpath = add_suffix(vcf2txt_fpath, suffix)
                with file_transaction(cnf.work_dir, mut_fpath) as fixed_mut_fpath_tx:
                    with open(mut_fpath) as f, open(fixed_mut_fpath_tx, 'w') as fixed_f_out:
                        for line_i, l in enumerate(f):
                            fs = l.replace('\n', '').split('\t')
                            if line_i == 0 and sample_i == 0:
                                out.write(l)
                            if line_i == 0:
                                fixed_f_out.write(l)
                                if status_col is not None and status_col != fs.index('Significance'):
                                    critical('Different format in ' + mut_fpath + ': status_col=' +
                                             str(fs.index('Significance')) + ', but the first sample was ' + str(status_col) +
                                             ', please rerun VarFilter from the beginning')
                                status_col = fs.index('Significance')
                                reason_col = status_col + 1
                                n_samples_col = fs.index('N_samples')
                                n_var_col = fs.index('N_Var')
                                pcnt_sample_col = fs.index('Pcnt_sample')
                                ave_af_col = fs.index('Ave_AF')
                                if 'Incidentalome' in fs:
                                    incidentalome_col = fs.index('Incidentalome')
                            if line_i > 0:
                                fs = l.replace('\n', '').split('\t')
                                chrom, pos, db_id, ref, alt = fs[1:6]
                                vark = ':'.join([chrom, pos, ref, alt])
                                assert len(fs) > reason_col, 'len(fs)=' + str(len(fs)) + ' > reason_col=' + str(reason_col) + \
                                                             ' in ' + sample.name + ', ' + vcf2txt_fpath + ' for line\n' + l

                                freq = freq_in_cohort_by_vark[vark]
                                cnt = count_in_cohort_by_vark[vark]
                                fs[n_samples_col] = str(len(samples))
                                fs[n_var_col] = str(cnt)
                                fs[pcnt_sample_col] = str(freq)
                                fs[ave_af_col] = ''
                                l = '\t'.join(fs) + '\n'

                                if do_cohort_filtering:
                                    if fs[status_col] in ['known', 'likely']:
                                        variants_count['not_filtered'] += 1
                                    elif freq >= cnf.variant_filtering.max_ratio and cnt > cnf.variant_filtering.max_sample_cnt:
                                        artefacts_samples[vark].append(sample.name)
                                        # if incidentalome_col:
                                        #     fs.remove(fs[incidentalome_col])
                                        artefacts_data[vark] = fs
                                        continue
                                variants_count['good_freq'] += 1
                                fixed_f_out.write(l)
                                out.write(l)
                                written_lines_count += 1
    return artefacts_samples, artefacts_data, variants_count, written_lines_count
