from collections import OrderedDict, defaultdict
import os
from os.path import basename, join, isfile, islink, splitext, isdir, dirname
from random import random
from time import sleep
import traceback
import shutil

from joblib import Parallel, delayed
import source
from source.calling_process import call, call_check_output
from source.config import defaults
from source.profiling import fn_timer
from source.tools_from_cnf import get_script_cmdline, get_system_path, get_java_tool_cmdline
from source.logger import err, warn, send_email, critical
from source.utils import is_us
from source.variants.tsv import make_tsv
from source.variants.vcf_processing import iterate_vcf, get_sample_column_index, bgzip_and_tabix, igvtools_index, \
    verify_vcf
from source.file_utils import safe_mkdir, add_suffix, verify_file, open_gzipsafe, \
    symlink_plus, file_transaction, num_lines
from source.logger import info
from source.webserver.exposing import convert_path_to_url


def combine_vcfs(cnf, vcf_fpath_by_sname, combined_vcf_fpath):
    gatk = get_java_tool_cmdline(cnf, 'gatk')
    cmdl = '{gatk} -T CombineVariants -R {cnf.genome.seq}'.format(**locals())
    for s_name, vcf_fpath in vcf_fpath_by_sname.items():
        cmdl += ' --variant:' + s_name + ' ' + vcf_fpath
    if cnf.reuse_intermediate and isfile(combined_vcf_fpath + '.gz') and verify_vcf(combined_vcf_fpath + '.gz'):
        info(combined_vcf_fpath + '.gz exists, reusing')
        return combined_vcf_fpath

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


def index_vcf(cnf, sample_name, pass_vcf_fpath, filt_vcf_fpath, caller_name=None):
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
        caller_name=None, sample_min_freq=None, threads_num=1):

    threads_num = min(len(var_samples), cnf.threads)
    info('Number of threads for filtering: ' + str(threads_num))

    safe_mkdir(output_dirpath)

    vcf_fpath_by_sample = {s.name: s.anno_vcf_fpath for s in var_samples}
    res = run_vcf2txt(cnf, vcf_fpath_by_sample, vcf2txt_out_fpath, sample_min_freq)
    if not res:
        err('vcf2txt run returned non-0')
        return None

    vardict2mut_py = get_script_cmdline(cnf, 'python', join('scripts', 'post', 'vardict2mut.py'))
    if not vardict2mut_py:
        critical('vardict2mut_py not found')

    info('Running vardict2mut')
    res = run_vardict2mut(cnf, vcf2txt_out_fpath,
        add_suffix(vcf2txt_out_fpath, source.mut_pass_suffix),
        sample_min_freq=sample_min_freq, vardict2mut_executable=vardict2mut_py)
    if not res:
        critical('vardict2mut.py run returned non-0')
    mut_fpath = res
    mut_fpath = convert_path_to_url(mut_fpath)
    info()

    info('Done filtering with vcf2txt/vardict2mut, saved to ' + str(mut_fpath))
    return mut_fpath


@fn_timer
def run_vardict2mut(cnf, vcf2txt_res_fpath, vardict2mut_res_fpath=None,
                    sample_min_freq=None, vardict2mut_executable=None):
    cmdline = None
    if vardict2mut_res_fpath is None:
        vardict2mut_res_fpath = add_suffix(vcf2txt_res_fpath, source.mut_pass_suffix)

    if not vardict2mut_executable:
        vardict2mut_executable = get_script_cmdline(cnf, 'python', join('scripts', 'post', 'vardict2mut.py'))

    c = cnf.variant_filtering
    min_freq = cnf.min_freq or c.min_freq
    if min_freq == 'bcbio' or min_freq is None:
        min_freq = sample_min_freq
    if min_freq is None:
        min_freq = defaults['default_min_freq']

    cmdline = '{vardict2mut_executable} {vcf2txt_res_fpath} -f {min_freq} '
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
        cmdline += ' --sys-cnf ' + cnf.sys_cnf
        cmdline += ' --run-cnf ' + cnf.run_cnf
        if cnf.project_name:
            cmdline += ' --project-name ' + cnf.project_name
        cmdline += (' --reuse ' if cnf.reuse_intermediate else '')
        cmdline += ' --genome ' + cnf.genome.name
        cmdline += ' -o ' + vardict2mut_res_fpath
        if cnf.cohort_freqs_fpath:
            cmdline += ' --cohort-freqs ' + cnf.cohort_freqs_fpath
        cmdline = cmdline.format(**locals())
        res = call(cnf, cmdline, output_fpath=vardict2mut_res_fpath,
                   stdout_to_outputfile=False, exit_on_error=False)

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
            and verify_file(var_sample.filt_vcf_fpath) \
            and verify_file(var_sample.pass_filt_vcf_fpath) \
            and verify_file(var_sample.filt_tsv_fpath):
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


def make_vcf2txt_cmdl_params(cnf, vcf_fpath_by_sample, sample_min_freq=None):
    c = cnf.variant_filtering
    min_freq = cnf.min_freq or c.min_freq
    if min_freq == 'bcbio' or min_freq is None:
        min_freq = sample_min_freq
    if min_freq is None:
        min_freq = defaults['default_min_freq']

    cmdline = \
        '-f {min_freq} -n {c.sample_cnt} -F {c.ave_freq} -p {c.min_p_mean} -q {c.min_q_mean} ' \
        '-r {c.fraction} -R {c.max_ratio} -P {c.filt_p_mean} -Q {c.filt_q_mean} -D {c.filt_depth} ' \
        '-M {c.min_mq} -V {c.min_vd} -G {c.maf} -o {c.signal_noise} -L'.format(**locals())

    if c.bias:
        cmdline += ' -b '

    if c.count_undetermined is False:
        cmdline += ' -u '

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


def run_vcf2txt(cnf, vcf_fpath_by_sample, vcf2txt_out_fpath, sample_min_freq=None):
    info()
    info('Running VarDict vcf2txt...')

    vcf2txt = get_script_cmdline(cnf, 'perl', 'vcf2txt', is_critical=True)

    cmdline = vcf2txt + ' ' + \
        make_vcf2txt_cmdl_params(cnf, vcf_fpath_by_sample, sample_min_freq=sample_min_freq)

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

            send_email(msg_other=msg,
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


def count_cohort_freqs(cnf, samples, cohort_freqs_fpath, max_ratio):
    info('Calculating frequences of varaints in the cohort')
    if cnf.reuse_intermediate and verify_file(cohort_freqs_fpath, silent=True):
        info(cohort_freqs_fpath + ' exists, reusing')
        return cohort_freqs_fpath

    freq_in_sample_by_vark = defaultdict(int)
    for varks in Parallel(n_jobs=len(samples))(delayed(get_counts)(s.anno_vcf_fpath) for s in samples):
        for vark in varks:
            freq_in_sample_by_vark[vark] += 1

    info('Counted ' + str(len(freq_in_sample_by_vark)) + ' variants in ' + str(len(samples)) + ' samples, '
       'writing to ' + cohort_freqs_fpath)

    with file_transaction(cnf.work_dir, cohort_freqs_fpath) as tx:
        with open(tx, 'w') as out:
            lines_written = 0
            for vark, count in freq_in_sample_by_vark.items():
                freq = float(count) / len(samples)
                if freq > max_ratio:
                    chrom, pos, ref, alt = vark
                    out.write('chr' + chrom + ':' + str(pos) + ':' + ref + ':' + alt +
                              '\t' + str(freq) + '\n')
                    lines_written += 1
    freq_in_sample_by_vark = None
    info('Done, written ' + str(lines_written) + ' varks with freq > '
         + str(max_ratio) + ' to ' + cohort_freqs_fpath)
    return verify_file(cohort_freqs_fpath)


def get_counts(vcf_fpath):
    varks = []

    with open_gzipsafe(vcf_fpath) as f:
        for l in f:
            if l.startswith('#'):
                continue
            fs = l.split()
            if len(fs) < 7:
                continue
            pass_ = fs[6]
            if pass_ == 'PASS':
                vark = (fs[0].replace('chr', ''), int(fs[1]), fs[3], fs[4])
                varks.append(vark)
    info('Counted ' + str(len(varks)) + ' variants in ' + vcf_fpath)
    return varks
