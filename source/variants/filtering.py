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
from source.variants.vcf_processing import iterate_vcf, get_sample_column_index, bgzip_and_tabix, igvtools_index
from source.file_utils import safe_mkdir, add_suffix, verify_file, open_gzipsafe, \
    symlink_plus, file_transaction, num_lines
from source.logger import info
from source.webserver.exposing import convert_path_to_url


def combine_vcfs(cnf, vcf_fpath_by_sname, combined_vcf_fpath):
    gatk = get_java_tool_cmdline(cnf, 'gatk')
    cmdl = '{gatk} -T CombineVariants -R {cnf.genome.seq}'.format(**locals())
    for s_name, vcf_fpath in vcf_fpath_by_sname.items():
        cmdl += ' --variant:' + s_name + ' ' + vcf_fpath
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
        bgzip_and_tabix(cnf, combined_vcf_fpath)
        return combined_vcf_fpath
    else:
        warn('Could not join VCFs')
        return None


def index_vcf(sample_name, pass_vcf_fpath, filt_vcf_fpath, caller_name=None):
    global glob_cnf
    cnf = glob_cnf

    info()
    info(sample_name + ((', ' + caller_name) if caller_name else '') + ': indexing')

    for fpath in [pass_vcf_fpath, filt_vcf_fpath]:
        if not cnf.reuse_intermediate and not verify_file(fpath, silent=True):
            err(fpath + ' does not exist - cannot IGV index')
        else:
            if cnf.reuse_intermediate and verify_file(fpath + '.idx', silent=True):
                info('Reusing existing ' + fpath + '.idx')
            else:
                igvtools_index(cnf, fpath)

    if not cnf.reuse_intermediate and not verify_file(filt_vcf_fpath, silent=True):
        err(filt_vcf_fpath + ' does not exist - cannot gzip and tabix')
    else:
        if cnf.reuse_intermediate and verify_file(filt_vcf_fpath + '.gz', silent=True) \
                and verify_file(filt_vcf_fpath + '.gz.tbi', silent=True):
            info(filt_vcf_fpath + '.gz and .gz.tbi exist; reusing')
        else:
            bgzip_and_tabix(cnf, filt_vcf_fpath)


def filter_with_vcf2txt(cnf, var_samples, output_dirpath, vcf2txt_out_fpath,
        caller_name=None, sample_min_freq=None, threads_num=1):

    threads_num = min(len(var_samples), cnf.threads)
    info('Number of threads for filtering: ' + str(threads_num))

    safe_mkdir(output_dirpath)

    vcf_fpath_by_sample = {s.name: s.anno_vcf_fpath for s in var_samples}
    res = run_vcf2txt(cnf, vcf_fpath_by_sample, vcf2txt_out_fpath, sample_min_freq)
    if not res:
        err('vcf2txt run returned non-0')
        return None

    vardict2mut_perl = get_script_cmdline(cnf, 'perl', join('VarDict', 'vardict2mut.pl'), is_critical=True)
    vardict2mut_py = get_script_cmdline(cnf, 'python', join('scripts', 'post', 'vardict2mut.py'))

    res = run_vardict2mut(cnf, vcf2txt_out_fpath,
                          add_suffix(vcf2txt_out_fpath, 'pl.' + source.mut_pass_suffix),
                          sample_min_freq=sample_min_freq,
                          vardict2mut_executable=vardict2mut_perl)
    if not res:
        err('vardict2mut.pl run returned non-0')
    pl_mut_fpath = res

    if not vardict2mut_py:
        critical('vardict2mut_py not found')

    info('Running python version ' + vardict2mut_py)
    res = run_vardict2mut(cnf, vcf2txt_out_fpath,
        add_suffix(vcf2txt_out_fpath, source.mut_pass_suffix),
        sample_min_freq=sample_min_freq, vardict2mut_executable=vardict2mut_py)
    if not res:
        critical('vardict2mut.py run returned non-0')
    mut_fpath = py_mut_fpath = res
    if pl_mut_fpath:
        pl_mut_url = convert_path_to_url(pl_mut_fpath)
    py_mut_url = convert_path_to_url(py_mut_fpath)

    # Compare results, send email
    msg = 'VarDict2mut comparison\n'
    if pl_mut_fpath:
        msg += 'Perl:\n\t' + pl_mut_fpath + '\n'
        msg += '\t' + pl_mut_url + '\n'
    msg += 'Python:\n\t' + py_mut_fpath + '\n'
    msg += '\t' + py_mut_url + '\n\n'

    if pl_mut_fpath:
        py_line_num = num_lines(py_mut_fpath)
        pl_line_num = num_lines(pl_mut_fpath)
        if py_line_num == pl_line_num:
            msg += 'Line numbers is equal: ' + str(pl_line_num) + '\n'
        else:
            msg += 'Line numbers differ: perl (' + str(pl_line_num) + '), py (' + str(py_line_num) + ')\n'
        try:
            if call(cnf, get_system_path(cnf, 'diff') + ' -q ' + pl_mut_fpath +
                    ' ' + py_mut_fpath, exit_on_error=False, return_err_code=True) != 0:
                msg += 'Differ found.\n'
            else:
                msg += 'Files equal.\n'
        except:
            traceback.print_exc()
            msg += 'Diff failed with exception.\n'

    info(msg)
    # send_email(msg_other=msg, only_me=True)
    info()

    write_vcfs(cnf, var_samples,
               join(output_dirpath, source.varfilter_name),
               caller_name, vcf2txt_out_fpath, mut_fpath, threads_num)
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
        min_freq = defaults.default_min_freq

    cmdline = '{vardict2mut_executable} {vcf2txt_res_fpath} -f {min_freq} '
    if vardict2mut_executable.endswith('.pl'):
        cmdline += ' --report_reason '
        if c.min_hotspot_freq is not None and c.min_hotspot_freq != 'default':
            cmdline += ' -F ' + str(c.min_hotspot_freq)
        if c.max_ratio_vardict2mut is not None:
            cmdline += ' -R ' + str(c.max_ratio_vardict2mut)
        if cnf.genome.ruledir: cmdline += ' --ruledir {cnf.genome.ruledir} '
        if cnf.genome.filter_common_snp: cmdline += ' --filter_common_snp {cnf.genome.filter_common_snp} '
        if cnf.genome.filter_common_artifacts: cmdline += ' --filter_common_artifacts {cnf.genome.filter_common_artifacts} '
        if cnf.genome.actionable: cmdline += ' --actionable {cnf.genome.actionable} '
        if cnf.genome.compendia_ms7_hotspot: cmdline += ' --compendia_ms7_hotspot {cnf.genome.compendia_ms7_hotspot} '
        if cnf.snpeffect_export_polymorphic: cmdline += ' --snpeffect_export_polymorphic {cnf.snpeffect_export_polymorphic} '
        if cnf.actionable_hotspot: cmdline += ' --actionable_hotspot {cnf.actionable_hotspot} '
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
        var_sample.filt_tsv_fpath = splitext(var_sample.filt_vcf_fpath)[0] + '.tsv'

    if cnf.reuse_intermediate \
            and verify_file(gz) \
            and verify_file(var_sample.pass_filt_vcf_fpath)\
            and verify_file(var_sample.filt_tsv_fpath):
        info(var_sample.filt_vcf_fpath + '.gz' + ' and ' + var_sample.pass_filt_vcf_fpath + ' exist; reusing.')

    else:
        safe_mkdir(dirname(var_sample.filt_vcf_fpath))
        safe_mkdir(dirname(var_sample.pass_filt_vcf_fpath))

        with open_gzipsafe(var_sample.anno_vcf_fpath) as vcf_f, \
             file_transaction(work_dir, var_sample.filt_vcf_fpath) as filt_tx, \
             file_transaction(work_dir, var_sample.pass_filt_vcf_fpath) as pass_tx:
            with open(filt_tx, 'w') as filt_f, open(pass_tx, 'w') as pass_f:
                info(var_sample.name + ((', ' + caller_name) if caller_name else '') + ': opened ' +
                     var_sample.anno_vcf_fpath + ', writing to ' +
                     var_sample.filt_vcf_fpath + ' and ' + var_sample.pass_filt_vcf_fpath)

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
             var_sample.filt_vcf_fpath + ' and ' + var_sample.pass_filt_vcf_fpath)

        info()
        info(var_sample.name + ((', ' + caller_name) if caller_name else '') + ': writing filtered TSVs')
        # Converting to TSV - saving .anno.filt.tsv
        if 'tsv_fields' in cnf.annotation and cnf.tsv:
            tmp_tsv_fpath = make_tsv(cnf, var_sample.filt_vcf_fpath, var_sample.name)
            if not tmp_tsv_fpath:
                err('TSV convertion didn\'t work')
            else:
                if isfile(var_sample.filt_tsv_fpath):
                    os.remove(var_sample.filt_tsv_fpath)
                shutil.copy(tmp_tsv_fpath, var_sample.filt_tsv_fpath)

            info(var_sample.name + ((', ' + caller_name) if caller_name else '') +
                 ': saved filtered TSV to ' + var_sample.filt_tsv_fpath)

    info('Done postprocessing filtered VCF.')


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
    postprocess_vcf(cnf, cnf.work_dir, sample, caller_name, variants_dict, mutations, vcf2txt_res_fpath)


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


def run_vcf2txt(cnf, vcf_fpath_by_sample, vcf2txt_out_fpath, sample_min_freq=None):
    info()
    info('Running VarDict vcf2txt...')

    vcf2txt = get_script_cmdline(cnf, 'perl', join('VarDict', 'vcf2txt.pl'), is_critical=True)

    c = cnf.variant_filtering
    min_freq = cnf.min_freq or c.min_freq
    if min_freq == 'bcbio' or min_freq is None:
        min_freq = sample_min_freq
    if min_freq is None:
        min_freq = defaults.default_min_freq

    cmdline = '{vcf2txt} ' \
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

    if cnf.is_wgs:
        info('WGS; running vcftxt separately for each sample to save memory.')
        vcf2txt_outputs_by_vcf_fpath = OrderedDict()
        for vcf_fpath in vcf_fpath_by_sample.values():
            sample_output_fpath = add_suffix(vcf2txt_out_fpath, splitext(basename(vcf_fpath))[0])
            res = __run_vcf2txt(cnf, cmdline + ' ' + '<(gunzip -c ' + vcf_fpath + ')', sample_output_fpath)
            if res:
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

    else:
        cmdline += ' ' + ' '.join('<(gunzip -c ' + vcf_fpath + ')' for vcf_fpath in vcf_fpath_by_sample.values())
        res = __run_vcf2txt(cnf, cmdline, vcf2txt_out_fpath)
        return res


def __run_vcf2txt(cnf, cmdline, output_fpath):
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