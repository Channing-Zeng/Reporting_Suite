#!/usr/bin/env python
import traceback
from genericpath import exists, isfile, getctime
import os
import sys
import shutil

from os.path import basename, join, expanduser, splitext, realpath, dirname, pardir, abspath
from collections import OrderedDict
##from memory_profiler import profile
import re

from source.variants import vcf_parser
from source.variants.vcf_parser.model import _Record

from source.calling_process import call_subprocess, call
from source.change_checking import check_file_changed
from source.config import join_parent_conf
from source.file_utils import iterate_file, verify_file, intermediate_fname, convert_file, adjust_path, splitext_plus, \
    add_suffix
from source.targetcov.bam_and_bed_utils import verify_bam
from source.tools_from_cnf import get_java_tool_cmdline, get_system_path, get_script_cmdline
from source.file_utils import file_transaction
from source.file_utils import open_gzipsafe, which, file_exists
from source.logger import step_greetings, info, critical, err, warn, debug
from source.variants.Effect import Effect


class Record(_Record):
    # noinspection PyMissingConstructor
    def __init__(self, _record, input_fpath, line_num):
        self.__dict__.update(_record.__dict__)
        self.line_num = line_num
        file_base_name = basename(input_fpath)
        self.sample_name_from_file = file_base_name.split('-')[0]
        self._bias = None
        self._variant = None
        self._af = None
        self._t_ref_count = None
        self._t_alt_count = None

    def is_rejected(self):
        if self.FILTER:
            assert '.' not in self.FILTER
        return self.FILTER and 'PASS' not in self.FILTER

    def get_variant(self):
        if self._variant is None:
            self._variant = ':'.join(map(str, [self.CHROM, self.POS, self.REF, self.ALT]))
        return self._variant


def verify_vcf(vcf_fpath, silent=False, is_critical=False):
    if not verify_file(vcf_fpath, silent=silent, is_critical=is_critical):
        return None
    debug('File ' + vcf_fpath + ' exists and not empty')
    vcf = open_gzipsafe(vcf_fpath)
    debug('File ' + vcf_fpath + ' opened')
    l = next(vcf, None)
    if l is None:
        (critical if is_critical else err)('Error: cannot read the VCF file ' + vcf_fpath)
        return None
    if not l.startswith('##fileformat=VCF'):
        (critical if is_critical else err)('Error: VCF must start with ##fileformat=VCF ' + vcf_fpath)
        return None

    try:
        reader = vcf_parser.Reader(vcf)
    except:
        err('Error: cannot open the VCF file ' + vcf_fpath)
        if is_critical: raise
    else:
        debug('File ' + vcf_fpath + ' opened as VCF')
        try:
            rec = next(reader)
        except IndexError:
            err('Error: cannot parse records in the VCF file ' + vcf_fpath)
            debug('IndexError parsing VCF file ' + vcf_fpath)
            if is_critical: raise
        except ValueError:
            err('Error: cannot parse records in the VCF file ' + vcf_fpath)
            debug('ValueError parsing VCF file ' + vcf_fpath)
            if is_critical: raise
        except StopIteration:
            debug('No records in the VCF file ' + vcf_fpath)
            if not silent:
                warn('VCF file ' + vcf_fpath + ' has no records.')
            return vcf_fpath
        except:
            err('Error: cannot parse records in the VCF file ' + vcf_fpath)
            debug('Other error parsing VCF file ' + vcf_fpath)
            if is_critical: raise
        else:
            debug('A record was read from the VCF file ' + vcf_fpath)
            return vcf_fpath
        # f = open_gzipsafe(output_fpath)
        # l = f.readline()
        # if 'Cannot allocate memory' in l:
        #     f.close()
        #     f = open_gzipsafe(output_fpath)
        #     contents = f.read()
        #     if not silent:
        #         if is_critical:
        #             critical('SnpSift failed with memory issue:\n' + contents)
        #         else:
        #             err('SnpSift failed with memory issue:\n' + contents)
        #             return None
        #     f.close()
        #     return None
        # return output_fpath
    finally:
        vcf.close()


def iterate_vcf(cnf, input_fpath, proc_rec_fun, suffix=None, check_result=True,
                overwrite=False, reuse_intermediate=True, *args, **kwargs):
    info('iterate_vcf: overwrite=' + str(overwrite))
    def _convert_vcf(inp_f, out_f):
        max_bunch_size = 100000
        written_records = 0
        bunch = []

        reader = vcf_parser.Reader(inp_f)
        writer = vcf_parser.Writer(out_f, reader)

        i = 0
        while True:
            rec = next(reader, None)
            if rec is None:
                break

            rec = proc_rec_fun(Record(rec, input_fpath, i), *args, **kwargs)
            if rec:
                bunch.append(rec)
                written_records += 1

            if len(bunch) >= max_bunch_size:
                writer.write_records(bunch)
                info('Written lines: ' + str(written_records))
                bunch = []
            i += 1

        writer.write_records(bunch)
        bunch = []
        info('Written lines: ' + str(written_records))

    info('before convert_file: overwrite=' + str(overwrite))
    res = convert_file(
        cnf, input_fpath, _convert_vcf, suffix,
        check_result=check_result, overwrite=overwrite, reuse_intermediate=reuse_intermediate)
    return verify_vcf(res, is_critical=check_result)


def vcf_is_empty(cnf, vcf_fpath):
    vcf = open_gzipsafe(vcf_fpath)
    reader = vcf_parser.Reader(vcf)
    result = True
    for rec in reader:
        result = False
    vcf.close()
    return result


def _get_qual_threshold(input_fpath):
    qual_threshold = None
    q_filter_regex = re.compile(r'##FILTER=<ID=q(\d+),Description="Mean Base Quality Below \d+">')
    with open_gzipsafe(input_fpath) as f:
        for l in f:
            if not l.startswith('##'):
                break
            m = q_filter_regex.match(l)
            if m:
                qual_threshold = int(m.group(1))
                break
    return qual_threshold


def remove_rejected(cnf, input_fpath, output_fpath=None):
    # if not input_fpath.endswith('.gz') or not file_exists(input_fpath + '.tbi'):
    #     input_fpath = bgzip_and_tabix(cnf, input_fpath)

    qual_threshold = _get_qual_threshold(input_fpath)
    info('VCF QUAL threshold is ' + str(qual_threshold))
    if qual_threshold > cnf.variant_filtering.min_q_mean:
        info('Requested QUAL threshold is ' + str(cnf.variant_filtering.min_q_mean) +
             ', which is higher than in VCF, so keeping records with FILTER=q' + str(qual_threshold))

    def fn(l, i):
        if l.startswith('#'):
            return l
        else:
            fs = l.split('\t')
            if fs[6] == 'q' + str(qual_threshold) and qual_threshold > cnf.variant_filtering.min_q_mean:
                fs[6] = 'PASS'
            if fs[6] == 'PASS':
                return l
            else:
                return None
    return iterate_file(cnf, input_fpath, fn, suffix='pass')

    # bcftools = get_system_path(cnf, 'bcftools')
    # output_fpath = output_fpath or add_suffix(input_fpath, 'pass')
    # if output_fpath.endswith('.gz'):
    #     output_fpath = splitext(output_fpath)[0]
    # cmdl = '{bcftools} filter -i \'FILTER="PASS"|FILTER="."\' {input_fpath}'.format(**locals())
    # return call(cnf, cmdl, output_fpath=output_fpath)


def read_samples_info_and_split(common_cnf, options, inputs):
    #TODO: _set_up_dirs(cnf) for each sample

    info('')
    info('Processing input details...')

    details = None
    for key in inputs:
        if options.get(key):
            common_cnf[key] = adjust_path(options[key])
            info('Using ' + common_cnf[key])
            details = [common_cnf]
    if not details:
        details = common_cnf.get('details')
    if not details:
        critical('Please, provide input ' + ', '.join(inputs) +
                 ' in command line or in run info yaml config.')

    all_samples = OrderedDict()

    for one_item_cnf in details:
        if 'vcf' not in one_item_cnf:
            critical('ERROR: A section in details does not contain field "var".')
        one_item_cnf['vcf'] = adjust_path(one_item_cnf['vcf'])
        verify_file(one_item_cnf['vcf'], 'Input file', is_critical=True)

        join_parent_conf(one_item_cnf, common_cnf)

        work_vcf = join(one_item_cnf['work_dir'], basename(one_item_cnf['vcf']))
        check_file_changed(one_item_cnf, one_item_cnf['vcf'], work_vcf)
        if not one_item_cnf.get('reuse_intermediate'):
            with open_gzipsafe(one_item_cnf['vcf']) as inp, open_gzipsafe(work_vcf, 'w') as out:
                out.write(inp.read())
        one_item_cnf['vcf'] = work_vcf

        vcf_header_samples = read_sample_names_from_vcf(one_item_cnf['vcf'])

        # MULTIPLE SAMPELS
        if ('samples' in one_item_cnf or one_item_cnf.get('split_samples')) and len(vcf_header_samples) == 0:
            sample_cnfs = _verify_sample_info(one_item_cnf, vcf_header_samples)

            for header_sample_name in vcf_header_samples:
                if header_sample_name not in sample_cnfs:
                    sample_cnfs[header_sample_name] = one_item_cnf.copy()

                if header_sample_name in all_samples:
                    critical('ERROR: duplicated sample name: ' + header_sample_name)

                cnf = all_samples[header_sample_name] = sample_cnfs[header_sample_name]
                cnf['name'] = header_sample_name
                if cnf.get('keep_intermediate'):
                    cnf['log'] = join(cnf['work_dir'], cnf['name'] + '.log')

                # cnf['vcf'] = extract_sample(cnf, one_item_cnf['vcf'], cnf['name'])
                info()

        # SINGLE SAMPLE
        else:
            cnf = one_item_cnf

            if 'bam' in cnf:
                cnf['bam'] = adjust_path(cnf['bam'])
                verify_bam(cnf['bam'], is_critical=True)

            cnf['name'] = splitext_plus(basename(cnf['vcf']))[0]

            if cnf.get('keep_intermediate'):
                cnf['log'] = join(cnf['work_dir'], cnf['name'] + '.log')

            cnf['vcf'] = work_vcf
            all_samples[cnf['name']] = cnf

    if not all_samples:
        info('No samples.')
    else:
        info('Using samples: ' + ', '.join(all_samples) + '.')

    return all_samples


def get_trasncripts_fpath(cnf):
    if cnf.transcripts_fpath:
        if verify_file(cnf.transcripts_fpath):
            return cnf.transcripts_fpath

        if isfile(cnf.transcripts_fpath):
            os.remove(cnf.transcripts_fpath)

    # custom_transcripts_fpath = cnf['snpeff'].get('only_transcripts')
    # if custom_transcripts_fpath:
    #     if verify_file(custom_transcripts_fpath, 'Transcripts for snpEff -onlyTr'):
    #         transcripts_fpath = custom_transcripts_fpath
    #
    # else:

    dump_transcript_fpath = join(cnf.work_dir, 'snpeff_transcripts.txt')
    if isfile(dump_transcript_fpath) and verify_file(dump_transcript_fpath):
        cnf.transcripts_fpath = dump_transcript_fpath
        return cnf.transcripts_fpath

    snpeff = get_java_tool_cmdline(cnf, 'snpeff')
    if not snpeff:
        critical('No snpeff or it is incorrect path in system config.')
    db_path = cnf['genome'].get('snpeff')
    if db_path:
        db_path_cmdline = ' -dataDir ' + db_path
    else:
        db_path_cmdline = ''
        # err('Please, provide a path to SnpEff data in '
        #     'the "genomes" section in the system config.')
    if isfile(dump_transcript_fpath):
        os.remove(dump_transcript_fpath)
    genome = cnf.genome.name
    cmdline = '{snpeff} dump {db_path_cmdline} -v -txt {genome}'.format(**locals())
    if call(cnf, cmdline, output_fpath=dump_transcript_fpath):
        cnf.transcripts_fpath = dump_transcript_fpath

    return cnf.transcripts_fpath


# def _prepare_fields_for_maf_converter(cnf, vcf_fpath, sample_name):
#     main_sample_index = get_sample_column_index(vcf_fpath, sample_name)
#
#     def proc_line(rec):
#         main_sample = rec.get_main_sample(main_sample_index)
#         if main_sample:
#             sample_data = main_sample.data._asdict()
#
#             for key in ['BIAS', 'QUAL', 'QMEAN', 'PMEAN', 'PSTD', 'QSTD',
#                         'SBF', 'ODDRATIO', 'MQ', 'SN']:
#                 if key in sample_data and key not in rec.INFO:
#                     rec.INFO[key] = sample_data[key]
#
#         if rec.ID:
#             rec.INFO['COSMIC_overlapping_mutations'] = ','.join(
#                 id for id in rec.ID.split(';') if 'COSM' in id) or None
#         return rec
#
#     return iterate_vcf(cnf, vcf_fpath, proc_line, 'maf')
#
#
# def convert_to_maf(cnf, vcf_fpath, tumor_sample_name, transcripts_fpath,
#                    bam_fpath=None, normal_sample_name=None):
#     step_greetings('Converting to MAF')
#
#     if transcripts_fpath:
#         transcripts_fpath_copy = join(cnf.work_dir, tumor_sample_name + '_' + basename(transcripts_fpath))
#         if isfile(transcripts_fpath_copy) and not file_exists(transcripts_fpath_copy):
#             os.remove(transcripts_fpath_copy)
#
#         if not file_exists(transcripts_fpath_copy):
#             info('Copying transcripts file ' + transcripts_fpath + ' to ' + transcripts_fpath_copy)
#             shutil.copyfile(transcripts_fpath, transcripts_fpath_copy)
#             transcripts_fpath = transcripts_fpath_copy
#
#     info('Preparing VCF for MAF conversion...')
#     vcf_fpath = _prepare_fields_for_maf_converter(cnf, vcf_fpath, tumor_sample_name)
#     info()
#
#     vcf_fpath = vcf_one_per_line(cnf, vcf_fpath)
#
#     #########################################################
#     transcripts = '--transcripts ' + transcripts_fpath if transcripts_fpath else ''
#
#     bam_fpath = bam_fpath or '.'
#     normal_sample_name = normal_sample_name or '.'
#     fname, _ = splitext_plus(basename(vcf_fpath))
#     maf_fpath = join(cnf['work_dir'], fname + '.maf')
#
#     perl = get_system_path(cnf, 'perl')
#     vcf2maf = join(dirname(realpath(__file__)), 'vcf2maf.pl')
#     cmdline = ('{perl} {vcf2maf} '
#                '--input-snpeff {vcf_fpath} '
#                '--bam-file {bam_fpath} '
#                '--tumor-id {tumor_sample_name} '
#                '--normal-id {normal_sample_name} '
#                '--output-maf {maf_fpath} '
#                '{transcripts} ').format(**locals())
#
#     res = call(cnf, cmdline, output_fpath=maf_fpath, stdout_to_outputfile=False, exit_on_error=False)
#     if not res:
#         return None
#
#     if verify_file(maf_fpath, 'MAF'):
#         warn('MAF file saved to ' + maf_fpath)
#     else:
#         err('Converting to MAF didn\'t generate output file ' + maf_fpath)
#         final_maf_fpath = None
#
#     return maf_fpath


def fix_chromosome_names(cnf, vcf_fpath):
    with open(vcf_fpath) as f:
        for l in f:
            if not l.startswith('#'):
                if l.startswith('chr'):
                    info('Chomosome names are hg19, no need to fix.')
                    return vcf_fpath

    step_greetings('Fixing chromosome names')

    def _proc_rec(rec):
        if not rec.CHROM.startswith('chr'):
            rec.CHROM = 'chr' + rec.CHROM
        return rec

    out_fpath = iterate_vcf(cnf, vcf_fpath, _proc_rec, 'chr')

    if not verify_file(out_fpath):
        err('Could not run fix_chromosome_names')

    return out_fpath


def vcf_one_per_line(cnf, vcf_fpath):
    info('Converting VCF to one-effect-per-line...')

    oneperline_vcf_fpath = intermediate_fname(cnf, vcf_fpath, 'opl')
    vcfoneperline_cmline = get_script_cmdline(cnf, 'perl', join('ext_tools', 'vcfOnePerLine.pl'))
    call(cnf, vcfoneperline_cmline, oneperline_vcf_fpath, stdin_fpath=vcf_fpath, exit_on_error=False)
    info()

    if not verify_file(oneperline_vcf_fpath):
        critical('Error: vcf_one_per_line didn\'t generate output file.')
    return oneperline_vcf_fpath


def read_sample_names_from_vcf(vcf_fpath):
    f = open_gzipsafe(vcf_fpath)
    basic_fields = next(
        (l.strip()[1:].split() for l in f
        if l.strip().startswith('#CHROM')), None)
    if not basic_fields:
        critical('Error: no VCF header in ' + vcf_fpath)
    if len(basic_fields) < 9:
        return []
    return basic_fields[9:]


def get_sample_column_index(vcf_fpath, samplename, suppress_warn=False):
    vcf_header_samples = read_sample_names_from_vcf(vcf_fpath)

    if len(vcf_header_samples) == 0:
        return None

    if len(vcf_header_samples) == 1:
        if vcf_header_samples[0].lower() == samplename.lower():
            return 0
        else:
            return None

    name = next((name for name in vcf_header_samples if name.lower() == samplename.lower()), None)
    if name is None:
        if not suppress_warn:
            warn('No sample ' + samplename + ' in header with samples ' + ', '.join(vcf_header_samples) + ' for ' + vcf_fpath)
        name = next((name for name in vcf_header_samples if name.lower() != 'none'), None)
        if name is None:
            err('All sample names are None.')
        return None
    else:
        return vcf_header_samples.index(name)


def leave_main_sample(cnf, vcf_fpath, samplename):
    index = get_sample_column_index(vcf_fpath, samplename)
    if index is None:
        return vcf_fpath

    # def _f1(rec):
    #     rec.samples = [sample_name]
    #     return rec
    #
    info('Keeping SAMPLE only for the first sample (' + samplename + ')')
    # vcf_fpath = iterate_vcf(cnf, vcf_fpath, _f1, suffix=sample_name)
    # out_fpath = extract_sample(cnf, vcf_fpath, sample_name)
    # info()

    def _f(line, i):
        if line and (line.startswith('#CHROM') or line[0] != '#'):
            ts = line.split('\t')
            return '\t'.join(ts[:9] + [ts[9 + index]])
        return line
    vcf_fpath = iterate_file(cnf, vcf_fpath, _f, suffix='1sm')

    if not verify_file(vcf_fpath):
        err('Error: leave_first_sample didnt generate output file.')
        return None

    return vcf_fpath


def _verify_sample_info(vcf_conf, vcf_header_samples):
    if 'samples' in vcf_conf:
        for header_sample_name, sample_conf in vcf_conf['samples'].items():
            join_parent_conf(sample_conf, vcf_conf)

            bam = sample_conf.get('bam')
            if bam and not verify_file(bam, 'Bam file'):
                exit()
            sample_conf['bam'] = adjust_path(bam)

    sample_cnfs = vcf_conf.get('samples') or OrderedDict()

    # compare input sample names to vcf header
    if sample_cnfs:
        for input_sample_name, sample_conf in sample_cnfs.items():
            if input_sample_name not in vcf_header_samples:
                critical('ERROR: sample ' + input_sample_name +
                         ' is not in VCF header ' + vcf_header_samples + '\n'
                         'Available samples: ' + ', '.join(vcf_header_samples))
    return sample_cnfs


def remove_prev_eff_annotation(cnf, input_fpath):
    fields_to_del = ['EFF', 'ANN']

    def proc_line(l, i):
        if l.startswith('##SnpEff'):
            return None

        elif any(f in l for f in fields_to_del):
            if l.startswith('##INFO='):
                try:
                    if l.split('=', 1)[1].split(',', 1)[0].split('=')[1] in fields_to_del:
                        return None
                except IndexError:
                    critical('Incorrect VCF at line: ' + l)

            elif not l.startswith('#'):
                fields = l.split('\t')
                info_line = fields[7]
                info_pairs = [attr.split('=') for attr in info_line.split(';')]
                info_pairs = filter(lambda pair: pair[0] not in fields_to_del, info_pairs)
                info_line = ';'.join('='.join(pair) if len(pair) == 2 and pair[0] not in fields_to_del
                                     else pair[0] for pair in info_pairs)
                fields = fields[:7] + [info_line] + fields[8:]
                return '\t'.join(fields)
        return l

    return iterate_file(cnf, input_fpath, proc_line, suffix='noEFF')


def igvtools_index(cnf, vcf_fpath):
    igvtools = get_system_path(cnf, 'igvtools')
    if not igvtools:
        err('Warning: no igvtools found, cannot index VCF.')
        return None
    if igvtools.endswith('.jar'):
        igvtools = get_java_tool_cmdline(cnf, 'igvtools')
        if igvtools is None:
            err('Warning: no jar igvtools found, cannot index VCF.')
            return None

    cmdline = '{igvtools} index {vcf_fpath}'.format(**locals())
    call(cnf, cmdline, exit_on_error=False)
    if exists('igv.log'):
        try:
            os.remove('igv.log')
        except OSError:
            pass
    return vcf_fpath + '.idx'


def bgzip_and_tabix(cnf, vcf_fpath, tabix_parameters='', **kwargs):
    gzipped_fpath = join(vcf_fpath + '.gz')
    tbi_fpath = gzipped_fpath + '.tbi'

    if cnf.reuse_intermediate and \
           file_exists(gzipped_fpath) and \
           file_exists(tbi_fpath) and getctime(tbi_fpath) >= getctime(gzipped_fpath):
        info('Actual compressed VCF and index exist, reusing')
        return gzipped_fpath

    info('Compressing and tabixing VCF file, writing ' + gzipped_fpath + '(.tbi)')
    bgzip = get_system_path(cnf, 'bgzip')
    tabix = get_system_path(cnf, 'tabix')
    if not bgzip:
        err('Cannot index VCF because bgzip is not found in PATH or ' + cnf.sys_cnf)
    if not tabix:
        err('Cannot index VCF because tabix is not found in PATH or ' + cnf.sys_cnf)
    if not bgzip and not tabix:
        return vcf_fpath

    retrying = False
    while True:
        if isfile(tbi_fpath): os.remove(tbi_fpath)
        if isfile(vcf_fpath):
            if isfile(gzipped_fpath):
                 os.remove(gzipped_fpath)
            info('BGzipping VCF')
            cmdline = '{bgzip} {vcf_fpath}'.format(**locals())
            call(cnf, cmdline, None, **kwargs)
        else:
            if not verify_file(gzipped_fpath):
                err('Neither uncompressed ' + vcf_fpath + ' nor ' + gzipped_fpath + ' exist')
                return None

        info('Tabixing VCF')
        cmdline = '{tabix} {tabix_parameters} {gzipped_fpath}'.format(**locals())

        exit_on_error = False
        if retrying:
            exit_on_error = True
        kwargs['exit_on_error'] = exit_on_error
        call(cnf, cmdline, **kwargs)
        if isfile(gzipped_fpath + '.tbi'):
            break
        if retrying:
            critical('Cannot tabix ' + vcf_fpath)
        if not isfile(vcf_fpath):
            call(cnf, 'gunzip ' + gzipped_fpath, None)
        retrying = True

    return gzipped_fpath


def vcf_merge(cnf, vcf_fpaths, combined_vcf_fpath):
    vcf_merge_cmdline = get_system_path(cnf, join('ext_tools', 'vcftools', 'scripts', 'vcf-merge'))
    if vcf_merge_cmdline is None:
        critical('No vcf_merge in path')

    cmdline = vcf_merge_cmdline + ' ' + ' '.join(vcf_fpaths)
    perl_module_dirpath = abspath(join(dirname(__file__), pardir, pardir, 'ext_modules', 'perl_modules'))
    os.environ['PERL5LIB'] = perl_module_dirpath

    res = call(cnf, cmdline, combined_vcf_fpath, exit_on_error=False)
    if not res:
        return None
