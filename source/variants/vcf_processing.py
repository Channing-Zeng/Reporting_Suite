#!/usr/bin/env python
from genericpath import exists, isfile
import os
import sys
import shutil

from os.path import basename, join, expanduser, splitext, realpath, dirname, pardir
from collections import OrderedDict
##from memory_profiler import profile

from ext_modules import vcf_parser
from ext_modules.vcf_parser.model import _Record

from source.calling_process import call_subprocess, call
from source.change_checking import check_file_changed
from source.config import join_parent_conf
from source.file_utils import iterate_file, verify_file, intermediate_fname, convert_file, adjust_path, splitext_plus
from source.tools_from_cnf import get_java_tool_cmdline, get_tool_cmdline
from source.file_utils import file_transaction
from source.file_utils import open_gzipsafe, which, file_exists
from source.logger import step_greetings, info, critical, err


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

    def get_main_sample(self):
        if len(self._sample_indexes) == 0:
            return None
        try:
            sample_index = [sname.lower() for sname in self._sample_indexes]\
                .index(self.sample_name_from_file.lower())
        except ValueError:
            return self.samples[0]
        else:
            return self.samples[sample_index]

    def bias(self):
        if self._bias is None:
            bias_ = self.get_val('BIAS')
            if bias_ is not None:
                if not isinstance(bias_, basestring) or len(bias_) != 3 or bias_[1] not in [';', ':', ',']:
                    err('BIAS: ' + str(bias_) + ' for variant ' + self.get_variant())
                if bias_ == 0.0:
                    bias_ = None
                    err('Warning: BIAS is 0 ' + ' for variant ' + self.get_variant())
                else:
                    bias_.replace(';', ':').replace('.', ':').replace(',', ':')
            self._bias = bias_
        return self._bias

    def af(self):
        if self._af is None:
            self._af = self.get_val('AF')
        return self._af

    def get_val(self, key, default=None):
        if key in self.INFO:
            val = self.INFO[key]
            if not isinstance(val, basestring):
                try:
                    val = val[0]
                except:
                    pass
        else:
            main_sample = self.get_main_sample()
            if main_sample:
                sample_data = main_sample.data._asdict()
                if key in sample_data:
                    val = sample_data[key]

                    if not isinstance(val, basestring):
                        try:
                            val = ':'.join(val)
                        except:
                            pass
                else:
                    return default
            else:
                return default

        if val is not None:
            try:
                return float(val)
            except:
                try:
                    return int(val)
                except:
                    return val
        else:
            return default

    def cls(self):
        cls = 'Novel'
        if self.ID and 'COSM' in self.ID:
            cls = 'COSMIC'
        elif self.ID and self.ID.startswith('rs'):
            if self.check_clnsig:
                cls = 'ClnSNP'
            else:
                cls = 'dbSNP'
        return cls

    def is_rejected(self):
        if self.FILTER:
            assert '.' not in self.FILTER
        return self.FILTER and 'PASS' not in self.FILTER

    def check_clnsig(self):
        clnsig = self.INFO.get('CLNSIG')
        if not clnsig:
            return 0

        for c in self.INFO.get('CLNSIG'):
            if 4 <= c <= 6 or c == 255:
                return 1
            # 0 - Uncertain significance,
            # 1 - not provided,
            # 2 - Benign,
            # 3 - Likely benign,
            # 4 - Likely pathogenic,
            # 5 - Pathogenic,
            # 6 - drug response,
            # 7 - histocompatibility,
            # 255 - other 

        return -1

    def sample_field(self):
        return self.sample_name_from_file

    def get_variant(self):
        if self._variant is None:
            self._variant = ':'.join(map(str, [self.CHROM, self.POS, self.REF, self.ALT]))
        return self._variant


def iterate_vcf(cnf, input_fpath, proc_rec_fun, suffix=None,
                overwrite=False, reuse_intermediate=True, *args, **kwargs):
    #@profile
    def _convert_vcf(inp_f, out_f):
        max_bunch_size = 10000
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

    return convert_file(cnf, input_fpath, _convert_vcf,
                        suffix, overwrite, reuse_intermediate)


def vcf_is_empty(cnf, vcf_fpath):
    with open(vcf_fpath) as vcf:
        reader = vcf_parser.Reader(vcf)
        for rec in reader:
            return False
    return True


def remove_rejected(cnf, input_fpath):
    step_greetings('Removing rejeted records...')

    def _proc_rec(rec):
        if rec.FILTER is None or rec.FILTER == [] or rec.FILTER == ['PASS']:
            return rec
        else:
            return None

    out_fpath = iterate_vcf(cnf, input_fpath, _proc_rec, 'no_rej')

    if not verify_file(out_fpath):
        exit(1)

    return out_fpath

    #def __iter_file(l, i):
    #    if l.startswith('#'):
    #        return l
    #
    #    try:
    #        filt = l.split('\t')[6]
    #    except:
    #        if len(l.split('\t')) < 6 and len(l.split()) >= 6:
    #            critical('Error at line number ' + str(i) + ': fields separated by spaces rather than tabs?')
    #    else:
    #        if filt in ['PASS', '.']:
    #            return l
    #        else:
    #            return None

    #output_fpath = iterate_file(cnf, input_fpath, __iter_file)
    #info('Saved to ' + output_fpath)
    #return output_fpath


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
        if not verify_file(one_item_cnf['vcf'], 'Input file'):
            sys.exit(1)

        join_parent_conf(one_item_cnf, common_cnf)

        work_vcf = join(one_item_cnf['work_dir'], basename(one_item_cnf['vcf']))
        check_file_changed(one_item_cnf, one_item_cnf['vcf'], work_vcf)
        if not one_item_cnf.get('reuse_intermediate'):
            with open_gzipsafe(one_item_cnf['vcf']) as inp, open(work_vcf, 'w') as out:
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

                cnf['vcf'] = extract_sample(cnf, one_item_cnf['vcf'], cnf['name'])
                info()

        # SINGLE SAMPLE
        else:
            cnf = one_item_cnf

            if 'bam' in cnf:
                cnf['bam'] = adjust_path(cnf['bam'])
                if not verify_file(cnf['bam']):
                    sys.exit(1)

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
    db_path = cnf['genome'].get('snpeff')
    if not db_path:
        err('Please, provide a path to SnpEff data in '
            'the "genomes" section in the system config.')
    else:
        if isfile(dump_transcript_fpath):
            os.remove(dump_transcript_fpath)
        genome = cnf.genome.name
        cmdline = '{snpeff} dump -dataDir {db_path} -v -txt {genome}'.format(**locals())
        if call(cnf, cmdline, output_fpath=dump_transcript_fpath):
            cnf.transcripts_fpath = dump_transcript_fpath

    return cnf.transcripts_fpath


def convert_to_maf(cnf, vcf_fpath, tumor_sample_name, transcripts_fpath,
                   bam_fpath=None, normal_sample_name=None):
    step_greetings('Converting to MAF')

    vcf_fpath = vcf_one_per_line(cnf, vcf_fpath)

    #########################################################
    transcripts = '--transcripts ' + transcripts_fpath if transcripts_fpath else ''

    bam_fpath = bam_fpath or '.'
    normal_sample_name = normal_sample_name or '.'
    fname, _ = splitext_plus(basename(vcf_fpath))
    maf_fpath = join(cnf['work_dir'], fname + '.maf')

    perl = get_tool_cmdline(cnf, 'perl')
    vcf2maf = join(dirname(realpath(__file__)), '../../external/vcf2maf-1.1.0/vcf2maf.pl')
    cmdline = ('{perl} {vcf2maf} '
               '--input-snpeff {vcf_fpath} '
               '--bam-file {bam_fpath} '
               '--tumor-id {tumor_sample_name} '
               '--normal-id {normal_sample_name} '
               '--output-maf {maf_fpath} '
               '{transcripts} ').format(**locals())

    call(cnf, cmdline, output_fpath=maf_fpath, stdout_to_outputfile=False)

    if verify_file(maf_fpath, 'MAF'):
        info('MAF file saved to ' + maf_fpath)
    else:
        err('Converting to MAF didn\'t generate output file ' + maf_fpath)
        final_maf_fpath = None

    return maf_fpath


def fix_chromosome_names(cnf, vcf_fpath):
    step_greetings('Fixing chromosome names...')

    def _proc_rec(rec):
        if not rec.CHROM.startswith('chr'):
            rec.CHROM = 'chr' + rec.CHROM
        return rec

    out_fpath = iterate_vcf(cnf, vcf_fpath, _proc_rec, 'chr')

    if not verify_file(out_fpath):
        critical()

    return out_fpath


def vcf_one_per_line(cnf, vcf_fpath):
    info('Converting VCF to one-effect-per-line...')

    if not which('perl'):
        exit('Perl executable required, maybe you need to run "module load perl"?')

    src_fpath = join(dirname(realpath(__file__)))
    perl = get_tool_cmdline(cnf, 'perl')
    external_fpath = join(src_fpath, pardir, pardir, 'external')
    vcfoneperline_cmline = perl + ' ' + join(external_fpath, 'vcfOnePerLine.pl')
    oneperline_vcf_fpath = intermediate_fname(cnf, vcf_fpath, 'opl')

    call(cnf, vcfoneperline_cmline, oneperline_vcf_fpath, stdin_fpath=vcf_fpath, exit_on_error=False)
    info()

    if not verify_file(oneperline_vcf_fpath):
        critical('Error: vcf_one_per_line didn\'t generate output file.')
    return oneperline_vcf_fpath


def read_sample_names_from_vcf(vcf_fpath):
    basic_fields = next(
        (l.strip()[1:].split() for l in open(vcf_fpath)
        if l.strip().startswith('#CHROM')), None)
    if not basic_fields:
        critical('Error: no VCF header in ' + vcf_fpath)
    if len(basic_fields) < 9:
        return []
    return basic_fields[9:]


def get_main_sample_index(cnf, vcf_fpath, samplename):
    vcf_header_samples = read_sample_names_from_vcf(vcf_fpath)

    if len(vcf_header_samples) == 0:
        return None

    if len(vcf_header_samples) == 1:
        return 0

    name = next((name for name in vcf_header_samples if name.lower() in samplename.lower()), None)
    if name is None:
        err('No sample ' + samplename + ' in header with samples ' + ', '.join(vcf_header_samples) + ' for ' + vcf_fpath)
        name = next((name for name in vcf_header_samples if name.lower() != 'none'), None)
        if name is None:
            err('All sample names are None. Skipping.')
            return None

    try:
        return vcf_header_samples.index(name)
    except ValueError:
        return None


def leave_main_sample(cnf, vcf_fpath, samplename):
    index = get_main_sample_index(cnf, vcf_fpath, samplename)
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


def extract_sample(cnf, input_fpath, samplename):
    work_dir = cnf['work_dir']

    vcf_header_samples = read_sample_names_from_vcf(input_fpath)

    if len(vcf_header_samples) == 0:
        return input_fpath

    if len(vcf_header_samples) > 0 and samplename not in vcf_header_samples:
            critical('Error: not sample ' + samplename + ' in the ' + input_fpath +
                     ' header: available only ' + ', '.join(vcf_header_samples))

    if len(vcf_header_samples) == 1:
        return input_fpath

    info('-' * 70)
    info('Extracting sample ' + samplename)
    info('-' * 70)

    executable = get_java_tool_cmdline(cnf, 'gatk')
    ref_fpath = cnf['genome']['seq']

    corr_samplename = ''.join([c if c.isalnum() else '_' for c in samplename])

    output_fname = splitext_plus(input_fpath)[0] + '.' + corr_samplename + '.vcf'
    output_fpath = join(work_dir, output_fname)

    cmd = '{executable} -nt 30 -R {ref_fpath} -T SelectVariants ' \
          '--variant {input_fpath} -sn {samplename} ' \
          '-o {output_fpath}'.format(**locals())
    call_subprocess(cnf, cmd, None, output_fpath, stdout_to_outputfile=False,
         to_remove=[input_fpath + '.idx'])

    return output_fpath


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
    field_to_del = 'EFF'

    def proc_line(l, i):
        if l.startswith('##SnpEff'):
            return None

        elif field_to_del in l:
            if l.startswith('##INFO='):
                try:
                    if l.split('=', 1)[1].split(',', 1)[0].split('=')[1] == field_to_del:
                        return None
                except IndexError:
                    critical('Incorrect VCF at line: ' + l)

            elif not l.startswith('#'):
                fields = l.split('\t')
                info_line = fields[7]
                info_pairs = [attr.split('=') for attr in info_line.split(';')]
                info_pairs = filter(lambda pair: pair[0] != field_to_del, info_pairs)
                info_line = ';'.join('='.join(pair) if len(pair) == 2 and pair[0] != field_to_del
                                     else pair[0] for pair in info_pairs)
                fields = fields[:7] + [info_line] + fields[8:]
                return '\t'.join(fields)
        return l

    return iterate_file(cnf, input_fpath, proc_line, suffix='no_' + field_to_del.lower())


def igvtools_index(cnf, vcf_fpath):
    igvtools = get_tool_cmdline(cnf, 'igvtools')
    if not igvtools:
        err('Warning: no igvtools found, cannot index VCF.')
        return None
    if igvtools.endswith('.jar'):
        igvtools = get_java_tool_cmdline(cnf, 'igvtools')

    cmdline = '{igvtools} index {vcf_fpath}'.format(**locals())
    call(cnf, cmdline)
    if exists('igv.log'):
        try:
            os.remove('igv.log')
        except OSError:
            pass
    return vcf_fpath + '.idx'


def tabix_vcf(cnf, vcf_fpath):
    bgzip = get_tool_cmdline(cnf, 'bgzip')
    tabix = get_tool_cmdline(cnf, 'tabix')

    gzipped_fpath = join(vcf_fpath + '.gz')
    tbi_fpath = gzipped_fpath + '.tbi'

    if not bgzip:
        err('Cannot index VCF because bgzip is not found in PATH or '  + cnf.sys_cnf)
    if not tabix:
        err('Cannot index VCF because tabix is not found in PATH or '  + cnf.sys_cnf)
    if not bgzip and not tabix:
        return None, None

    info('BGzipping VCF')
    cmdline = '{bgzip} -c {vcf_fpath}'.format(**locals())
    call(cnf, cmdline, gzipped_fpath, overwrite=True)

    info('Tabixing VCF')
    cmdline = '{tabix} -f -p vcf {gzipped_fpath}'.format(**locals())
    call(cnf, cmdline, tbi_fpath, overwrite=True)

    return gzipped_fpath, tbi_fpath