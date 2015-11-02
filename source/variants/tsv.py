from collections import OrderedDict
import os
import shutil
from os.path import dirname, realpath, join, basename, isfile, pardir
from ext_modules import vcf_parser as vcf

from source.calling_process import call
from source.file_utils import intermediate_fname, verify_file, splitext_plus
from source.tools_from_cnf import get_java_tool_cmdline
from source.file_utils import file_transaction
from source.file_utils import file_exists
from source.logger import step_greetings, info, err, critical, warn
from source.variants.vcf_processing import read_sample_names_from_vcf, leave_main_sample, vcf_one_per_line, \
    get_sample_column_index


def make_tsv(cnf, vcf_fpath, samplename, main_sample_index=None):
    step_greetings('Exporting to TSV...')

    vcf_fpath = vcf_one_per_line(cnf, vcf_fpath)

    if main_sample_index is None:
        main_sample_index = get_sample_column_index(vcf_fpath, samplename)

    tsv_fpath = _extract_fields_new(cnf, vcf_fpath, samplename, main_sample_index)
    if not tsv_fpath:
        return tsv_fpath

    # manual_tsv_fields = cnf.annotation.get('tsv_fields')
    # if manual_tsv_fields:
    #     field_map = dict((rec.keys()[0], rec.values()[0]) for rec in manual_tsv_fields)
    #     if cnf.get('keep_intermediate'):
    #         info('Saved TSV file to ' + tsv_fpath)
    #     tsv_fpath = _rename_fields(cnf, tsv_fpath, field_map)
    #     if cnf.get('keep_intermediate'):
    #         info('Saved TSV file with nice names to ' + tsv_fpath)

    return tsv_fpath


def _extract_fields_new(cnf, vcf_fpath, samplename, main_sample_index=0):
    fname, _ = splitext_plus(basename(vcf_fpath))
    tsv_fpath = join(cnf.work_dir, fname + '.tsv')

    if cnf.get('reuse_intermediate'):
        if file_exists(tsv_fpath):
            info(tsv_fpath + ' exists, reusing')
            return tsv_fpath

    manual_tsv_fields = cnf.annotation['tsv_fields']
    if not manual_tsv_fields:
        return None

    all_fields = []
    basic_fields = []
    info_fields = []
    eff_fields = []
    gt_fields = []
    tumor_gt = 'GEN[' + str(main_sample_index) + '].'
    normal_gt = 'GEN[' + str(1 - main_sample_index) + '].'

    lines = []

    with open(vcf_fpath) as inp:
        reader = vcf.Reader(inp)

        info('TSV saver: Building field list')
        for f in [rec.keys()[0] for rec in manual_tsv_fields]:
            if f.startswith('GEN'):
                _f = f.split('.')[1]
                if len(reader.samples) > 0:
                    if _f in reader.formats:
                        gt_fields.append(_f)
                        all_fields.append(f.replace('GEN[*].', tumor_gt))
                        if len(reader.samples) > 1:
                            all_fields.append(f.replace('GEN[*].', normal_gt))
                else:
                    warn('Warning: ' + f + ' is not in VCF header FORMAT records')

            elif f in ['CHROM', 'POS', 'REF', 'ALT', 'ID', 'FILTER', 'QUAL']:
                all_fields.append(f)
                basic_fields.append(f)

            elif any(f.startswith(af) and af in reader.infos for af in ['EFF', 'ANN']):
                all_fields.append(f)
                eff_fields.append(f)

            else:
                if f in reader.infos:
                    info_fields.append(f)
                    all_fields.append(f)
                elif f == 'SAMPLE':
                    all_fields.append(f)
                else:
                    warn('Warning: ' + f + ' is not in VCF header INFO records')

        info('TSV saver: Iterating over records...')
        d = OrderedDict()
        for rec in reader:
            for f in basic_fields:
                d[f] = rec.__dict__[f]

            for f in info_fields:
                d[f] = rec.INFO[f] if f in rec.INFO else ''

            if 'SAMPLE' not in d:
                d['SAMPLE'] = samplename

            if eff_fields:
                eff = rec.INFO.get(eff_fields[0][:3])
                if not eff:
                    for f in eff_fields:
                        d[f] = ''
                else:
                    eff_fs = eff[0].split('|')
                    eff_d = dict()
                    for val, header in zip(eff_fs, ['ALLELE', 'EFFECT', 'IMPACT', 'GENE', 'GENEID', 'FEATURE', 'FEATUREID', 'BIOTYPE', 'RANK', 'HGVS_C', 'HGVS_P', 'CDNA_POSLEN', 'CDS_POSLEN', 'AA_POSLEN', 'DISTANCE', 'LOG']):
                        if 'POSLEN' in header:
                            eff_d[header.split('_')[0] + '_POS'] = val.split('/')[0] if val else ''
                            eff_d[header.split('_')[0] + '_LEN'] = val.split('/')[1] if val else ''
                        else:
                            eff_d[header] = val
                    #ANN=GA |3_prime_UTR_variant|MODIFIER|RPL22|RPL22|transcript|NM_000983.3|Coding|4/4|c.*173dupT|||||173|;
                    #Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'
                    for f in eff_fields:
                        d[f] = eff_d[f.split('.')[1]]

            if rec.FORMAT:
                for _f in gt_fields:
                    if _f in rec.FORMAT:
                        d[tumor_gt + _f] = rec.samples[main_sample_index][_f]
                        if len(rec.samples) > 1 - main_sample_index:
                            d[normal_gt + _f] = rec.samples[1 - main_sample_index][_f]
                        else:
                            d[normal_gt + _f] = ''
                    else:
                        d[tumor_gt + _f] = ''
                        d[normal_gt + _f] = ''

            fs = []
            for f in all_fields:
                v = d[f]
                fs.append(v if v != '.' else '')
            lines.append(fs)

    info('TSV saver: Adding GEN[*] fields both for sample and for matched normal...')
    field_map = dict()
    for rec in manual_tsv_fields:
        k = rec.keys()[0]
        v = rec.values()[0]
        if k.startswith('GEN[*].'):
            _f = k.split('.')[1]
            field_map[tumor_gt + _f] = v
            field_map[normal_gt + _f] = 'Matched_' + v
        else:
            field_map[k] = v

    info('TSV saver: Writing TSV to ' + tsv_fpath)
    with file_transaction(cnf.work_dir, tsv_fpath) as tx:
        with open(tx, 'w') as out:
            out.write('\t'.join(field_map[f] for f in all_fields) + '\n')
            for fs in lines:
                new_fs = []
                for f in fs:
                    if isinstance(f, list):
                        new_fs.append(','.join(map(str, f)))
                    elif f is None:
                        new_fs.append('')
                    else:
                        new_fs.append(str(f))
                out.write('\t'.join(new_fs) + '\n')

    info('TSV saver: saved ' + tsv_fpath)
    return tsv_fpath


# def _extract_fields(cnf, vcf_fpath, main_sample_index=0):
#     fname, _ = splitext_plus(basename(vcf_fpath))
#     tsv_fpath = join(cnf['work_dir'], fname + '.tsv')
#
#     if cnf.get('reuse_intermediate'):
#         if file_exists(tsv_fpath):
#             info(tsv_fpath + ' exists, reusing')
#             return tsv_fpath
#
#     all_format_fields = set()
#
#     # # Split FORMAT field and sample fields
#     # def proc_line(l, i):
#     #     if l.startswith('#'):
#     #         return l
#     #     vals = l.strip().split('\t')
#     #     if len(vals) <= 9:
#     #         return l
#     #     info_field = vals[7]
#     #     format_fields = vals[8].split(':')
#     #     sample_fields = vals[9].split(':')
#     #     for f, s in zip(format_fields, sample_fields):
#     #         if f not in ['DP', 'MQ']:
#     #             if f == 'GT':
#     #                 s = '"' + s + '"'
#     #             f = 'GEN[*].' + f
#     #             info_field += ';' + f + '=' + s
#     #             all_format_fields.add(f)
#     #     l = '\t'.join(vals[:7] + [info_field])
#     #     return l
#
#     # broken_format_column_vcf_fpath = iterate_file(cnf, vcf_fpath, proc_line, 'split_format_fields')
#
#     fields = None
#
#     _manual_tsv_fields = cnf.annotation['tsv_fields']
#     if _manual_tsv_fields:
#         fields = []
#
#         with open(vcf_fpath) as inp:
#             reader = vcf.Reader(inp)
#             infos =
#             formats = reader.formats
#
#             for f in [rec.keys()[0] for rec in _manual_tsv_fields]:
#                 if f.startswith('GEN'):
#                     _f = f.split('.')[1]
#                     if _f in formats:
#                         fields.append(f)
#                     else:
#                         warn('Warning: ' + f + ' is not in VCF header FORMAT records')
#
#                 elif f.startswith('EFF') or f in ['CHROM', 'POS', 'REF', 'ALT', 'ID', 'FILTER', 'QUAL']:
#                     fields.append(f)
#
#                 else:
#                     if f in infos:
#                         fields.append(f)
#                     else:
#                         warn('Warning: ' + f + ' is not in VCF header INFO records')
#
#     # else:
#         # first_line = next(l.strip()[1:].split() for l in open(vcf_fpath)
#         #   if l.strip().startswith('#CHROM'))
#         # basic_fields = [f for f in first_line[:9] if f != 'INFO'
#         #   and f != 'FORMAT' and f != sample_name]
#         # manual_annots = filter(lambda f: f and f != 'ID', all_fields)
#         # fields = (basic_fields + all_format_fields + manual_annots +
#         #    self.run_cnf.get('additional_tsv_fields', []))
#     else:
#         return None
#
#     if main_sample_index is None:  # No sample names or no main sample is ambiguous
#         fields = [f for f in fields if not f.startswith('GEN[*]')]
#     else:
#         fields = [f.replace('[*]', '[' + str(main_sample_index) + ']') for f in fields]
#
#     # info_fields = [f for f in fields if f not in ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'ID'm ]
#     # with open(broken_format_column_vcf_fpath) as vcf_f:
#     #     fields = filter_info_tsv_fileds(vcf_f, fields)
#
#     # fields = [f for f in fields if '[*]' not in f or 'EFF[*]' in f]
#
#     anno_line = ' '.join([f for f in fields if f != 'SAMPLE'])
#     snpsift = get_java_tool_cmdline(cnf, 'snpsift')
#     if not snpsift:
#         err('No SnpSIFT installed, skipping converting to TSV.')
#         return None
#     snpsift_cmdline = snpsift + ' extractFields ' + vcf_fpath + ' ' + anno_line
#
#     res = call(cnf, snpsift_cmdline, tsv_fpath, exit_on_error=False, print_stderr=True)
#     if res is None:
#         return None
#     info()
#
#     # REMOVE EMPTY, ADD SAMPLE COLUMN
#     with open(tsv_fpath) as tsv:
#         names, col_counts = None, None
#         for i, l in enumerate(tsv):
#             if i == 0:
#                 names = [v for v in l.split('\t') if v != '\n']
#                 col_counts = [0 for _ in names]
#             else:
#                 values = [v for v in l.split('\t') if v != '\n']
#
#                 for i, v in enumerate(values):
#                     if v:
#                         if len(col_counts) <= i:
#                             err('TSV file may be incorrectly generated (number of '
#                                 'column names (%d) is less than number of fields '
#                                 'in the record (%d)).' %
#                                 (len(names), len(values)))
#                         # while len(col_counts) <= i:
#                         #     col_counts.append(0)
#                         col_counts[i] += 1
#
#     with file_transaction(cnf.work_dir, tsv_fpath) as tx:
#         with open(tx, 'w') as out, open(tsv_fpath) as f:
#             for i, l in enumerate(f):
#                 values = [v for v in l.split('\t')]
#
#                 if i == 0:
#                     values[0] = values[0][1:]
#
#                 if fields[0] == 'SAMPLE':
#                     if i == 0:
#                         out.write('SAMPLE')
#                     else:
#                         out.write(cnf['name'])
#                     out.write('\t')
#
#                 # values = [v.replace('\n', '') for v in values]
#                 out.write('\t'.join([v.replace('\n', '') for j, v in enumerate(values)
#                                      if j < len(col_counts) and col_counts[j]]) + '\n')
#
#     # with file_transaction(cnf['tmp_dir'], tsv_fpath) as tx_tsv_fpath:
#     #     info(cmdline + ' < ' + (splitted_FORMAT_column_vcf_fpath
#     #                                         or vcf_fpath) + ' > ' + tx_tsv_fpath)
#     #     res = subprocess.call(cmdline,
#     #                           stdin=open(splitted_FORMAT_column_vcf_fpath or vcf_fpath),
#     #                           stdout=open(tx_tsv_fpath, 'w'), shell=True)
#
#     # if splitted_FORMAT_column_vcf_fpath:
#     #     os.remove(splitted_FORMAT_column_vcf_fpath)
#
#     # info('')
#     # if res != 0:
#     #     critical('Command returned status ' + str(res) +
#     #              ('. Log in ' + cnf['log'] if 'log' in cnf else '.'))
#
#     if not verify_file(tsv_fpath):
#         critical('Error: "snpsift extract" did not generate output file ' + tsv_fpath)
#     return tsv_fpath


def _rename_fields(cnf, inp_tsv_fpath, field_map):
    if cnf.get('keep_intermediate'):
        step_greetings('Renaming fields.')

    with open(inp_tsv_fpath) as f:
        first_line = f.readline()
    fields = first_line.split()
    new_fields = [field_map.get(f) or f for f in fields]
    new_first_line = '\t'.join(new_fields)

    if cnf.get('keep_intermediate'):
        out_tsv_fpath = intermediate_fname(cnf, inp_tsv_fpath, 'renamed')
    else:
        out_tsv_fpath = inp_tsv_fpath

    with file_transaction(cnf.work_dir, out_tsv_fpath) as tx_out_fpath:
        with open(tx_out_fpath, 'w') as out:
            out.write(new_first_line + '\n')
            with open(inp_tsv_fpath) as f:
                for i, l in enumerate(f):
                    if i >= 1:
                        out.write(l)

    if not cnf.get('keep_intermediate'):
        shutil.move(out_tsv_fpath, inp_tsv_fpath)
        return inp_tsv_fpath
    else:
        return out_tsv_fpath