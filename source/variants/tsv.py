import os
import shutil
from os.path import dirname, realpath, join, basename, isfile, pardir
from ext_modules import vcf_parser as vcf

from source.calling_process import call
from source.file_utils import intermediate_fname, verify_file, splitext_plus
from source.tools_from_cnf import get_java_tool_cmdline
from source.file_utils import file_transaction
from source.file_utils import file_exists
from source.logger import step_greetings, info, err, critical
from source.variants.vcf_processing import read_sample_names_from_vcf, leave_first_sample, vcf_one_per_line


def make_tsv(cnf, vcf_fpath):
    step_greetings('Exporting to TSV...')

    vcf_fpath = leave_first_sample(cnf, vcf_fpath)

    vcf_fpath = vcf_one_per_line(cnf, vcf_fpath)

    tsv_fpath = _extract_fields(cnf, vcf_fpath)
    if not tsv_fpath:
        return tsv_fpath

    manual_tsv_fields = cnf.get('tsv_fields')
    if manual_tsv_fields:
        field_map = dict((rec.keys()[0], rec.values()[0]) for rec in manual_tsv_fields)
        if cnf.get('keep_intermediate'):
            info('Saved TSV file to ' + tsv_fpath)
        tsv_fpath = _rename_fields(cnf, tsv_fpath, field_map)
        if cnf.get('keep_intermediate'):
            info('Saved TSV file with nice names to ' + tsv_fpath)

    return tsv_fpath


def filter_info_tsv_fileds(inp_f, fields):
    reader = vcf.Reader(inp_f)
    new_fields = fields[:]

    for i, rec in enumerate(reader):
        for field in fields:
            if field.split('[')[0] in rec.INFO:
                new_fields.append(field)

    return new_fields


def _extract_fields(cnf, vcf_fpath):
    fname, _ = splitext_plus(basename(vcf_fpath))
    tsv_fpath = join(cnf['work_dir'], fname + '.tsv')

    if cnf.get('reuse_intermediate'):
        if file_exists(tsv_fpath):
            info(tsv_fpath + ' exists, reusing')
            return tsv_fpath

    all_format_fields = set()

    # # Split FORMAT field and sample fields
    # def proc_line(l, i):
    #     if l.startswith('#'):
    #         return l
    #     vals = l.strip().split('\t')
    #     if len(vals) <= 9:
    #         return l
    #     info_field = vals[7]
    #     format_fields = vals[8].split(':')
    #     sample_fields = vals[9].split(':')
    #     for f, s in zip(format_fields, sample_fields):
    #         if f not in ['DP', 'MQ']:
    #             if f == 'GT':
    #                 s = '"' + s + '"'
    #             f = 'GEN[*].' + f
    #             info_field += ';' + f + '=' + s
    #             all_format_fields.add(f)
    #     l = '\t'.join(vals[:7] + [info_field])
    #     return l

    # broken_format_column_vcf_fpath = iterate_file(cnf, vcf_fpath, proc_line, 'split_format_fields')

    _manual_tsv_fields = cnf['tsv_fields']
    if _manual_tsv_fields:
        fields = [
            rec.keys()[0] for rec
            in _manual_tsv_fields]
    # else:
        # first_line = next(l.strip()[1:].split() for l in open(vcf_fpath)
        #   if l.strip().startswith('#CHROM'))
        # basic_fields = [f for f in first_line[:9] if f != 'INFO'
        #   and f != 'FORMAT' and f != sample_name]
        # manual_annots = filter(lambda f: f and f != 'ID', all_fields)
        # fields = (basic_fields + all_format_fields + manual_annots +
        #    self.run_cnf.get('additional_tsv_fields', []))
    else:
        return None

    sample_names = read_sample_names_from_vcf(vcf_fpath)
    if len(sample_names) == 0:
        fields = [f for f in fields if not f.startswith('GEN[*]')]

    # info_fields = [f for f in fields if f not in ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'ID'm ]
    # with open(broken_format_column_vcf_fpath) as vcf_f:
    #     fields = filter_info_tsv_fileds(vcf_f, fields)

    # fields = [f for f in fields if '[*]' not in f or 'EFF[*]' in f]


    anno_line = ' '.join([f for f in fields if f != 'SAMPLE'])
    snpsift = get_java_tool_cmdline(cnf, 'snpsift')
    snpsift_cmdline = snpsift + ' extractFields ' + vcf_fpath + ' ' + anno_line

    res = call(cnf, snpsift_cmdline, tsv_fpath, exit_on_error=False, print_stderr=True)
    if res is None:
        return None
    info()

    # REMOVE EMPTY, ADD SAMPLE COLUMN
    with open(tsv_fpath) as tsv:
        names, col_counts = None, None
        for i, l in enumerate(tsv):
            if i == 0:
                names = [v for v in l.split('\t') if v != '\n']
                col_counts = [0 for _ in names]
            else:
                values = [v for v in l.split('\t') if v != '\n']

                for i, v in enumerate(values):
                    if v:
                        if len(col_counts) <= i:
                            err('TSV file may be incorrectly generated (number of '
                                'column names (%d) is less than number of fields '
                                'in the record (%d)).' %
                                (len(names), len(values)))

                        col_counts[i] += 1

    with file_transaction(cnf, tsv_fpath) as tx:
        with open(tx, 'w') as out, open(tsv_fpath) as f:
            for i, l in enumerate(f):
                values = [v for v in l.split('\t')]

                if i == 0:
                    values[0] = values[0][1:]

                if fields[0] == 'SAMPLE':
                    if i == 0:
                        out.write('SAMPLE')
                    else:
                        out.write(cnf['name'])
                    out.write('\t')

                # values = [v.replace('\n', '') for v in values]
                out.write('\t'.join([v.replace('\n', '') for j, v in enumerate(values)
                                     if j < len(col_counts) and col_counts[j]]) + '\n')

    # with file_transaction(cnf['tmp_dir'], tsv_fpath) as tx_tsv_fpath:
    #     info(cmdline + ' < ' + (splitted_FORMAT_column_vcf_fpath
    #                                         or vcf_fpath) + ' > ' + tx_tsv_fpath)
    #     res = subprocess.call(cmdline,
    #                           stdin=open(splitted_FORMAT_column_vcf_fpath or vcf_fpath),
    #                           stdout=open(tx_tsv_fpath, 'w'), shell=True)

    # if splitted_FORMAT_column_vcf_fpath:
    #     os.remove(splitted_FORMAT_column_vcf_fpath)

    # info('')
    # if res != 0:
    #     critical('Command returned status ' + str(res) +
    #              ('. Log in ' + cnf['log'] if 'log' in cnf else '.'))

    if not verify_file(tsv_fpath):
        critical('Error: "snpsift extract" did not generate output file ' + tsv_fpath)
    return tsv_fpath


def _rename_fields(cnf, inp_tsv_fpath, field_map):
    if cnf.get('keep_intermediate'):
        step_greetings('Renaming fields.')

    with open(inp_tsv_fpath) as f:
        first_line = f.readline()
    fields = first_line.split()
    new_fields = [field_map.get(f, f) for f in fields]
    new_first_line = '\t'.join(new_fields)

    if cnf.get('keep_intermediate'):
        out_tsv_fpath = intermediate_fname(cnf, inp_tsv_fpath, 'renamed')
    else:
        out_tsv_fpath = inp_tsv_fpath

    with file_transaction(cnf, out_tsv_fpath) as tx_out_fpath:
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