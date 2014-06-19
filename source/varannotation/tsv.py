import os
import shutil
from os.path import dirname, realpath, join, basename, isfile
from source.logger import step_greetings

from source.transaction import file_transaction
from source.utils_from_bcbio import which, splitext_plus, file_exists
from source.utils import iterate_file, get_java_tool_cmdline, \
    intermediate_fname, info, call_subprocess


def make_tsv(cnf, vcf_fpath):
    tsv_fpath = _extract_fields(cnf, vcf_fpath, cnf['work_dir'], cnf['name'])

    manual_tsv_fields = cnf.get('tsv_fields')
    if manual_tsv_fields:
        field_map = dict((rec.keys()[0], rec.values()[0]) for rec in manual_tsv_fields)
        if cnf.get('keep_intermediate'):
            info('Saved TSV file to ' + tsv_fpath)
        tsv_fpath = _rename_fields(cnf, tsv_fpath, field_map)
        if cnf.get('keep_intermediate'):
            info('Saved TSV file with nice names to ' + tsv_fpath)

    # Copying final TSV
    final_tsv_fname = splitext_plus(basename(cnf['vcf']))[0] + '.anno.tsv'
    final_tsv_fpath = join(cnf['output_dir'], final_tsv_fname)
    if isfile(final_tsv_fpath):
        os.remove(final_tsv_fpath)
    shutil.copyfile(tsv_fpath, final_tsv_fpath)

    return final_tsv_fpath


def _extract_fields(cnf, vcf_fpath, work_dir, sample_name=None):
    step_greetings('Extracting fields')

    name, _ = splitext_plus(basename(vcf_fpath))
    tsv_fpath = join(work_dir, name + '.tsv')

    if cnf.get('reuse_intermediate'):
        if file_exists(tsv_fpath):
            info(tsv_fpath + ' exists, reusing')
            return tsv_fpath

    all_format_fields = set()

    # Split FORMAT field and sample fields
    def proc_line(l):
        if l.startswith('#'):
            return l
        vals = l.strip().split('\t')
        if len(vals) <= 9:
            return l
        info_field = vals[7]
        format_fields = vals[8].split(':')
        sample_fields = vals[9].split(':')
        for f, s in zip(format_fields, sample_fields):
            if f not in ['DP', 'MQ']:
                if f == 'GT':
                    s = '"' + s + '"'
                f = 'gt_' + f
                info_field += ';' + f + '=' + s
                all_format_fields.add(f)
        l = '\t'.join(vals[:7] + [info_field])
        return l
    splitted_FORMAT_column_vcf_fpath = iterate_file(cnf, vcf_fpath, proc_line,
        'split_format_fields', keep_original_if_not_keep_intermediate=True)

    manual_tsv_fields = cnf['tsv_fields']
    print str(manual_tsv_fields)
    if manual_tsv_fields:
        fields_line = [
            rec.keys()[0] for rec
            in manual_tsv_fields
            if rec.keys()[0] != 'SAMPLE']
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

    anno_line = ' '.join(fields_line)
    snpsift_cmline = get_java_tool_cmdline(cnf, 'snpsift')

    if not which('perl'):
        exit('Perl executable required, maybe you need to run "module load perl"?')
    src_fpath = join(dirname(realpath(__file__)))
    vcfoneperline_cmline = 'perl ' + join(src_fpath, 'vcfOnePerLine.pl')

    cmdline = vcfoneperline_cmline + ' | ' + snpsift_cmline + \
              ' extractFields - ' + anno_line

    call_subprocess(cnf, cmdline, None, tsv_fpath,
         stdin_fpath=(splitted_FORMAT_column_vcf_fpath or vcf_fpath),
         to_remove=[splitted_FORMAT_column_vcf_fpath])

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
                        col_counts[i] += 1

    with file_transaction(cnf.tmp_dir, tsv_fpath) as tx:
        with open(tx, 'w') as out, open(tsv_fpath) as f:
            for i, l in enumerate(f):
                values = [v for v in l.split('\t')]

                if i == 0:
                    values[0] = values[0][1:]

                if manual_tsv_fields[0].keys()[0] == 'SAMPLE':
                    if i == 0:
                        out.write('SAMPLE')
                    else:
                        out.write(cnf['name'])
                    out.write('\t')

                out.write('\t'.join([v for j, v in enumerate(values) if col_counts[j] and '\n' not in v]) + '\n')

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

    with file_transaction(cnf['tmp_dir'], out_tsv_fpath) as tx_out_fpath:
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