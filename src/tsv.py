from genericpath import isfile
import subprocess
import os
from os.path import dirname, realpath, join, basename
import shutil

from src.transaction import file_transaction
from src.utils import which, splitext_plus, file_exists
from src.my_utils import critical, iterate_file, step_greetings, get_java_tool_cmdline, \
    intermediate_fname, info


def make_tsv(cnf, vcf_fpath):
    tsv_fpath = _extract_fields(cnf, vcf_fpath, cnf['work_dir'], cnf['name'])

    manual_tsv_fields = cnf.get('tsv_fields')
    if manual_tsv_fields:
        field_map = dict((rec.keys()[0], rec.values()[0]) for rec in manual_tsv_fields)
        if cnf.get('keep_intermediate'):
            info(cnf['log'], 'Saved TSV file to ' + tsv_fpath)
        tsv_fpath = _rename_fields(cnf, tsv_fpath, field_map, cnf['work_dir'])
        if cnf.get('keep_intermediate'):
            info(cnf['log'], 'Saved TSV file with nice names to ' + tsv_fpath)

    # Copying final TSV
    final_tsv_fname = splitext_plus(basename(cnf['vcf']))[0] + '.anno.tsv'
    final_tsv_fpath = join(cnf['output_dir'], final_tsv_fname)
    if isfile(final_tsv_fpath):
        os.remove(final_tsv_fpath)
    shutil.copyfile(tsv_fpath, final_tsv_fpath)

    return final_tsv_fpath


def _extract_fields(cnf, vcf_fpath, work_dir, sample_name=None):
    step_greetings(cnf, 'Extracting fields')

    name, _ = splitext_plus(vcf_fpath)
    tsv_fpath = name + '.tsv'

    if cnf.get('reuse_intermediate'):
        if file_exists(tsv_fpath):
            info(cnf['log'], tsv_fpath + ' exists, reusing')
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
    split_format_fields_vcf = iterate_file(cnf, vcf_fpath, proc_line, work_dir,
        'split_format_fields', keep_original_if_not_keep_intermediate=True)

    manual_tsv_fields = cnf.get('tsv_fields')
    if manual_tsv_fields:
        fields = [rec.keys()[0] for rec in manual_tsv_fields]
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

    anno_line = ' '.join(fields)
    snpsift_cmline = get_java_tool_cmdline(cnf, 'snpsift')

    if not which('perl'):
        exit('Perl executable required, maybe you need to run "module load perl"?')
    src_fpath = join(dirname(realpath(__file__)))
    vcfoneperline_cmline = 'perl ' + join(src_fpath, 'vcfOnePerLine.pl')
    cmdline = vcfoneperline_cmline + ' | ' + snpsift_cmline + ' extractFields - ' + anno_line

    with file_transaction(tsv_fpath) as tx_tsv_fpath:
        info(cnf['log'], cmdline + ' < ' + (split_format_fields_vcf
                                            or vcf_fpath) + ' > ' + tx_tsv_fpath)
        res = subprocess.call(cmdline,
                              stdin=open(split_format_fields_vcf or vcf_fpath),
                              stdout=open(tx_tsv_fpath, 'w'), shell=True)

    if split_format_fields_vcf:
        os.remove(split_format_fields_vcf)

    info(cnf['log'], '')
    if res != 0:
        critical(cnf['log'], 'Command returned status ' + str(res) +
                 ('. Log in ' + cnf['log'] if 'log' in cnf else '.'))
    return tsv_fpath


def _rename_fields(cnf, inp_tsv_fpath, field_map, work_dir):
    if cnf.get('keep_intermediate'):
        step_greetings(cnf, 'Renaming fields.')

    with open(inp_tsv_fpath) as f:
        first_line = f.readline()[1:]
    fields = first_line.split()
    new_fields = [field_map.get(f, f) for f in fields]
    new_first_line = '\t'.join(new_fields)

    if cnf.get('keep_intermediate'):
        out_tsv_fpath = intermediate_fname(work_dir, inp_tsv_fpath, 'renamed')
    else:
        out_tsv_fpath = inp_tsv_fpath

    with file_transaction(out_tsv_fpath) as tx_out_fpath:
        with open(tx_out_fpath, 'w') as out:
            out.write(new_first_line + '\n')
            with open(inp_tsv_fpath) as f:
                for i, l in enumerate(f):
                    if i >= 1:
                        out.write(l)

    if not cnf.get('keep_intermediate'):
        os.rename(out_tsv_fpath, inp_tsv_fpath)
        return inp_tsv_fpath
    else:
        return out_tsv_fpath