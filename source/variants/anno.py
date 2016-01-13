from collections import OrderedDict
import shutil
import os
from os.path import splitext, basename, join, isfile, relpath
import socket
import sys
import re

import source
from source.calling_process import call_subprocess, call
from source.file_utils import iterate_file, intermediate_fname, verify_file, splitext_plus, add_suffix, file_transaction, \
    safe_mkdir
from source.logger import step_greetings, critical, info, err, warn
from source.tools_from_cnf import get_system_path, get_java_tool_cmdline, get_snpeff_type
from source.file_utils import file_exists, code_base_path
from source.variants import qc
from source.variants.vcf_processing import iterate_vcf, remove_prev_eff_annotation, bgzip_and_tabix, igvtools_index


def intersect_vcf(cnf, input_fpath, db_fpath, key):
    vcf_fpath = input_fpath

    db_fpath = verify_file(db_fpath)
    if not db_fpath:
        return None

    info('Intersecting with ' + db_fpath + ', writing key ' + str(key))

    info('Preparing db...')
    def _add_info_flag(l, i):
        if l.startswith('#'):
            return l
        fs = l.split('\t')
        info_col, ft_keys, ft_vals = fs[-3], fs[-2], fs[-1]
        ft_dict = dict(zip(ft_keys.split(':'), ft_vals.split(':')))
        for ann in ['DP', 'MQ']:
            val = ft_dict.get(ann, None)
            if val:
                # ft_keys[key.replace('.', '_') + '_' + ann] = val
                # del ft_keys[ann]
                info_col += ';' + key.replace('.', '_') + '_' + ann + '=' + val
        # ft_items = ft_dict.items()
        # ft_keys = [k for k, v in ft_items]
        # ft_vals = [v for k, v in ft_items]
        # return '\t'.join(fs[:-2]) + '\t' + ':'.join(ft_keys) + '\t' + ':'.join(ft_vals)
        return '\t'.join(fs[:-3]) + '\t' + info_col + '\t' + '\t'.join(fs[-2:])
        # rec.FORMAT[key.replace('.', '_') + '_DP'] = rec.genotype(key.split('.')[0])['DP']
        # rec.INFO[key.replace('.', '_') + '_MQ'] = rec.genotype(key.split('.')[0])['MQ']
        # return rec
    # db_fpath = iterate_vcf(cnf, db_fpath, _add_info_flag, suffix='INFO_FLAGS')
    db_fpath = iterate_file(cnf, db_fpath, _add_info_flag, suffix='INFO_FLAGS')

    info('Adding header meta info...')
    def _add_header(l, i):
        if l.startswith('#CHROM'):
            ext_l = ''
            for ann in ['DP', 'MQ']:
                ext_l += '##INFO=<ID=' + key.replace('.', '_') + '_' + ann + ',Number=1,Type=Integer,Description="description">\n'
            return ext_l + l
        return l
    db_fpath = iterate_file(cnf, db_fpath, _add_header, suffix='INFO_HEADER')

    # out_fpath = add_suffix(db_fpath, 'HEADERS')
    # if cnf.reuse_intermediate and verify_file(out_fpath, silent=True):
    #     info(out_fpath + ' exists, reusing')
    # else:
    #     reader = vcf_parser.Reader(open(db_fpath))
    #     for k in 'DP', 'MQ':
    #         k = k + '_' + key.replace('.', '_')
    #         reader.infos[k] = _Info(id=k, num=1, type='Integer', desc=k + ' ' + key)
    #
    #     with file_transaction(cnf.work_dir, out_fpath) as tx:
    #         recs = []
    #         cnt = 0
    #         with open(tx, 'w') as f:
    #             writer = vcf_parser.Writer(f, reader)
    #             while True:
    #                 cnt += 1
    #                 rec = next(reader, None)
    #                 if rec is None:
    #                     break
    #                 recs.append(rec)
    #                 if cnt % 1000000 == 0:
    #                     info('Written ' + str(cnt) + ' lines')
    #                     writer.write_records(recs)
    #                     recs = []
    #             writer.write_records(recs)
    #     db_fpath = out_fpath
    db_fpath = bgzip_and_tabix(cnf, db_fpath)

    info('Annotating using this db...')
    vcf_conf = {
        'path': db_fpath,
        'annotations': [key.replace('.', '_') + '_DP', key.replace('.', '_') + '_MQ']}
    vcf_fpath = _snpsift_annotate(cnf, vcf_conf, key, vcf_fpath)

    info('Moving INFO to FORMAT...')
    def _move_info_to_format(l, i):
        if l.startswith('#'):
            return l
        fs = l.split('\t')
        info_dict = dict([kv.split('=') if '=' in kv else (kv, True) for kv in fs[7].split(';')])
        ft_keys = fs[8].split(':')
        all_ft_vals = [ft_vals.split(':') for ft_vals in fs[9:]]
        ft_dicts = [OrderedDict(zip(ft_keys, ft_vals)) for ft_vals in all_ft_vals]
        for ann in ['DP', 'MQ']:
            k = key.replace('.', '_') + '_' + ann
            for ft_dict in ft_dicts:
                ft_dict[k] = info_dict.get(k, '.')
        all_ft_vals = []
        for ft_dict in ft_dicts:
            ft_items = ft_dict.items()
            ft_keys = [k for k, v in ft_items]
            all_ft_vals.append([v for k, v in ft_items])
        l = '\t'.join(fs[:8]) + '\t' + ':'.join(ft_keys)
        for ft_vals in all_ft_vals:
            l += '\t' + ':'.join(ft_vals)
        return l
        # rec.FORMAT[key.replace('.', '_') + '_DP'] = rec.genotype(key.split('.')[0])['DP']
        # rec.INFO[key.replace('.', '_') + '_MQ'] = rec.genotype(key.split('.')[0])['MQ']
        # return rec
    # db_fpath = iterate_vcf(cnf, db_fpath, _add_info_flag, suffix='INFO_FLAGS')
    vcf_fpath = iterate_file(cnf, vcf_fpath, _move_info_to_format, suffix='FORMAT_FLAGS')

    info('Adding FORMAT header meta info...')
    def _add_format_header(l, i):
        if l.startswith('#CHROM'):
            ext_l = ''
            ext_l += '##FORMAT=<ID=' + key.replace('.', '_') + '_DP,Number=1,Type=Integer,Description="Number of high-quality bases">\n'
            ext_l += '##FORMAT=<ID=' + key.replace('.', '_') + '_MQ,Number=1,Type=Integer,Description="Average mapping quality">\n'
            return ext_l + l
        return l
    vcf_fpath = iterate_file(cnf, vcf_fpath, _add_format_header, suffix='FORMAT_HEADER')

    info()
    if vcf_fpath:
        info('Renaming ' + vcf_fpath + ' -> ' + input_fpath)
        os.rename(vcf_fpath, input_fpath)
    else:
        warn('Intersection with ' + key + ' didn\'t work')
    return input_fpath


def run_annotators(cnf, vcf_fpath, bam_fpath):
    info('run_annotators')
    annotated = False
    original_vcf = cnf.vcf

    db_section_by_name = OrderedDict((dbname, cnf.annotation[dbname])
        for dbname in ['dbsnp', 'clinvar', 'cosmic', 'oncomine']
        if dbname in cnf.annotation and not cnf.annotation[dbname].get('skip-annotation'))

    if not cnf.no_check:
        to_delete_id_ref = []
        if 'dbsnp' in db_section_by_name.keys():
            info('Removing IDs from dbsnp as rs*')
            to_delete_id_ref.append('rs')
        if 'cosmic' in db_section_by_name.keys():
            info('Removing IDs from dbsnp as COS*')
            to_delete_id_ref.append('COS')

        def delete_ids(rec):  # deleting existing dbsnp and cosmic ID annotations
            if rec.ID:
                if isinstance(rec.ID, basestring):
                    if any(rec.ID.startswith(pref) for pref in to_delete_id_ref):
                        rec.ID = None
                else:
                    rec.ID = [id_ for id_ in rec.ID if not any(id_.startswith(pref) for pref in to_delete_id_ref)]

            if not rec.FILTER:
                rec.FILTER = 'PASS'

            return rec

        info('Removing previous rs* and COS* IDs')
        vcf_fpath = iterate_vcf(cnf, vcf_fpath, delete_ids, suffix='delID')

    for dbname, dbconf in db_section_by_name.items():
        res = _snpsift_annotate(cnf, dbconf, dbname, vcf_fpath)
        if res:
            vcf_fpath = res
            annotated = True

    if 'custom_vcfs' in cnf.annotation:
        for dbname, custom_conf in cnf.annotation['custom_vcfs'].items():
            res = _snpsift_annotate(cnf, custom_conf, dbname, vcf_fpath)
            if res:
                vcf_fpath = res
                annotated = True

    if 'dbnsfp' in cnf.annotation:
        res = _snpsift_db_nsfp(cnf, vcf_fpath)
        if res:
            vcf_fpath = res
            annotated = True

    if 'snpeff' in cnf.annotation:
        res, summary_fpath, genes_fpath = _snpeff(cnf, vcf_fpath)
        if res:
            vcf_fpath = res
            annotated = True

            final_summary_fpath = join(cnf.output_dir, basename(summary_fpath))
            final_genes_fpath = join(cnf.output_dir, basename(genes_fpath))
            if isfile(final_summary_fpath): os.remove(final_summary_fpath)
            if isfile(final_genes_fpath): os.remove(final_genes_fpath)
            if file_exists(summary_fpath): shutil.move(summary_fpath, final_summary_fpath)
            if file_exists(genes_fpath): shutil.move(genes_fpath, final_genes_fpath)

    if 'tracks' in cnf.annotation and cnf.annotation['tracks'] and cnf.annotation['tracks']:
        track_fapths = []
        for track_name in cnf.annotation['tracks']:
            if isfile(track_name) and verify_file(track_name):
                track_fapths.append(track_name)
            else:
                if 'tracks' in cnf['genome'] and cnf['genome']['tracks'] and track_name in cnf['genome']['tracks']:
                    track_fpath = cnf['genome']['tracks'][track_name]
                    if verify_file(track_fpath):
                        track_fapths.append(track_fpath)
        for track_fapth in track_fapths:
            res = _tracks(cnf, track_fapth, vcf_fpath)
            if res:
                annotated = True
                vcf_fpath = res

    step_greetings('Intersection with database VCFs...')
    if 'intersect_with' in cnf.annotation:
        for key, db_fpath in cnf.annotation['intersect_with'].items():
            res = intersect_vcf(cnf, input_fpath=vcf_fpath, db_fpath=db_fpath, key=key)
            if res:
                annotated = True
                vcf_fpath = res

    if 'mongo' in cnf.annotation:
        res = _mongo(cnf, vcf_fpath)
        if res:
            vcf_fpath = res
            annotated = True

    if not annotated:
        warn('Warning: No annotations were applied to ' + original_vcf)

    return annotated, vcf_fpath


def finialize_annotate_file(cnf, vcf_fpath, sample, callername):
    # vcf_fpath = leave_first_sample(cnf, vcf_fpath)

    if not cnf.no_check:
        vcf_fpath = _filter_malformed_fields(cnf, vcf_fpath)

    if not cnf.no_check:
        info()
        info('Adding SAMPLE=' + sample.name + ' annotation...')
        vcf_fpath = _add_annotation(cnf, vcf_fpath, 'SAMPLE', sample.name, number='1', type_='String', description='Sample name')

    final_vcf_fname = sample.name + (('-' + callername) if callername else '') + '.anno.vcf'
    final_vcf_fpath = join(cnf.output_dir, final_vcf_fname)

    info('Moving final VCF ' + vcf_fpath + ' to ' + final_vcf_fpath)
    if isfile(final_vcf_fpath):
        os.remove(final_vcf_fpath)
    if isfile(final_vcf_fpath + '.gz'):
        os.remove(final_vcf_fpath + '.gz')
    shutil.copy(vcf_fpath, final_vcf_fpath)

    if cnf.qc:
        report = qc.make_report(cnf, vcf_fpath, sample)
        qc_dirpath = join(cnf.output_dir, 'qc')
        safe_mkdir(qc_dirpath)
        qc.save_report(cnf, report, sample, callername, qc_dirpath, source.varqc_name)

    # Indexing
    info()
    info('Compressing and indexing with bgzip+tabix ' + final_vcf_fpath)
    final_vcf_fpath = bgzip_and_tabix(cnf, final_vcf_fpath)

    return [final_vcf_fpath]


def _mongo(cnf, input_fpath):
    step_greetings('Annotating from Mongo')

    if 'mongo' not in cnf.annotation:
        return None

    executable = get_java_tool_cmdline(cnf, join('ext_tools', 'mongo_loader', 'VCFStore.jar'))
    output_fpath = intermediate_fname(cnf, input_fpath, 'mongo')
    project_name = cnf.project_name

    cmdline = (
        '{executable} -module annotation -inputFile {input_fpath} ' ''
        '-outputFile {output_fpath} -project {project_name} ').format(**locals())
    if call_subprocess(cnf, cmdline, input_fpath, output_fpath,
                       stdout_to_outputfile=False, exit_on_error=False):
        return output_fpath
    else:
        return None


def _snpsift_annotate(cnf, vcf_conf, dbname, input_fpath):
    if not vcf_conf:
        err('No database for ' + dbname + ', skipping.')
        return None

    step_greetings('Annotating with ' + dbname)

    executable = get_java_tool_cmdline(cnf, 'snpsift')
    java = get_system_path(cnf, 'java')
    info('Java version:')
    call(cnf, java + ' -version')
    info()

    db_path = cnf['genome'].get(dbname)
    if not db_path:
        db_path = vcf_conf.get('path')
        if not db_path:
            err('Please, privide a path to ' + dbname + ' in the "genomes" section in the system config. The config is: ' + str(cnf['genome']))
            return
        verify_file(db_path, is_critical=True)

    annotations = vcf_conf.get('annotations')

    if not cnf.no_check:
        info('Removing previous annotations...')
        def delete_annos(rec):
            for anno in annotations:
                if anno in rec.INFO:
                    del rec.INFO[anno]
            return rec
        if annotations:
            input_fpath = iterate_vcf(cnf, input_fpath, delete_annos, suffix='d')

    anno_line = ''
    if annotations:
        anno_line = '-info ' + ','.join(annotations)

    cmdline = '{executable} annotate -v {anno_line} {db_path} {input_fpath}'.format(**locals())
    output_fpath = intermediate_fname(cnf, input_fpath, dbname)
    if output_fpath.endswith('.gz'):
        output_fpath = output_fpath[:-3]
    output_fpath = call_subprocess(cnf, cmdline, input_fpath, output_fpath,
        stdout_to_outputfile=True, exit_on_error=False)
    if not output_fpath:
        err('Error: snpsift resulted ' + str(output_fpath) + ' for ' + dbname)
        return output_fpath

    if not cnf.no_check:
        info_pattern = re.compile(r'''\#\#INFO=<
            ID=(?P<id>[^,]+),\s*
            Number=(?P<number>-?\d+|\.|[AG]),\s*
            Type=(?P<type>Integer|Float|Flag|Character|String),\s*
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)

        def _fix_after_snpsift(line, i, ctx):
            if not line.startswith('#'):
                if not ctx['met_CHROM']:
                    return None
                line = line.replace(' ', '_')
                assert ' ' not in line

            # elif line.startswith('##INFO=<ID=om'):
            #     line = line.replace(' ', '')

            elif not ctx['met_CHROM'] and line.startswith('#CHROM'):
                ctx['met_CHROM'] = True

            elif line.startswith('##INFO'):
                m = info_pattern.match(line)
                if m:
                    line = '##INFO=<ID={0},Number={1},Type={2},Description="{3}">'.format(
                        m.group('id'), m.group('number'), m.group('type'), m.group('desc'))
            return line

        output_fpath = iterate_file(cnf, output_fpath, _fix_after_snpsift, suffix='fx', ctx=dict(met_CHROM=False))

    return output_fpath


def _snpsift_db_nsfp(cnf, input_fpath):
    if 'dbnsfp' not in cnf.annotation or 'dbnsfp' not in cnf.genome:
        return None

    step_greetings('DB SNFP')

    executable = get_java_tool_cmdline(cnf, 'snpsift')

    db_path = cnf['genome']['dbnsfp']
    if not verify_file(db_path, 'DB NSFP file'):
        err('DB NSFP file is incorrect. Skipping.')
        return None

    annotations = cnf.annotation['dbnsfp'].get('annotations') or []

    # all_fields.extend(['dbNSFP_' + ann for ann in annotations])

    ann_line = ('-f ' + ','.join(annotations)) if annotations else ''

    cmdline = '{executable} dbnsfp {ann_line} -v -db {db_path} ' \
              '{input_fpath}'.format(**locals())
    output_fpath = intermediate_fname(cnf, input_fpath, 'db_nsfp')
    if output_fpath.endswith('.gz'):
        output_fpath = output_fpath[:-3]
    if call_subprocess(cnf, cmdline, input_fpath, output_fpath, stdout_to_outputfile=True, exit_on_error=False):
        return output_fpath
    else:
        return None


def _snpeff(cnf, input_fpath):
    if 'snpeff' not in cnf.annotation or 'snpeff' not in cnf.genome:
        return None, None, None

    step_greetings('SnpEff')

    snpeff = get_java_tool_cmdline(cnf, 'snpeff')

    stats_fpath = join(cnf.work_dir, cnf.sample + '-' + cnf.caller + '.snpEff_summary.csv')

    ref_name = cnf.genome.snpeff.reference or cnf.genome.name
    # if ref_name == 'GRCh37': ref_name += '.75'
    # if ref_name.startswith('hg38'): ref_name = 'GRCh38.78'

    opts = ''
    if cnf.annotation.snpeff.cancer: opts += ' -cancer'

    # db_path = cnf.genome.snpeff.reference
    # if db_path and isinstance(db_path, basestring): opts += ' -dataDir ' + db_path

    if cnf.annotation.snpeff.clinical_reporting or cnf.annotation.snpeff.canonical:
        opts += ' -canon '
    else:
        custom_transcripts = cnf.genome.snpeff.transcripts or cnf.snpeff_transcripts
        if custom_transcripts:
            if not verify_file(custom_transcripts, 'Transcripts for snpEff -onlyTr'):
                return None, None, None
            opts += ' -onlyTr ' + custom_transcripts + ' '

    if cnf.resources.snpeff.config:
        conf = get_system_path(cnf, cnf.resources.snpeff.config)
        if conf:
            opts += ' -c ' + conf + ' '
        else:
            err('Cannot find snpEff config file ' + str(cnf.resources.snpeff.config))

    if cnf.annotation.snpeff.extra_options:
        opts += ''

    if not cnf.no_check:
        info('Removing previous snpEff annotations...')
        res = remove_prev_eff_annotation(cnf, input_fpath)
        if not res:
            err('Could not remove preivous snpEff annotations')
            return None, None, None
        input_fpath = res

    snpeff_type = get_snpeff_type(snpeff)
    if snpeff_type == "old":
        opts += ' -stats ' + stats_fpath + ' -csvStats'
    else:
        opts += ' -csvStats ' + stats_fpath

    cmdline = ('{snpeff} eff {opts} -noLog -i vcf -o vcf {ref_name} '
               '{input_fpath}').format(**locals())

    output_fpath = intermediate_fname(cnf, input_fpath, 'snpEff')
    if output_fpath.endswith('.gz'):
        output_fpath = output_fpath[:-3]
    res = call_subprocess(cnf, cmdline, input_fpath, output_fpath,
                          exit_on_error=False, stdout_to_outputfile=True)
    if res:
        return output_fpath, stats_fpath, splitext(stats_fpath)[0] + '.genes.txt'
    else:
        return None, None, None


def _tracks(cnf, track_fpath, input_fpath):
    if not verify_file(track_fpath):
        return None

    field_name = splitext_plus(basename(track_fpath))[0]

    step_greetings('Intersecting with ' + field_name)

    toolpath = get_system_path(cnf, 'vcfannotate')
    if not toolpath:
        err('WARNING: Skipping annotation with tracks: vcfannotate '
            'executable not found, you probably need to specify path in system_config, or '
            'run load bcbio:  . /group/ngs/bin/bcbio-prod.sh"')
        return None

    # self.all_fields.append(field_name)

    cmdline = '{toolpath} -b {track_fpath} -k {field_name} {input_fpath}'.format(**locals())

    assert input_fpath
    output_fpath = intermediate_fname(cnf, input_fpath, field_name)
    if output_fpath.endswith('.gz'):
        output_fpath = output_fpath[:-3]
    output_fpath = call_subprocess(cnf, cmdline, input_fpath, output_fpath,
                                   stdout_to_outputfile=True)
    if not output_fpath:
        err('Error: tracks resulted ' + str(output_fpath) + ' for ' + track_fpath)
        return output_fpath

    # Set TRUE or FALSE for tracks
    def proc_line(line, i):
        if field_name in line:
            if not line.startswith('#'):
                fields = line.split('\t')
                info_line = fields[7]
                info_pairs = [attr.split('=') for attr in info_line.split(';')]
                info_pairs = [[pair[0], ('TRUE' if pair[1] else 'FALSE')]
                              if pair[0] == field_name and len(pair) > 1
                              else pair for pair in info_pairs]
                info_line = ';'.join('='.join(pair) if len(pair) == 2
                                     else pair[0] for pair in info_pairs)
                fields = fields[:7] + [info_line] + fields[8:]
                return '\t'.join(fields)
        return line

    assert output_fpath
    return iterate_file(cnf, output_fpath, proc_line, suffix='trk')


def _add_annotation(cnf, input_fpath, key, value, number, type_, description):
    step_greetings('Adding annotation...')
    def proc_rec(rec):
        rec.INFO[key] = value
        return rec
    output_fpath = iterate_vcf(cnf, input_fpath, proc_rec)

    info('Adding header meta info...')
    def _add_format_header(l, i):
        if l.startswith('#CHROM'):
            ext_l = ''
            ext_l += '##INFO=<ID={key},Number={number},Type={type_},Description="{desc}">\n'.format(
                key=key, number=number, type_=type_, desc=description)
            return ext_l + l
        return l
    output_fpath = iterate_file(cnf, output_fpath, _add_format_header)
    return output_fpath


def _filter_malformed_fields(cnf, input_fpath):
    step_greetings('Correcting malformed fields...')

    def proc_rec(rec):
        for k, v in rec.INFO.items():
            if isinstance(v, list):
                if v[-1] == '.':
                    rec.INFO[k] = rec.INFO[k][:-1]
                if v[0] == '.':
                    rec.INFO[k] = rec.INFO[k][1:]
        return rec

    def proc_line(line, i):
        if line.startswith('#'):
            return line.replace("\' \">", "\'\">")  # For vcf-merge
        return line

        # else:
            # if ',.' in line or '.,' in line:
            #     fields = line.split('\t')
            #     info_line = fields[7]
            #     info_pairs = [attr.split('=') for attr in info_line.split(';')]
            #     new_info_pairs = []
            #     for p in info_pairs:
            #         if len(p) == 2:
            #             if p[1].endswith(',.'):
            #                 p[1] = p[1][:-2]
            #             if p[1].startswith('.,'):
            #                 p[1] = p[1][2:]
            #             new_info_pairs.append('='.join(p))
            #     info_line = ';'.join(new_info_pairs)
            #     fields = fields[:7] + [info_line] + fields[8:]
            #     return '\t'.join(fields)

    info('Correcting INFO fields...')
    output_fpath = iterate_vcf(cnf, input_fpath, proc_rec, suffix='corr')
    info('')
    info('Correcting headers for vcf-merge...')
    output_fpath = iterate_file(cnf, output_fpath, proc_line, suffix='corr_headr')

    return output_fpath