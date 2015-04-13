import shutil
import os
from os.path import splitext, basename, join, isfile
import socket
import sys

from source.calling_process import call_subprocess, call
from source.file_utils import iterate_file, intermediate_fname, verify_file
from source.logger import step_greetings, critical, info, err, warn
from source.tools_from_cnf import get_system_path, get_java_tool_cmdline
from source.file_utils import file_exists, code_base_path
from source.variants.tsv import make_tsv
from source.variants.vcf_processing import iterate_vcf, remove_prev_eff_annotation, bgzip_and_tabix, igvtools_index


def run_annotators(cnf, vcf_fpath, bam_fpath):
    annotated = False
    original_vcf = cnf.vcf

    def delete_ids(rec):  # deleting existing dbsnp and cosmic ID annotations
        if rec.ID:
            if isinstance(rec.ID, basestring):
                if rec.ID.startswith('rs') or rec.ID.startswith('COS'):
                    rec.ID = None
            else:
                rec.ID = [id for id in rec.ID if not id.startswith('rs') and not id.startswith('COS')]

        if not rec.FILTER:
            rec.FILTER = 'PASS'

        return rec

    vcf_fpath = iterate_vcf(cnf, vcf_fpath, delete_ids, suffix='delID')

    dbs = [(dbname, cnf.annotation[dbname])
           for dbname in ['dbsnp', 'clinvar', 'cosmic', 'oncomine']
           if dbname in cnf.annotation]

    if 'custom_vcfs' in cnf.annotation:
        dbs.extend(cnf.annotation['custom_vcfs'].items())

    for dbname, dbconf in dbs:
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

    if cnf.annotation.get('tracks'):
        for track in cnf.annotation['tracks']:
            res = _tracks(cnf, track, vcf_fpath)
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


def finialize_annotate_file(cnf, vcf_fpath, samplename, callername):
    # vcf_fpath = leave_first_sample(cnf, vcf_fpath)

    if not cnf.get('no_correct_vcf'):
        vcf_fpath = _filter_malformed_fields(cnf, vcf_fpath)

    info('Adding SAMPLE=' + samplename + ' annotation...')
    vcf_fpath = _add_annotation(cnf, vcf_fpath, 'SAMPLE', samplename)

    final_vcf_fname = samplename + '-' + callername + '.anno.vcf'
    final_tsv_fname = samplename + '-' + callername + '.anno.tsv'
    final_vcf_fpath = join(cnf.output_dir, final_vcf_fname)
    final_tsv_fpath = join(cnf.output_dir, final_tsv_fname)

    # Moving final VCF
    if isfile(final_vcf_fpath):
        os.remove(final_vcf_fpath)
    if isfile(final_vcf_fpath + '.gz'):
        os.remove(final_vcf_fpath + '.gz')
    shutil.copy(vcf_fpath, final_vcf_fpath)

    # Converting to TSV
    if 'tsv_fields' in cnf.annotation:
        tsv_fpath = make_tsv(cnf, vcf_fpath, samplename)
        if not tsv_fpath:
            err('TSV convertion didn\'t work')
        else:
            if isfile(final_tsv_fpath):
                os.remove(final_tsv_fpath)
            shutil.copy(tsv_fpath, final_tsv_fpath)
    else:
        final_tsv_fpath = None

    # Indexing
    info()
    info('Compressing and indexing with bgzip+tabix ' + final_vcf_fpath)
    final_vcf_fpath = bgzip_and_tabix(cnf, final_vcf_fpath)

    return final_vcf_fpath, final_tsv_fpath


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

    info('Removing previous annotations...')
    annotations = vcf_conf.get('annotations')

    def delete_annos(rec):
        for anno in annotations:
            if anno in rec.INFO:
                del rec.INFO[anno]
        return rec
    if annotations:
        input_fpath = iterate_vcf(cnf, input_fpath, delete_annos, suffix='delANNO')

    executable = get_java_tool_cmdline(cnf, 'snpsift')
    java = get_system_path(cnf, 'java')
    info('Java version:')
    call(cnf, java + ' -version')
    info()

    db_path = cnf['genome'].get(dbname)
    if not db_path:
        db_path = vcf_conf.get('path')
        if not db_path:
            err('Please, privide a path to ' + dbname + ' in the run config '
                '("path:" field), or in the "genomes" section in the system config')
            return
        verify_file(db_path, is_critical=True)

    anno_line = ''
    if annotations:
        anno_line = '-info ' + ','.join(annotations)

    cmdline = '{executable} annotate -v {anno_line} {db_path} {input_fpath}'.format(**locals())
    output_fpath = intermediate_fname(cnf, input_fpath, dbname)
    output_fpath = call_subprocess(cnf, cmdline, input_fpath, output_fpath,
                                   stdout_to_outputfile=True, exit_on_error=False)
    if not output_fpath:
        err('Error: snpsift resulted ' + str(output_fpath) + ' for ' + dbname)
        return output_fpath

    # all_fields.extend(annotations)

    def proc_line(line, i):
        if not line.startswith('#'):
            line = line.replace(' ', '_')
            assert ' ' not in line
        return line

    output_fpath = iterate_file(cnf, output_fpath, proc_line, suffix='f')
    return output_fpath


def _snpsift_db_nsfp(cnf, input_fpath):
    if 'dbnsfp' not in cnf:
        return None

    step_greetings('DB SNFP')

    executable = get_java_tool_cmdline(cnf, 'snpsift')

    db_path = cnf['genome'].get('dbnsfp')
    if not db_path:
        critical('Please, provide a path to DB NSFP file in '
                 'the "genomes" section in the system config.')

    annotations = cnf.annotation['dbnsfp'].get('annotations') or []

    # all_fields.extend(['dbNSFP_' + ann for ann in annotations])

    ann_line = ('-f ' + ','.join(annotations)) if annotations else ''

    cmdline = '{executable} dbnsfp {ann_line} -v -db {db_path} ' \
              '{input_fpath}'.format(**locals())
    output_fpath = intermediate_fname(cnf, input_fpath, 'db_nsfp')
    if call_subprocess(cnf, cmdline, input_fpath, output_fpath,
                       stdout_to_outputfile=True, exit_on_error=False):
        return output_fpath
    else:
        return None


def _snpeff(cnf, input_fpath):
    step_greetings('SnpEff')

    info('Removing previous EFF annotations...')

    res = remove_prev_eff_annotation(cnf, input_fpath)
    if res:
        input_fpath = res
    info('')

    if 'snpeff' not in cnf.annotation:
        return None, None, None

    snpeff = get_java_tool_cmdline(cnf, 'snpeff')

    stats_fpath = join(cnf.work_dir, cnf.name + '-' + cnf.caller + '.snpEff_summary.html')

    ref_name = cnf.genome.name
    if ref_name == 'GRCh37': ref_name += '.75'

    opts = ''
    if cnf.annotation.snpeff.cancer: opts += ' -cancer'

    db_path = cnf.genome.snpeff
    if db_path: opts += ' -dataDir ' + db_path

    if cnf.annotation.snpeff.clinical_reporting or cnf.annotation.snpeff.canonical:
        opts += ' -canon '
    else:
        custom_transcripts = cnf.genome.snpeff_transcripts
        if custom_transcripts:
            if not verify_file(custom_transcripts, 'Transcripts for snpEff -onlyTr'):
                return None, None, None
            opts += ' -onlyTr ' + custom_transcripts + ' '

    if cnf.resources.snpeff.config:
        opts += ' -c ' + cnf.resources.snpeff.config + ' '

    if cnf.annotation.snpeff.extra_options:
        opts += ''

    cmdline = ('{snpeff} eff {opts} -stats {stats_fpath} -formatEff '
               '-csvStats -noLog -i vcf -o vcf {ref_name} '
               '{input_fpath}').format(**locals())

    output_fpath = intermediate_fname(cnf, input_fpath, 'snpEff')
    res = call_subprocess(cnf, cmdline, input_fpath, output_fpath,
                          exit_on_error=False, stdout_to_outputfile=True)
    if res:
        return output_fpath, stats_fpath, splitext(stats_fpath)[0] + '.genes.txt'
    else:
        return None, None, None


def _tracks(cnf, track_path, input_fpath):
    field_name = splitext(basename(track_path))[0]

    step_greetings('Intersecting with ' + field_name)

    toolpath = get_system_path(cnf, 'vcfannotate')
    if not toolpath:
        err('WARNING: Skipping annotation with tracks: vcfannotate '
            'executable not found, you probably need to specify path in system_config, or '
            'run load bcbio:  . /group/ngs/bin/bcbio-prod.sh"')
        return

    # self.all_fields.append(field_name)

    cmdline = '{toolpath} -b {track_path} -k {field_name} {input_fpath}'.format(**locals())

    assert input_fpath
    output_fpath = intermediate_fname(cnf, input_fpath, field_name)
    output_fpath = call_subprocess(cnf, cmdline, input_fpath, output_fpath,
                                   stdout_to_outputfile=True)
    if not output_fpath:
        err('Error: tracks resulted ' + str(output_fpath) + ' for ' + track_path)
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
    return iterate_file(cnf, output_fpath, proc_line, 'trk')


def _add_annotation(cnf, input_fpath, key, value):
    step_greetings('Adding annotation...')

    def proc_rec(rec):
        rec.INFO[key] = value
        return rec

    output_fpath = iterate_vcf(cnf, input_fpath, proc_rec, 'plus')
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
    output_fpath = iterate_vcf(cnf, input_fpath, proc_rec, 'corr')
    info('')
    info('Correcting headers for vcf-merge...')
    output_fpath = iterate_file(cnf, output_fpath, proc_line, 'corr_headr')

    return output_fpath