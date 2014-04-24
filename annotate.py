#!/usr/bin/env python
import sys
import subprocess
import shutil
import os
from os.path import join, splitext, basename, realpath, isdir, isfile, dirname
from collections import OrderedDict
import tempfile
from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from src.utils import which, splitext_plus, add_suffix, file_exists
from src.transaction import file_transaction
from src.my_utils import get_gatk_type, err, critical, info, \
    intermediate_fname, verify_file, safe_mkdir, remove_quotes, iterate_file, verify_dir


def annotate(samples, parallel=False):
    if len(samples) == 1:
        sample_name, sample_cnf = samples.items()[0]
        annotate_one(sample_name, sample_cnf)
    else:
        results = []
        if parallel:
            try:
                from joblib import Parallel, delayed
            except ImportError:
                critical(
                    '\nERROR: Joblib not found. You may want samples to be processed '
                    'in parallel, in this case, make sure python joblib intalled. '
                    '(pip install joblib).')
            else:
                for sample_name, sample_cnf in samples.items():
                    sample_cnf['verbose'] = False

                results = Parallel(n_jobs=len(samples)) \
                    (delayed(annotate_one) \
                            (sample_name, sample_cnf, multiple_samples=True)
                        for sample_name, sample_cnf in samples.items())
        else:
            for sample_name, sample_cnf in samples.items():
                results.append(
                    annotate_one(sample_name, sample_cnf, multiple_samples=True))

        info('')
        info('*' * 70)
        info('Results for each samples:')
        for (sample_name, cnf), (vcf, tsv, qc_txt, qc_plots) in zip(samples.items(), results):
            info(cnf['log'], sample_name + ':')
            info(cnf['log'], '  ' + vcf)
            info(cnf['log'], '  ' + tsv)
            if qc_txt:
                info(cnf['log'], '  ' + qc_txt)
            if qc_plots:
                info(cnf['log'], '  ' + qc_plots)

    for name, data in samples.items():
        work_dirpath = data['work_dir']
        tx_dirpath = join(work_dirpath, 'tx')

        if isdir(tx_dirpath):
            shutil.rmtree(tx_dirpath)

        if not data.get('keep_intermediate') \
                and isdir(work_dirpath):
            shutil.rmtree(work_dirpath)


def _parse_gatk_report(report_filename, databases, novelty, metrics):
    database_col_name = 'CompRod'
    database_col_id = None
    novelty_col_name = 'Novelty'
    novelty_col_id = None

    report = dict()
    comments_section = False
    cur_header = []
    cur_metrics_ids = []
    for line in open(report_filename):
        if not line.strip():
            continue
        if line.startswith('#'): # comment line
            comments_section = True
            continue
        elif comments_section:
            comments_section = False
            cur_header = line.split()
            cur_metrics_ids = []
            database_col_id = cur_header.index(database_col_name)
            novelty_col_id = cur_header.index(novelty_col_name)
            for metric in metrics:
                if metric in cur_header:
                    cur_metrics_ids.append(cur_header.index(metric))
                    if metric not in report:
                        report[metric] = dict()
        elif cur_metrics_ids:  # process lines only if there are metrics in current section
            values = line.split()
            cur_database = values[database_col_id]
            cur_novelty = values[novelty_col_id]
            if (cur_database not in databases) or (cur_novelty not in novelty):
                continue
            for metric_id in cur_metrics_ids:
                if cur_database not in report[cur_header[metric_id]]:
                    report[cur_header[metric_id]][cur_database] = dict()
                report[cur_header[metric_id]][cur_database][cur_novelty] = values[metric_id]
    return report


def _make_final_report(report_dict, report_filename, sample_name,
                        databases, novelty, metrics):
    header = ['Metric', 'Novelty'] + databases + ['Average']
    full_report = [header]
    for cur_metric in metrics:
        for cur_novelty in novelty:
            cur_row = [cur_metric, cur_novelty]
            sum = 0.0
            for cur_database in databases:
                if cur_metric == 'variantRatePerBp': # confusing name and value format
                    cur_row[0] = 'basesPerVariant'
                    cur_row.append("%.2f" % float(report_dict[cur_metric][cur_database][cur_novelty]))
                else:
                    cur_row.append(report_dict[cur_metric][cur_database][cur_novelty])
                sum += float(cur_row[-1])
            average = sum / len(databases)
            cur_row.append("%.2f" % average)
            full_report.append(cur_row)

    col_widths = [0] * len(header)
    for row in full_report:
        for id, value in enumerate(row):
            col_widths[id] = max(len(value), col_widths[id])

    out = open(report_filename, 'w')
    out.write('Sample name: ' + sample_name + '\n\n')
    for row in full_report:
        out.write('  '.join('%-*s' % (col_width, value) for col_width, value in zip(col_widths, row)) + "\r\n")
    out.close()


def _get_plots_from_bcftools(bcftools_report_dir, output_dir):
    original_plots_names = ['indels.0.png', 'substitutions.0.png']
    final_plots_names = ['indels.png', 'substitution.png']
    for id, original_plot in enumerate(original_plots_names):
        plot_src_filename = join(bcftools_report_dir, original_plot)
        plot_dst_filename = join(output_dir, final_plots_names[id])
        if os.path.exists(plot_src_filename):
            shutil.copyfile(plot_src_filename, plot_dst_filename)


def bcftools_qc(cnf, vcf_fpath):
    step_greetings(cnf, 'Quality control plots')

    work_dir = cnf['work_dir']

    gzipped_fpath = join(work_dir, os.path.basename(vcf_fpath) + '.gz')
    cmdline = 'bgzip %s -c ' % vcf_fpath
    with file_transaction(gzipped_fpath) as tx:
        subprocess.call(cmdline.split(), stdout=tx)

    tbi_fpath = gzipped_fpath + '.tbi'
    cmdline = 'tabix -p vcf %s ' % gzipped_fpath
    with file_transaction(tbi_fpath) as tx:
        subprocess.call(cmdline.split(), stdout=tx)

    bcftools = _get_tool_cmdline(cnf, 'bcftools')
    plot_vcfstats = _get_tool_cmdline(cnf, 'plot-vcfstats')
    text_report_fpath = join(work_dir, 'bcftools.report')
    viz_report_dir = join(work_dir, 'viz_report/')

    cmdline = '%s stats %s ' % (bcftools, gzipped_fpath)
    with file_transaction(text_report_fpath) as tx:
        subprocess.call(cmdline.split(), stdout=tx)

    cmdline = '%s -s %s -p %s --no-PDF ' % (plot_vcfstats, text_report_fpath, viz_report_dir)
    subprocess.call(cmdline.split())

    _get_plots_from_bcftools(viz_report_dir, cnf['output_dir'])
    return viz_report_dir


def gatk_qc(cnf, vcf_fpath):
    step_greetings(cnf, 'Quality control reports')

    log = cnf['log']
    work_dir = cnf['work_dir']

    qc_cnf = cnf['quality_control']
    databases = qc_cnf.get('database_vcfs')
    novelty = qc_cnf.get('novelty')
    metrics = qc_cnf.get('metrics')

    executable = _get_java_tool_cmdline(cnf, 'gatk')
    ref_fpath = cnf['genome']['seq']
    report_fpath = join(work_dir, cnf['name'] + '_gatk.report')

    cmdline = ('{executable} -nt 20 -R {ref_fpath} -T VariantEval'
               ' --eval:tmp {vcf_fpath} -o {report_fpath}').format(**locals())

    if 'dbsnp' in databases:
        cmdline += ' -D ' + databases['dbsnp']
        del databases['dbsnp']
    for db_name, db_path in databases.items():
        cmdline += ' -comp:' + db_name + ' ' + db_path

    _call_and_remove(
        cnf, cmdline, None, report_fpath, stdout_to_outputfile=False,
        to_remove=[vcf_fpath + '.idx'],
        keep_original_if_not_keep_intermediate=True)

    report = _parse_gatk_report(report_fpath, databases.keys(), novelty, metrics)

    final_report_fpath = join(cnf['output_dir'], cnf['name'] + '.qc.txt')

    _make_final_report(report, final_report_fpath, cnf['name'],
                       databases.keys(), novelty, metrics)
    return final_report_fpath



def annotate_one(sample_name, cnf, multiple_samples=False):
    if cnf.get('keep_intermediate'):
        cnf['log'] = join(cnf['work_dir'], sample_name + '_log.txt')
        if isfile(cnf['log']):
            os.remove(cnf['log'])

        with open(join(cnf['work_dir'], sample_name + '_config.yaml'), 'w') as f:
            f.write(dump(cnf, Dumper=Dumper))
    else:
        cnf['log'] = None

    if multiple_samples:
        info('')
        info('*' * 70)
        msg = '*' * 3 + ' Sample ' + sample_name + ' '
        info(cnf.get('log'), msg + ('*' * (70 - len(msg)) if len(msg) < 70 else ''))
        info(cnf.get('log'), 'VCF: ' + cnf['vcf'])
        if cnf.get('bam'):
            info(cnf.get('log'), 'BAM: ' + cnf['bam'])

    vcf_fpath = original_vcf_fpath = cnf['vcf']
    # sample['fields'] = []
    work_dir = cnf['work_dir']

    if cnf.get('split_genotypes'):
        vcf_fpath = split_genotypes(cnf, vcf_fpath)

    if cnf.get('ensemble'):
        vcf_fpath = filter_ensemble(cnf, vcf_fpath)

    if 'gatk' in cnf:
        vcf_fpath = gatk(cnf, vcf_fpath, cnf.get('bam'), work_dir)

    if 'dbsnp' in cnf:
        vcf_fpath = snpsift_annotate(cnf, cnf['dbsnp'],
                                     'dbsnp', vcf_fpath, work_dir)
    if 'cosmic' in cnf:
        vcf_fpath = snpsift_annotate(cnf, cnf['cosmic'],
                                     'cosmic', vcf_fpath, work_dir)
    if 'custom_vcfs' in cnf:
        for dbname, vcf_conf in cnf['custom_vcfs'].items():
            vcf_fpath = snpsift_annotate(cnf, vcf_conf, dbname,
                                         vcf_fpath, work_dir)

    if 'dbnsfp' in cnf:
        vcf_fpath = snpsift_db_nsfp(cnf, vcf_fpath, work_dir)

    if 'snpeff' in cnf:
        remove_annotation(cnf, 'EFF', vcf_fpath, work_dir)
        vcf_fpath = snpeff(cnf, vcf_fpath, work_dir)

    if cnf.get('tracks'):
        for track in cnf['tracks']:
            vcf_fpath = tracks(cnf, track, vcf_fpath, work_dir)

    vcf_fpath = filter_fields(cnf, vcf_fpath, work_dir)

    # Copying final VCF
    final_fname = add_suffix(basename(original_vcf_fpath), 'anno')
    final_vcf_fpath = join(cnf['output_dir'], final_fname)
    if isfile(final_vcf_fpath):
        os.remove(final_vcf_fpath)
    shutil.copyfile(vcf_fpath, final_vcf_fpath)

    # making TSV
    tsv_fpath = extract_fields(cnf, vcf_fpath, work_dir, sample_name)

    manual_tsv_fields = cnf.get('tsv_fields')
    if manual_tsv_fields:
        field_map = dict((rec.keys()[0], rec.values()[0]) for rec in manual_tsv_fields)
        if cnf.get('keep_intermediate'):
            info(cnf['log'], 'Saved TSV file to ' + tsv_fpath)
        tsv_fpath = rename_fields(cnf, tsv_fpath, field_map, work_dir)
        if cnf.get('keep_intermediate'):
            info(cnf['log'], 'Saved TSV file with nice names to ' + tsv_fpath)

    # Copying final TSV
    final_tsv_fpath = splitext_plus(final_vcf_fpath)[0] + '.tsv'
    if isfile(final_tsv_fpath):
        os.remove(final_tsv_fpath)
    shutil.copyfile(tsv_fpath, final_tsv_fpath)

    txt_report_fpath = None
    plots_dirpath = None
    if 'quality_control' in cnf:
        txt_report_fpath = gatk_qc(cnf, final_vcf_fpath)
        plots_dirpath = bcftools_qc(cnf, final_vcf_fpath)
        info(cnf['log'], 'Saved quality control text report to ' + txt_report_fpath)
        info(cnf['log'], 'Saved quality control plots to ' + plots_dirpath)

    info(cnf['log'], '')
    info(cnf['log'], 'Saved final VCF to ' + final_vcf_fpath)
    info(cnf['log'], 'Saved final TSV to ' + final_tsv_fpath)

    return final_vcf_fpath, final_tsv_fpath, txt_report_fpath, plots_dirpath


def remove_annotation(cnf, field_to_del, input_fpath, work_dir):
    def proc_line(l):
        if field_to_del in l:
            if l.startswith('##INFO='):
                try:
                    if l.split('=', 1)[1].split(',', 1)[0].split('=')[1] == field_to_del:
                        return None
                except IndexError:
                    critical(cnf['log'], 'Incorrect VCF at line: ' + l)
            elif not l.startswith('#'):
                fields = l.split('\t')
                info_line = fields[7]
                info_pairs = [attr.split('=') for attr in info_line.split(';')]
                info_pairs = filter(lambda pair: pair[0] != field_to_del, info_pairs)
                info_line = ';'.join('='.join(pair) if len(pair) == 2
                                     else pair[0] for pair in info_pairs)
                fields = fields[:7] + [info_line] + fields[8:]
                return '\t'.join(fields)
        return l
    return iterate_file(cnf, input_fpath, proc_line, work_dir)


def extract_sample(cnf, input_fpath, samplename, work_dir):
    step_greetings(cnf, 'Separating out sample ' + samplename)

    executable = _get_java_tool_cmdline(cnf, 'gatk')
    ref_fpath = cnf['genome']['seq']

    corr_samplename = ''.join([c if c.isalnum() else '_' for c in samplename])
    output_fpath = intermediate_fname(work_dir, input_fpath, suf=corr_samplename)

    cmd = '{executable} -nt 30 -R {ref_fpath} -T SelectVariants ' \
          '--variant {input_fpath} -sn {samplename} -o {output_fpath}'.format(**locals())
    _call_and_remove(cnf, cmd, input_fpath, output_fpath,
                     stdout_to_outputfile=False,
                     keep_original_if_not_keep_intermediate=True)
    return output_fpath


def _call_and_remove(cnf, cmdline, input_fpath, output_fpath,
                     stdout_to_outputfile=True, to_remove=None,
                     keep_original_if_not_keep_intermediate=False):
    to_remove = to_remove or []

    if cnf.get('keep_intermediate') and cnf.get('reuse_intermediate'):
        if file_exists(output_fpath):
            info(cnf.get('log'), output_fpath + ' exists, reusing')
            return output_fpath

    err_fpath = None
    if cnf.get('work_dir'):
        _, err_fpath = tempfile.mkstemp(dir=cnf.get('work_dir'), prefix='err_tmp')
        to_remove.append(err_fpath)

    with file_transaction(output_fpath) as tx_out_fpath:
        if stdout_to_outputfile:
            info(cnf.get('log'), cmdline + ' > ' + tx_out_fpath)
        else:
            cmdline = cmdline.replace(output_fpath, tx_out_fpath)
            info(cnf.get('log'), cmdline)

        if cnf['verbose']:
            proc = subprocess.Popen(
                cmdline.split(),
                stdout=open(tx_out_fpath, 'w') if stdout_to_outputfile else subprocess.PIPE,
                stderr=subprocess.STDOUT if not stdout_to_outputfile else subprocess.PIPE)

            if proc.stdout:
                for line in iter(proc.stdout.readline, ''):
                    info(cnf.get('log'), '   ' + line.strip())
            elif proc.stderr:
                for line in iter(proc.stderr.readline, ''):
                    info(cnf.get('log'), '   ' + line.strip())

            ret_code = proc.wait()
            if ret_code != 0:
                for fpath in to_remove:
                    if fpath and isfile(fpath):
                        os.remove(fpath)
                critical(cnf.get('log'), 'Command returned status ' + str(ret_code) +
                         ('. Log in ' + cnf['log'] if 'log' in cnf else '.'))
        else:
            res = subprocess.call(
                cmdline.split(),
                stdout=open(tx_out_fpath, 'w') if stdout_to_outputfile else
                    (open(err_fpath, 'a') if err_fpath else None),
                stderr=open(err_fpath, 'a') if err_fpath else None)
            if res != 0:
                with open(err_fpath) as err_f:
                    info(cnf.get('log'), '')
                    info(cnf.get('log'), err_f.read())
                    info(cnf.get('log'), '')
                for fpath in to_remove:
                    if fpath and isfile(fpath):
                        os.remove(fpath)
                critical(cnf.get('log'), 'Command returned status ' + str(res) +
                         ('. Log in ' + cnf['log'] if 'log' in cnf else '.'))
            else:
                if 'log' in cnf and err_fpath:
                    with open(err_fpath, 'a') as err_f, \
                         open(cnf.get('log'), 'a') as log_f:
                        log_f.write('')
                        log_f.write(err_f.read())
                        log_f.write('')

    for fpath in to_remove:
        if fpath and isfile(fpath):
            os.remove(fpath)

    if (not cnf.get('keep_intermediate') and
        not keep_original_if_not_keep_intermediate and
            input_fpath):
        os.remove(input_fpath)
    info(cnf.get('log'), 'Saved to ' + output_fpath)
    return output_fpath


def _get_java_tool_cmdline(cnf, name):
    cmdline_template = _get_script_cmdline_template(cnf, 'java', name)
    jvm_opts = cnf['resources'][name].get('jvm_opts', []) + ['']
    return cmdline_template % (' '.join(jvm_opts) + ' -jar')


def _get_script_cmdline_template(cnf, executable, script_name):
    if not which(executable):
        exit(executable + ' executable required, maybe you need '
             'to run "module load ' + executable + '"?')
    if 'resources' not in cnf:
        critical(cnf['log'], 'System config yaml must contain resources section with '
                 + script_name + ' path.')
    if script_name not in cnf['resources']:
        critical(cnf['log'], 'System config resources section must contain '
                 + script_name + ' info (with a path to the tool).')
    tool_config = cnf['resources'][script_name]
    if 'path' not in tool_config:
        critical(script_name + ' section in the system config must contain a path to the tool.')
    tool_path = tool_config['path']
    if not verify_file(tool_path, script_name):
        exit(1)
    return executable + ' %s ' + tool_path


def snpsift_annotate(cnf, vcf_conf, dbname, input_fpath, work_dir):
    step_greetings(cnf, 'Annotate with ' + dbname)

    executable = _get_java_tool_cmdline(cnf, 'snpsift')

    db_path = cnf['genome'].get(dbname)
    if not db_path:
        db_path = vcf_conf.get('path')
        if not db_path:
            critical(cnf['log'], 'Please, privide a path to ' + dbname + ' in the run config '
                     '("path:" field), or in the "genomes" section in the system config')
        if not verify_file(db_path):
            exit()

    annotations = cnf[dbname].get('annotations')
    # all_fields.extend(annotations)
    anno_line = ('-info ' + ','.join(annotations)) if annotations else ''
    cmdline = '{executable} annotate -v {anno_line} {db_path} {input_fpath}'.format(**locals())
    output_fpath = intermediate_fname(work_dir, input_fpath, dbname)
    output_fpath = _call_and_remove(cnf, cmdline, input_fpath, output_fpath,
                                    stdout_to_outputfile=True)
    def proc_line(line):
        if not line.startswith('#'):
            line = line.replace(' ', '_')
            assert ' ' not in line
        return line
    output_fpath = iterate_file(cnf, output_fpath, proc_line, work_dir)
    return output_fpath


def snpsift_db_nsfp(cnf, input_fpath, work_dir):
    if 'dbnsfp' not in cnf:
        return input_fpath

    step_greetings(cnf, 'DB SNFP')

    executable = _get_java_tool_cmdline(cnf, 'snpsift')

    db_path = cnf['genomes'].get('dbnsfp')
    if not db_path:
        critical(cnf['log'], 'Please, provide a path to DB NSFP file in '
                 'the "genomes" section in the system config.')

    annotations = cnf['dbnsfp'].get('annotations', [])
    # self.all_fields.extend(['dbNSFP_' + ann for ann in annotations])
    ann_line = ('-f ' + ','.join(annotations)) if annotations else ''

    cmdline = '{executable} dbnsfp {ann_line} -v {db_path} {input_fpath}'.format(**locals())
    output_fpath = intermediate_fname(work_dir, input_fpath, 'db_nsfp')
    return _call_and_remove(cnf, cmdline, input_fpath, output_fpath,
                            stdout_to_outputfile=True)


def snpeff(cnf, input_fpath, work_dir):
    if 'snpeff' not in cnf:
        return input_fpath

    step_greetings(cnf, 'SnpEff')

    # self.all_fields.extend([
    #     "EFF[*].EFFECT", "EFF[*].IMPACT", "EFF[*].FUNCLASS", "EFF[*].CODON",
    #     "EFF[*].AA", "EFF[*].AA_LEN", "EFF[*].GENE", "EFF[*].CODING",
    #     "EFF[*].TRID", "EFF[*].RANK"])

    executable = _get_java_tool_cmdline(cnf, 'snpeff')
    ref_name = cnf['genome']['name']
    db_path = cnf['genome'].get('snpeff')
    if not db_path:
        critical(cnf['log'], 'Please, provide a path to SnpEff data in '
                 'the "genomes" section in the system config.')

    cmdline = ('{executable} eff -dataDir {db_path} -noStats -noLog -1 '
               '-i vcf -o vcf {ref_name} {input_fpath}').format(**locals())

    if cnf['snpeff'].get('clinical_reporting') or \
            cnf['snpeff'].get('canonical'):
        cmdline += ' -canon -hgvs '

    if cnf['snpeff'].get('cancer'):
        cmdline += ' -cancer '

    output_fpath = intermediate_fname(work_dir, input_fpath, 'snpEff')
    return _call_and_remove(cnf, cmdline, input_fpath, output_fpath,
                            stdout_to_outputfile=True)


def tracks(cnf, track_path, input_fpath, work_dir):
    field_name = splitext(basename(track_path))[0]

    step_greetings(cnf, 'Intersecting with ' + field_name)

    toolpath = _get_tool_cmdline(cnf, 'vcfannotate')
    if not toolpath:
        err(cnf['log'], 'WARNING: Skipping annotation with tracks: vcfannotate '
            'executable not found, you probably need to '
            'run the commandline:  . /group/ngs/bin/bcbio-prod.sh"')
        return

    # self.all_fields.append(field_name)

    cmdline = 'vcfannotate -b {track_path} -k {field_name} {input_fpath}'.format(**locals())

    output_fpath = intermediate_fname(work_dir, input_fpath, field_name)
    output_fpath = _call_and_remove(cnf, cmdline, input_fpath, output_fpath,
                                    stdout_to_outputfile=True)

    # Set TRUE or FALSE for tracks
    def proc_line(line):
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
    return iterate_file(cnf, output_fpath, proc_line, work_dir)


def gatk(cnf, input_fpath, bam_fpath, work_dir):
    if 'gatk' not in cnf:
        return input_fpath

    step_greetings(cnf, 'GATK')

    executable = _get_java_tool_cmdline(cnf, 'gatk')

    output_fpath = intermediate_fname(work_dir, input_fpath, 'gatk')

    ref_fpath = cnf['genome']['seq']

    cmdline = ('{executable} -nt 20 -R {ref_fpath} -T VariantAnnotator'
               ' --variant {input_fpath} -o {output_fpath}').format(**locals())
    if bam_fpath:
        cmdline += ' -I ' + bam_fpath

    gatk_annos_dict = {
        'Coverage': 'DP',
        'BaseQualityRankSumTest': 'BaseQRankSum',
        'FisherStrand': 'FS',
        'GCContent': 'GC',
        'HaplotypeScore': 'HaplotypeScore',
        'HomopolymerRun': 'HRun',
        'RMSMappingQuality': 'MQ',
        'MappingQualityRankSumTest': 'MQRankSum',
        'MappingQualityZero': 'MQ0',
        'QualByDepth': 'QD',
        'ReadPosRankSumTest': 'ReadPosRankSum'
    }
    annotations = cnf['gatk'].get('annotations', [])

    # self.all_fields.extend(gatk_annos_dict.get(ann) for ann in annotations)
    gatk_type = get_gatk_type(_get_java_tool_cmdline(cnf, 'gatk'))
    for ann in annotations:
        if ann == 'DepthOfCoverage' and gatk_type == 'restricted':
            info(cnf['log'], 'Notice: in the restricted Gatk version, DepthOfCoverage '
                 'is renamed to Coverage. Using the name Coverage.\n')
            ann = 'Coverage'
        if ann == 'Coverage' and gatk_type == 'lite':
            info(cnf['log'], 'Notice: in the lite Gatk version, the Coverage annotation '
                 'goes by name of DepthOfCoverage. '
                 'In the system config, the lite version of Gatk is '
                 'specified; using DepthOfCoverage.\n')
            ann = 'DepthOfCoverage'
        cmdline += " -A " + ann

    output_fpath = intermediate_fname(work_dir, input_fpath, 'gatk')
    return _call_and_remove(cnf, cmdline, input_fpath, output_fpath,
                            stdout_to_outputfile=False,
                            to_remove=[output_fpath + '.idx',
                                       input_fpath + '.idx'])


def rename_fields(cnf, inp_tsv_fpath, field_map, work_dir):
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


def extract_fields(cnf, vcf_fpath, work_dir, sample_name=None):
    step_greetings(cnf, 'Extracting fields')

    name, _ = splitext_plus(vcf_fpath)
    tsv_fpath = name + '.tsv'

    if cnf.get('keep_intermediate') and cnf.get('reuse_intermediate'):
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
                                           'split_format_fields',
                                           keep_original_if_not_keep_intermediate=True)

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
    snpsift_cmline = _get_java_tool_cmdline(cnf, 'snpsift')

    if not which('perl'):
        exit('Perl executable required, maybe you need to run "module load perl"?')
    src_fpath = join(dirname(realpath(__file__)), 'src')
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


def filter_ensemble(cnf, input_fpath):
    step_greetings(cnf, 'Extracting dataset by filename, filtering ensemble reject line.')
    return iterate_file(cnf, input_fpath, lambda l: 'REJECT' not in l, 'pass')


def split_genotypes(cnf, input_fpath):
    step_greetings(cnf, 'Splitting genotypes.')

    def proc_line(line):
        if line.startswith('#'):
            return line
        else:
            tokens = line.split()
            alt_field = remove_quotes(tokens[4])
            alts = alt_field.split(',')
            if len(alts) > 1:
                for alt in set(alts):
                    line = '\t'.join(tokens[:2] + ['.'] + [tokens[3]] + [alt] + tokens[5:8])
                    return line
            else:
                line = '\t'.join(tokens[:2] + ['.'] + tokens[3:8])
                return line

    return iterate_file(cnf, input_fpath, proc_line, 'split_gt')


def filter_fields(cnf, input_fpath, work_dir):
    step_greetings(cnf, 'Filtering incorrect fields.')

    def proc_line(line):
        if not line.startswith('#'):
            if ',.' in line or '.,' in line:
                fields = line.split('\t')
                info_line = fields[7]
                info_pairs = [attr.split('=') for attr in info_line.split(';')]
                new_info_pairs = []
                for p in info_pairs:
                    if len(p) == 2:
                        if p[1].endswith(',.'):
                            p[1] = p[1][:-2]
                        if p[1].startswith('.,'):
                            p[1] = p[1][2:]
                        new_info_pairs.append('='.join(p))
                info_line = ';'.join(new_info_pairs)
                fields = fields[:7] + [info_line] + fields[8:]
                return '\t'.join(fields)
        return line

    return iterate_file(cnf, input_fpath, proc_line, work_dir)


def step_greetings(cnf, name):
    info(cnf.get('log'), '')
    info(cnf.get('log'), '-' * 70)
    info(cnf.get('log'), name)
    info(cnf.get('log'), '-' * 70)


def _check_system_resources(cnf):
    to_exit = False
    if not which('java'):
        err(cnf['log'], '\n* Warning: Java not found. You may want to run "module load java", '
            'or better ". /group/ngs/bin/bcbio-prod.sh"\n')
        to_exit = True

    if not which('perl'):
        err(cnf['log'], '\n* Warning: Perl not found. You may want to run "module load perl", '
            'or better ". /group/ngs/bin/bcbio-prod.sh"\n')
    if not _get_tool_cmdline(cnf, 'vcfannotate',
                             extra_warn='You may want to load BCBio '
                                        'with ". /group/ngs/bin/bcbio-prod.sh"'):
        err(cnf['log'], '\n* Warning: skipping annotation with bed tracks.\n')

    # print ''
    # print 'In Waltham, run this as well:'
    # print '   export PATH=$PATH:/group/ngs/src/snpEff/snpEff3.5/scripts'
    # print '   export PERL5LIB=$PERL5LIB:/opt/az/local/bcbio-nextgen/'
    #       'stable/0.7.6/tooldir/lib/perl5/site_perl'

    resources = cnf.get('resources', None)
    if not resources:
        critical(cnf['log'], 'No "resources" section in system config.')

    for name, data in resources.items():
        if 'path' in data:
            if not verify_file(data['path'], name):
                to_exit = True

    if to_exit:
        exit()


def _get_tool_cmdline(sys_cnf, tool_name, extra_warn=''):
    tool_path = which(tool_name) or None

    if not 'resources' in sys_cnf \
            or tool_name not in sys_cnf['resources'] \
            or 'path' not in sys_cnf['resources'][tool_name]:
        if tool_path:
            return tool_path
        else:
            err(tool_name + ' executable was not found. '
                'You can either specify path in the system config, or load into your '
                'PATH environment variable.')
            if extra_warn:
                err(extra_warn)
            return None

    tool_path = sys_cnf['resources'][tool_name]['path']
    if verify_file(tool_path, tool_name):
        return tool_path
    else:
        err(tool_path + ' for ' + tool_name + ' does not exist or is not a file.')
        return None


def join_parent_conf(child_conf, parent_conf):
    bc = parent_conf.copy()
    bc.update(child_conf)
    child_conf.update(bc)
    return child_conf


def _read_samples_info(cnf):
    if 'details' not in cnf:
        critical(cnf['log'], 'ERROR: Run config does not contain "details" section.')

    cnf['output_dir'] = cnf.get('output_dir', os.getcwd())
    cnf['filter_reject'] = cnf.get('filter_reject', False)
    cnf['split_samples'] = cnf.get('split_samples', False)

    all_samples = OrderedDict()

    details = cnf['details']
    del cnf['details']

    for vcf_conf in details:
        if 'vcf' not in vcf_conf:
            critical(cnf['log'], 'ERROR: A section in details does not contain field "vcf".')
        orig_vcf = vcf_conf['vcf']
        if not verify_file(orig_vcf, 'Input file'):
            exit(1)

        join_parent_conf(vcf_conf, cnf)

        # multiple samples
        if 'samples' in vcf_conf or vcf_conf.get('split_samples'):
            # check bams if exist
            if 'samples' in vcf_conf:
                for header_sample_name, sample_conf in vcf_conf['samples'].items():
                    join_parent_conf(sample_conf, vcf_conf)

                    bam = sample_conf.get('bam')
                    if bam and not verify_file(bam, 'Bam file'):
                        exit()

            rec_samples = vcf_conf.get('samples') or OrderedDict()
            vcf_header_samples = _read_sample_names_from_vcf(orig_vcf)

            # compare input sample names to vcf header
            if rec_samples:
                for input_sample_name, sample_conf in rec_samples.items():
                    if input_sample_name not in vcf_header_samples:
                        critical(cnf['log'], 'ERROR: sample ' + input_sample_name +
                                 ' is not in VCF header ' + vcf_header_samples + '\n'
                                 'Available samples: ' + ', '.join(vcf_header_samples))
            # split vcf
            for header_sample_name in vcf_header_samples:
                if header_sample_name not in rec_samples:
                    rec_samples[header_sample_name] = dict(vcf_conf)
                if header_sample_name in all_samples:
                    critical(cnf['log'], 'ERROR: duplicated sample name: '
                                         + header_sample_name)

                data = rec_samples[header_sample_name]
                all_samples[header_sample_name] = data

                output_dirpath = realpath(data['output_dir'])
                safe_mkdir(output_dirpath, 'output_dir for ' + header_sample_name)
                work_dirpath = join(output_dirpath, 'intermediate')
                safe_mkdir(work_dirpath, 'intermediate directory')
                data['work_dir'] = work_dirpath

                extracted_vcf = extract_sample(cnf, orig_vcf, header_sample_name, work_dirpath)
                rec_samples[header_sample_name]['vcf'] = extracted_vcf

        else:  # single sample
            if 'bam' in vcf_conf:
                if not verify_file(vcf_conf['bam']):
                    exit()

            sample_name = splitext(basename(orig_vcf))[0]

            output_dirpath = realpath(vcf_conf['output_dir'])
            safe_mkdir(output_dirpath, 'output_dir for ' + sample_name)
            work_dirpath = join(output_dirpath, 'intermediate')
            safe_mkdir(work_dirpath, 'intermediate directory')
            vcf_conf['work_dir'] = work_dirpath

            all_samples[sample_name] = vcf_conf

    for name, cnf in all_samples.items():
        cnf['name'] = name

    return all_samples


# def _copy_vcf_to_final_dir(self, inp_fpath):
#     fname = basename(inp_fpath)
#
#     # if self.output_dir != realpath(dirname(inp_fpath)):
#         # if self.run_cnf.get('keep_intermediate'):
#         # new_fname = fname
#         # else:
#         #     new_fname = self.make_intermediate(fname, 'anno')
#     new_fpath = join(self.work_dir, fname)
#     if os.path.exists(new_fpath):
#         os.remove(new_fpath)
#     shutil.copyfile(inp_fpath, new_fpath)
#     return new_fpath
#     # else:
#     #     if not self.run_cnf.get('keep_intermediate'):
#     #         new_fpath = join(self.output_dir, self.make_intermediate(fname, 'anno'))
#     #         shutil.copyfile(inp_fpath, new_fpath)
#     #         return new_fpath
#     #     else:
#     #         return inp_fpath


def _read_sample_names_from_vcf(vcf_fpath):
    basic_fields = next(l.strip()[1:].split() for l in open(vcf_fpath)
                        if l.strip().startswith('#CHROM'))
    return basic_fields[9:]


def _load_genome_resources(cnf):
    if 'genome' not in cnf:
        critical('"genome" is not specified in run config.')
    if 'genomes' not in cnf:
        critical('"genomes" section is not specified in system config.')
    genome_name = cnf['genome']
    if genome_name not in cnf['genomes']:
        critical(genome_name + ' is not in "genomes section" of system config.')

    genome_cnf = cnf['genomes'][genome_name].copy()

    to_exit = False

    if 'seq' not in genome_cnf:
        err('Please, provide path to the reference file (seq).')
        to_exit = True

    if not verify_file(genome_cnf['seq'], 'Reference seq'):
        to_exit = True

    for f in 'dbsnp', 'cosmic', 'dbsnfp', '1000genomes':
        if f in genome_cnf:
            if not verify_file(genome_cnf[f], f):
                to_exit = True
    if 'snpeff' in genome_cnf:
        if not verify_dir(genome_cnf['snpeff'], 'snpeff'):
            to_exit = True

    if to_exit:
        exit(1)

    cnf['genome'] = genome_cnf
    del cnf['genomes']
    genome_cnf['name'] = genome_name

    info('Loaded resources for ' + genome_cnf['name'])


def _check_quality_control_config(cnf):
    qc_cnf = cnf.get('quality_control')
    if not qc_cnf:
        return

    if 'databases' not in qc_cnf:
        qc_cnf['databases'] = ['dbsnp']
        info('Warning: not databases for quality control, using [dbsnp]')

    if 'novelty' not in qc_cnf:
        qc_cnf['novelty'] = ['all', 'known', 'novel']
        info('Warning: no novelty specified for quality control, '
             'using default ' + ', '.join(qc_cnf['novelty']))

    if 'metrics' not in qc_cnf:
        qc_cnf['metircs'] = [
           'nEvalVariants', 'nSNPs', 'nInsertions', 'nDeletions',
           'nVariantsAtComp', 'compRate', 'nConcordant', 'concordantRate',
           'variantRate', 'variantRatePerBp', 'hetHomRatio', 'tiTvRatio']
        info('Warning: no metrics for quality control, using '
             'default ' + ', '.join(qc_cnf['metircs']))

    to_exit = False
    dbs_dict = {}
    for db in qc_cnf['databases']:
        if not db:
            err('Empty field for quality_control databases')
            to_exit = True
        elif file_exists(db):
            if not verify_file(db):
                to_exit = True
            dbs_dict[basename(db)] = db
        elif db not in cnf['genome']:
            to_exit = True
            err(cnf['log'], db + ' for variant qc is not in ')
        else:
            dbs_dict[db] = cnf['genome'][db]

    if to_exit:
        exit()

    qc_cnf['database_vcfs'] = dbs_dict
    print str(dbs_dict)


def process_config(system_config_path, run_config_path):
    sys_cnf = load(open(system_config_path), Loader=Loader)
    run_cnf = load(open(run_config_path), Loader=Loader)

    info('Loaded system config ' + system_config_path)
    info('Loaded run config ' + run_config_path)
    config = dict(run_cnf.items() + sys_cnf.items())

    info('')
    info('Checking configs...')

    _check_system_resources(config)

    _load_genome_resources(config)

    if 'quality_control' in config:
        _check_quality_control_config(config)

    samples = _read_samples_info(config)

    info('Configs checked.')
    return config, samples


def main(args):
    if len(args) < 1:
        exit('Usage: python ' + __file__ + ' run_info.yaml\n'
             '    or python ' + __file__ + ' system_info.yaml run_info.yaml')

    if (2, 5) > sys.version_info[:2] >= (3, 0):
        exit('Python version 2.5 and higher is supported (you are running ' +
             '.'.join(map(str, sys.version_info[:3])) + ')\n')

    if len(args) == 1:
        run_config_path = args[0]
        system_config_path = join(dirname(realpath(__file__)),
                                  'system_info_rask.yaml')
        sys.stderr.write('Notice: using system_info_rask.yaml as a default'
                         ' tools configutation file.\n\n')
    else:
        system_config_path = args[0]
        run_config_path = args[1]

    if not os.path.isfile(system_config_path):
        exit(system_config_path + ' does not exist or is a directory.\n')
    if not os.path.isfile(run_config_path):
        exit(run_config_path + ' does not exist or is a directory.\n')

    to_exit = False
    if not system_config_path.endswith('.yaml'):
        sys.stderr.write(system_config_path + ' does not end with .yaml, '
                                              'maybe incorrect parameter?\n')
        to_exit = True
    if not run_config_path.endswith('.yaml'):
        sys.stderr.write(run_config_path + ' does not end with .yaml,'
                                           ' maybe incorrect parameter?\n')
        to_exit = True
    if to_exit:
        exit()

    config, samples = process_config(system_config_path, run_config_path)
    try:
        annotate(samples, config.get('parallel'))
    except KeyboardInterrupt:
        exit()


if __name__ == '__main__':
    main(sys.argv[1:])