from genericpath import isfile
from os.path import splitext, basename, join, dirname, realpath
import os
import shutil
from source.logger import step_greetings

from source.utils import critical, iterate_file, \
    get_java_tool_cmdline, verify_file, intermediate_fname, call_subprocess, \
    get_tool_cmdline, err, get_gatk_type, info, index_bam
from source.utils_from_bcbio import add_suffix, file_exists


def run_annotators(cnf, vcf_fpath, bam_fpath=None):
    annotated = False
    original_vcf = vcf_fpath

    if 'gatk' in cnf:
        annotated = True
        vcf_fpath = _gatk(cnf, vcf_fpath, bam_fpath)

    if 'dbsnp' in cnf:
        annotated = True
        vcf_fpath = _snpsift_annotate(cnf, cnf['dbsnp'], 'dbsnp', vcf_fpath)
    if 'cosmic' in cnf:
        annotated = True
        vcf_fpath = _snpsift_annotate(cnf, cnf['cosmic'], 'cosmic', vcf_fpath)
    if 'custom_vcfs' in cnf:
        for dbname, custom_conf in cnf['custom_vcfs'].items():
            annotated = True
            vcf_fpath = _snpsift_annotate(
                cnf, custom_conf, dbname, vcf_fpath)

    if 'dbnsfp' in cnf:
        annotated = True
        vcf_fpath = _snpsift_db_nsfp(cnf, vcf_fpath)

    if 'snpeff' in cnf:
        annotated = True
        _remove_annotation(cnf, 'EFF', vcf_fpath)
        vcf_fpath, summary_fpath, genes_fpath = _snpeff(cnf, vcf_fpath)
        if isfile(join(cnf['output_dir'], summary_fpath)):
            os.remove(join(cnf['output_dir'], summary_fpath))
        if isfile(join(cnf['output_dir'], genes_fpath)):
            os.remove(join(cnf['output_dir'], genes_fpath))
        if file_exists(summary_fpath):
            shutil.move(summary_fpath, cnf['output_dir'])
        if file_exists(genes_fpath):
            shutil.move(genes_fpath, cnf['output_dir'])

    next_vcf_fpath = None
    if cnf.get('tracks'):
        for track in cnf['tracks']:
            next_vcf_fpath = _tracks(cnf, track, vcf_fpath)
            if next_vcf_fpath:
                annotated = True
                vcf_fpath = next_vcf_fpath

    if annotated:
        vcf_fpath = _filter_fields(cnf, vcf_fpath)

        # Copying final VCF
        final_vcf_fname = add_suffix(basename(cnf['vcf']), 'anno')
        final_vcf_fpath = join(cnf['output_dir'], final_vcf_fname)
        if isfile(final_vcf_fpath):
            os.remove(final_vcf_fpath)
        shutil.copyfile(vcf_fpath, final_vcf_fpath)

        # Converting to MAF
        final_maf_fpath = _convert_to_maf(cnf, final_vcf_fpath)

        return final_vcf_fpath, final_maf_fpath
    else:
        info('No annotations were run on ' + original_vcf + '. Please, specify some in run_info.')
        return None, None


def _remove_annotation(cnf, field_to_del, input_fpath):
    def proc_line(l):
        if field_to_del in l:
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
                info_line = ';'.join('='.join(pair) if len(pair) == 2
                                     else pair[0] for pair in info_pairs)
                fields = fields[:7] + [info_line] + fields[8:]
                return '\t'.join(fields)
        return l
    return iterate_file(cnf, input_fpath, proc_line)


def _convert_to_maf(cnf, final_vcf_fpath):
    step_greetings('Converting to MAF')

    final_vcf_fname = basename(final_vcf_fpath)

    final_maf_fpath = join(cnf['output_dir'], splitext(final_vcf_fname)[0] + '.maf')
    perl = get_tool_cmdline(cnf, 'perl')
    vcf2maf = join(dirname(realpath(__file__)), '../../external/vcf2maf-1.1.0/vcf2maf.pl')
    cmdline = '{perl} {vcf2maf} --input-snpeff {final_vcf_fpath} ' \
              '--output-maf {final_maf_fpath}'.format(**locals())
    call_subprocess(cnf, cmdline, None, None)
    return final_maf_fpath


def _snpsift_annotate(cnf, vcf_conf, dbname, input_fpath):
    step_greetings('Annotate with ' + dbname)

    executable = get_java_tool_cmdline(cnf, 'snpsift')

    db_path = cnf['genome'].get(dbname)
    if not db_path:
        db_path = vcf_conf.get('path')
        if not db_path:
            critical('Please, privide a path to ' + dbname + ' in the run config '
                     '("path:" field), or in the "genomes" section in the system config')
        if not verify_file(db_path):
            exit()

    anno_line = ''
    annotations = vcf_conf.get('annotations')
    if annotations:
        anno_line = '-info ' + ','.join(annotations)

    cmdline = '{executable} annotate -v {anno_line} {db_path} {input_fpath}'.format(**locals())
    output_fpath = intermediate_fname(cnf, input_fpath, dbname)
    output_fpath = call_subprocess(cnf, cmdline, input_fpath, output_fpath,
                        stdout_to_outputfile=True)

    # all_fields.extend(annotations)

    def proc_line(line):
        if not line.startswith('#'):
            line = line.replace(' ', '_')
            assert ' ' not in line
        return line
    output_fpath = iterate_file(cnf, output_fpath, proc_line)
    return output_fpath


def _snpsift_db_nsfp(cnf, input_fpath):
    if 'dbnsfp' not in cnf:
        return input_fpath

    step_greetings('DB SNFP')

    executable = get_java_tool_cmdline(cnf, 'snpsift')

    db_path = cnf['genome'].get('dbnsfp')
    if not db_path:
        critical('Please, provide a path to DB NSFP file in '
                 'the "genomes" section in the system config.')

    annotations = cnf['dbnsfp'].get('annotations', [])

    # all_fields.extend(['dbNSFP_' + ann for ann in annotations])

    ann_line = ('-f ' + ','.join(annotations)) if annotations else ''

    cmdline = '{executable} dbnsfp {ann_line} -v {db_path} ' \
              '{input_fpath}'.format(**locals())
    output_fpath = intermediate_fname(cnf, input_fpath, 'db_nsfp')
    return call_subprocess(cnf, cmdline, input_fpath, output_fpath,
                stdout_to_outputfile=True)


def _snpeff(cnf, input_fpath):
    if 'snpeff' not in cnf:
        return input_fpath, None, None

    step_greetings('SnpEff')

    # self.all_fields.extend([
    #     "EFF[*].EFFECT", "EFF[*].IMPACT", "EFF[*].FUNCLASS", "EFF[*].CODON",
    #     "EFF[*].AA", "EFF[*].AA_LEN", "EFF[*].GENE", "EFF[*].CODING",
    #     "EFF[*].TRID", "EFF[*].RANK"])

    executable = get_java_tool_cmdline(cnf, 'snpeff')
    ref_name = cnf['genome']['name']
    db_path = cnf['genome'].get('snpeff')
    if not db_path:
        critical('Please, provide a path to SnpEff data in '
                 'the "genomes" section in the system config.')

    stats_fpath = join(cnf['name'] + '.snpEff_summary.html')
    cmdline = ('{executable} eff -dataDir {db_path} -stats {stats_fpath} '
               '-csvStats -noLog -1 -i vcf -o vcf {ref_name} '
               '{input_fpath}').format(**locals())

    if cnf['snpeff'].get('clinical_reporting') or \
            cnf['snpeff'].get('canonical'):
        cmdline += ' -canon -hgvs '

    if cnf['snpeff'].get('cancer'):
        cmdline += ' -cancer '

    output_fpath = intermediate_fname(cnf, input_fpath, 'snpEff')
    return call_subprocess(cnf, cmdline, input_fpath, output_fpath,
                stdout_to_outputfile=True), \
        stats_fpath, splitext(stats_fpath)[0] + '.genes.txt'


def _tracks(cnf, track_path, input_fpath):
    field_name = splitext(basename(track_path))[0]

    step_greetings('Intersecting with ' + field_name)

    toolpath = get_tool_cmdline(cnf, 'vcfannotate')
    if not toolpath:
        err('WARNING: Skipping annotation with tracks: vcfannotate '
            'executable not found, you probably need to specify path in system_config, or '
            'run load bcbio:  . /group/ngs/bin/bcbio-prod.sh"')
        return

    # self.all_fields.append(field_name)

    cmdline = '{toolpath} -b {track_path} -k {field_name} {input_fpath}'.format(**locals())

    output_fpath = intermediate_fname(cnf, input_fpath, field_name)
    output_fpath = call_subprocess(cnf, cmdline, input_fpath, output_fpath,
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
    return iterate_file(cnf, output_fpath, proc_line)


def _gatk(cnf, input_fpath, bam_fpath):
    if 'gatk' not in cnf:
        return input_fpath

    step_greetings('GATK')

    if bam_fpath:
        index_bam(cnf, bam_fpath)

    executable = get_java_tool_cmdline(cnf, 'gatk')

    gatk_opts_line = ''
    if cnf.gatk:
        if 'options' in cnf.gatk:
            gatk_opts_line = ' '.join(cnf.gatk['options'])

    if 'threads' in cnf and ' -nt ' not in gatk_opts_line:
        gatk_opts_line += ' -nt ' + str(cnf.threads)

    output_fpath = intermediate_fname(cnf, input_fpath, 'gatk')

    # duplicating this from "call" function to avoid calling "gatk --version"
    if output_fpath and cnf.get('reuse_intermediate'):
        if file_exists(output_fpath):
            info(output_fpath + ' exists, reusing')
            return output_fpath

    # Avoid issues with incorrectly created empty GATK index files.
    # Occurs when GATK cannot lock shared dbSNP database on previous run
    # (not using dbSNP here anymore, but removing old inx just in case)
    idx_fpath = input_fpath + '.idx'
    if os.path.exists(idx_fpath) and not file_exists(idx_fpath):
        os.remove(idx_fpath)

    ref_fpath = cnf['genome']['seq']

    cmdline = ('{executable} {gatk_opts_line} -R {ref_fpath} -T VariantAnnotator'
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

    gatk_type = get_gatk_type(get_java_tool_cmdline(cnf, 'gatk'))
    for ann in annotations:
        if ann == 'DepthOfCoverage' and gatk_type == 'restricted':
            info('Notice: in the restricted Gatk version, DepthOfCoverage '
                 'is renamed to Coverage. Using the name Coverage.')
            info()
            ann = 'Coverage'
        if ann == 'Coverage' and gatk_type == 'lite':
            info('Notice: in the lite Gatk version, the Coverage annotation '
                 'goes by name of DepthOfCoverage. '
                 'In the system config, the lite version of Gatk is '
                 'specified; using DepthOfCoverage.')
            info()
            ann = 'DepthOfCoverage'
        cmdline += " -A " + ann

    output_fpath = intermediate_fname(cnf, input_fpath, 'gatk')
    return call_subprocess(cnf, cmdline, input_fpath, output_fpath,
                stdout_to_outputfile=False,
                to_remove=[output_fpath + '.idx',
                           input_fpath + '.idx'])


def _filter_fields(cnf, input_fpath):
    step_greetings('Filtering incorrect fields.')

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

    output_fpath = iterate_file(cnf, input_fpath, proc_line)

    info('Saved to ' + output_fpath)
    return output_fpath