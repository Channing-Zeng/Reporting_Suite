# Annotation script that takes 1-3 inputs, first being the vcf file name,
# second being an indicator if the vcf is from bcbio's ensemble pipeline ('true' if true) and
# third being 'RNA' if the vcf is from the rna-seq mutect pipeline
from genericpath import isfile, getsize
import os
from os.path import join, splitext
import subprocess
import sys
import shutil
from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def log_print(msg=''):
    print(msg)
    if 'log' in run_config:
        open(run_config['log'], 'a').write(msg + '\n')


def error(msg):
    sys.stderr.write(msg + '\n')
    exit(1)


def check_executable(program):
    if not which(program):
        error(program + ' executable required.')


def check_existence(file):
    if not isfile(file) or getsize(file) <= 0:
        error(file + ' does not exist or is empty.')


run_config = {}
system_config = {}


def _call_and_rename(cmdline, input_fpath, suffix, to_stdout=True):
    basepath, ext = splitext(input_fpath)
    output_fpath = basepath + '.' + suffix + ext

    if run_config.get('reuse') and isfile(output_fpath) and getsize(output_fpath) > 0:
        log_print(output_fpath + ' exists, reusing')
        return output_fpath

    log_print('')
    log_print('*' * 70)
    log_print(cmdline)

    res = subprocess.call(
        cmdline.split(),
        stdout=open(output_fpath, 'w') if to_stdout else open(run_config['log'], 'a'),
        stderr=open(run_config['log'] + '_err', 'w'))
    if res != 0:
        with open(run_config['log'] + '_err') as err:
            log_print('')
            log_print(err.read())
            log_print('')
        log_print('Command returned status ' + str(res) + ('. Log in ' + run_config['log']))
        exit(1)
    else:
        with open(run_config['log'] + '_err') as err, open(run_config['log'], 'a') as log:
            log.write('')
            log.write(err.read())
            log.write('')
        log_print('Saved to ' + output_fpath)
        print('Log in ' + run_config['log'])

    if not run_config.get('save_intermediate'):
        os.remove(input_fpath)
    log_print('Now processing ' + output_fpath)
    return output_fpath


def _get_java_tool_cmdline(name):
    cmdline_pattern = _get_tool_cmdline('java', name)
    jvm_opts = system_config['resources'][name].get('jvm_opts', [])
    return cmdline_pattern % (' '.join(jvm_opts) + ' -jar')


def _get_tool_cmdline(executable, name):
    check_executable(executable)
    if 'resources' not in system_config:
        error('System config yaml must contain resources section with ' + name + ' path.')
    if name not in system_config['resources']:
        error('System config resources section must contain ' + name + ' info (with a path to the tool).')
    tool_config = system_config['resources'][name]
    if 'path' not in tool_config:
        error(name + ' section in the system config must contain a path to the tool.')
    tool_path = tool_config['path']
    check_existence(tool_path)
    return executable + ' %s ' + tool_path


def snpsift_annotate(db_name, input_fpath):
    db_path = run_config['vcfs'][db_name].get('path')
    annotations = run_config['vcfs'][db_name].get('annotations')

    cmdline = _get_java_tool_cmdline('snpsift') + ' annotate -v %s %s' % (db_path, input_fpath)
    return _call_and_rename(cmdline, input_fpath, db_name, to_stdout=True)


def snpsift_db_nsfp(input_fpath):
    if 'db_nsfp' not in run_config:
        return input_fpath

    db_path = run_config['db_nsfp'].get('path')
    assert db_path, 'Please, provide a path to db nsfp file in run_config.'

    annotations = run_config['db_nsfp'].get('annotations', [])
    ann_line = ','.join(annotations)

    cmdline = _get_java_tool_cmdline('snpsift') + ' dbnsfp -f %s -v %s %s' % (ann_line, db_path, input_fpath)
    return _call_and_rename(cmdline, input_fpath, 'db_nsfp', to_stdout=True)


def snpeff(input_fpath):
    if 'snpeff' not in run_config:
        return input_fpath

    ref_name = run_config['genome_build']

    db_path = run_config['snpeff'].get('path')
    assert db_path, 'Please, provide a path to db nsfp file in run_config.'

    cmdline = _get_java_tool_cmdline('snpeff') + ' eff -dataDir %s -noStats -cancer ' \
              '-noLog -1 -i vcf -o vcf %s %s' % (db_path, ref_name, input_fpath)
    return _call_and_rename(cmdline, input_fpath, 'snpEff', to_stdout=True)


def rna_editing_sites(db, input_fpath):
    check_executable('vcfannotate')

    cmdline = 'vcfannotate -b %s -k RNA_editing_site %s' % (db, input_fpath)
    return _call_and_rename(cmdline, input_fpath, 'edit', to_stdout=True)


def gatk(input_fpath):
    if 'gatk' not in run_config:
        return input_fpath

    base_name, ext = os.path.splitext(input_fpath)
    output_fpath = base_name + '.gatk' + ext

    ref_fpath = run_config['reference']

    cmdline = _get_java_tool_cmdline('gatk') + ' -R %s -T VariantAnnotator -o %s --variant %s' % (ref_fpath, output_fpath, input_fpath)

    annotations = run_config['db_nsfp'].get('annotations', [])
    for ann in annotations:
        cmdline += " -A " + ann

    return _call_and_rename(cmdline, input_fpath, 'gatk', to_stdout=False)


def extract_fields(input_fpath):
    snpeff_cmline = _get_java_tool_cmdline('snpeff')
    vcfoneperline_cmline = _get_tool_cmdline('perl', 'vcfEffOnePerLine.pl') % ''

    cmdline = vcfoneperline_cmline + ' | ' + \
              snpeff_cmline + ' extractFields - ' \
              'CHROM POS ID CNT GMAF REF ALT QUAL FILTER TYPE ' \
              '"EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].CODON" ' \
              '"EFF[*].AA" "EFF[*].AA_LEN" "EFF[*].GENE" ' \
              '"EFF[*].FUNCLASS" "EFF[*].BIOTYPE" "EFF[*].CODING" ' \
              '"EFF[*].TRID" "EFF[*].RANK" ' \
              'dbNSFP_SIFT_score dbNSFP_Polyphen2_HVAR_score ' \
              'dbNSFP_Polyphen2_HVAR_pred dbNSFP_LRT_score dbNSFP_LRT_pred ' \
              'dbNSFP_MutationTaster_score dbNSFP_MutationTaster_pred ' \
              'dbNSFP_MutationAssessor_score dbNSFP_MutationAssessor_pred ' \
              'dbNSFP_FATHMM_score dbNSFP_Ensembl_geneid dbNSFP_Ensembl_transcriptid ' \
              'dbNSFP_Uniprot_acc dbNSFP_1000Gp1_AC dbNSFP_1000Gp1_AF ' \
              'dbNSFP_ESP6500_AA_AF dbNSFP_ESP6500_EA_AF KGPROD PM PH3 ' \
              'AB AC AF DP FS GC HRun HaplotypeScore ' \
              'G5 CDA GMAF GENEINFO OM DB GENE AA CDS ' \
              'MQ0 QA QD ReadPosRankSum '

    basepath, ext = os.path.splitext(input_fpath)
    output_fpath = basepath + '.extract' + ext

    if run_config.get('reuse') and isfile(output_fpath) and getsize(output_fpath):
        log_print(output_fpath + ' exists, reusing')
    else:
        log_print('')
        log_print('*' * 70)
        log_print(cmdline)
        res = subprocess.call(cmdline,
                              stdin=open(sample_fpath),
                              stdout=open(output_fpath, 'w'),
                              stderr=open(run_config['log'], 'a'),
                              shell=True)
        log_print('')
        if res != 0:
            log_print('Command returned status ' + str(res) + ('. Log in ' + run_config['log']))
            exit(1)
            # return input_fpath
        else:
            log_print('Saved to ' + output_fpath)
            print('Log in ' + run_config['log'])

        if not run_config.get('save_intermediate'):
            os.remove(sample_fpath)

    os.rename(sample_fpath, splitext(sample_fpath)[0] + '.tsv')


def process_rna(sample_fpath):
    sample_fname = os.path.basename(sample_fpath)
    sample_basename, ext = os.path.splitext(sample_fname)
    check_executable('vcf-subset')
    cmdline = 'vcf-subset -c %s -e %s' % (sample_basename.replace('-ensemble', ''), sample_fpath)

    return _call_and_rename(cmdline, sample_fpath, '.ensm', to_stdout=True)


def process_ensemble(sample_fpath):
    sample_basepath, ext = os.path.splitext(sample_fpath)
    pass_sample_fpath = sample_basepath + '.pass' + ext
    with open(sample_fpath) as sample, open(pass_sample_fpath, 'w') as pass_sample:
        for line in sample.readlines():
            if 'REJECT' not in line:
                pass_sample.write(line)
    if run_config.get('save_intermediate'):
        return pass_sample_fpath
    else:
        os.remove(sample_fpath)
        os.rename(pass_sample_fpath, sample_fpath)
        return sample_fpath


def annotate(sample_fpath):
    assert 'resources' in system_config

    if run_config.get('rna'):
        sample_fpath = process_rna(sample_fpath)

    if run_config.get('ensemble'):
        sample_fpath = process_ensemble(sample_fpath)

    assert 'genome_build' in run_config, 'Please, provide genome build (genome_build).'
    assert 'reference' in run_config, 'Please, provide path to the reference file (reference).'
    check_existence(run_config['reference'])

    if 'vcfs' in run_config:
        for vcf in run_config['vcfs']:
            sample_fpath = snpsift_annotate(vcf, sample_fpath)

    sample_fpath = snpsift_db_nsfp(sample_fpath)

    #if run_config.get('rna'):
    #    sample_fpath = rna_editing_sites(annot_track, sample_fpath, save_intermediate)

    sample_fpath = gatk(sample_fpath)
    sample_fpath = snpeff(sample_fpath)
    extract_fields(sample_fpath)


def remove_quotes(str):
    if str and str[0] == '"':
        str = str[1:]
    if str and str[-1] == '"':
        str = str[:-1]
    return str


def split_genotypes(sample_fpath, result_fpath):
    with open(sample_fpath) as vcf, open(result_fpath, 'w') as out:
        for i, line in enumerate(vcf):
            clean_line = line.strip()
            if not clean_line or clean_line[0] == '#':
                out.write(line)
            else:
                tokens = line.split()
                alt_field = remove_quotes(tokens[4])
                alts = alt_field.split(',')
                if len(alts) > 1:
                    for alt in alts:
                        line = '\t'.join(tokens[:2] + ['.'] + [tokens[3]] + [alt] + tokens[5:]) + '\n'
                        out.write(line)
                else:
                    line = '\t'.join(tokens[:2] + ['.'] + tokens[3:]) + '\n'
                    out.write(line)

    if run_config.get('save_intermediate'):
        return result_fpath
    else:
        os.remove(sample_fpath)
        os.rename(result_fpath, sample_fpath)
        return sample_fpath


if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) < 2:
        sys.stderr.write('Usage: python ' + __file__ + ' system_info_local.yaml run_info.yaml\n')
        exit(1)

    assert os.path.isfile(args[0]), args[0] + ' does not exist of is a directory.'
    assert os.path.isfile(args[1]), args[1] + ' does not exist of is a directory.'
    system_config = load(open(args[0]), Loader=Loader)
    run_config = load(open(args[1]), Loader=Loader)
    print('Loaded system config ' + args[0])
    print('Loaded run config ' + args[1])

    sample_fpath = os.path.realpath(run_config.get('file', None))
    assert os.path.isfile(sample_fpath), sample_fpath + ' does not exists or is not a file.'

    result_dir = os.path.realpath(run_config.get('output_dir', os.getcwd()))
    assert os.path.isdir(result_dir), result_dir + ' does not exists or is not a directory'

    sample_fname = os.path.basename(sample_fpath)
    sample_basename, ext = os.path.splitext(sample_fname)

    if result_dir != os.path.realpath(os.path.dirname(sample_fpath)):
        new_sample_fpath = os.path.join(result_dir, sample_fname)
        if os.path.exists(new_sample_fpath):
            os.remove(new_sample_fpath)
        shutil.copyfile(sample_fpath, new_sample_fpath)
        sample_fpath = new_sample_fpath

    if 'log' not in run_config:
        run_config['log'] = os.path.join(os.path.dirname(sample_fpath), sample_basename + '.log')
    if os.path.isfile(run_config['log']):
        os.remove(run_config['log'])

    log_print('Writing into ' + result_dir)
    log_print('Logging to ' + run_config['log'])
    log_print('')

    print('Note: please, load modules before start:')
    print('   source /etc/profile.d/modules.sh')
    print('   module load java')
    print('   module load perl')
    # print ''
    # print 'In Waltham, run this as well:'
    # print '   export PATH=$PATH:/group/ngs/src/snpEff/snpEff3.5/scripts'
    # print '   export PERL5LIB=$PERL5LIB:/opt/az/local/bcbio-nextgen/stable/0.7.6/tooldir/lib/perl5/site_perl'

    if run_config.get('split_genotypes'):
        sample_basepath, ext = os.path.splitext(sample_fpath)
        result_fpath = sample_basepath + '.split' + ext
        log_print('')
        log_print('*' * 70)
        log_print('Splitting genotypes.')
        sample_fpath = split_genotypes(sample_fpath, result_fpath)
        log_print('Saved to ' + result_fpath)
        log_print('')

    annotate(sample_fpath)