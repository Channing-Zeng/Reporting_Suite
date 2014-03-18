# Annotation script that takes 1-3 inputs, first being the vcf file name,
# second being an indicator if the vcf is from bcbio's ensemble pipeline ('true' if true) and
# third being 'RNA' if the vcf is from the rna-seq mutect pipeline
from genericpath import isfile, getsize
import os
from os.path import join, splitext, basename, realpath
import subprocess
import sys
import shutil
from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


def which(program):
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


def error(msg):
    sys.stderr.write(msg + '\n')
    exit(1)


def check_executable(program):
    if not which(program):
        error(program + ' executable required.')


def check_existence(file):
    if not file or not isfile(file) or getsize(file) <= 0:
        error(file + ' does not exist, is not a file, or is empty.')


class Annotator:
    def __init__(self, system_config, run_config):
        self.system_config = system_config
        self.run_config = run_config

        sample_fpath = self.run_config.get('file', None)
        assert sample_fpath, 'Run config does not contain field "file".'
        self.sample_fpath = realpath(sample_fpath)
        check_existence(self.sample_fpath)

        result_dir = realpath(run_config.get('output_dir', os.getcwd()))
        assert os.path.isdir(result_dir), result_dir + ' does not exists or is not a directory'

        sample_fname = os.path.basename(self.sample_fpath)
        sample_basename, ext = os.path.splitext(sample_fname)

        if result_dir != os.path.realpath(os.path.dirname(self.sample_fpath)):
            new_sample_fpath = os.path.join(result_dir, sample_fname)
            if os.path.exists(new_sample_fpath):
                os.remove(new_sample_fpath)
            shutil.copyfile(self.sample_fpath, new_sample_fpath)
            self.sample_fpath = new_sample_fpath

        if 'log' not in run_config:
            run_config['log'] = os.path.join(os.path.dirname(self.sample_fpath), sample_basename + '.log')
        if os.path.isfile(run_config['log']):
            os.remove(run_config['log'])

        self.log_print('Writing into ' + result_dir)
        self.log_print('Logging to ' + run_config['log'])
        self.log_print('')


    def log_print(self, msg=''):
        print(msg)
        if 'log' in self.run_config:
            open(self.run_config['log'], 'a').write(msg + '\n')


    def _call_and_rename(self, cmdline, input_fpath, suffix, to_stdout=True):
        basepath, ext = splitext(input_fpath)
        output_fpath = basepath + '.' + suffix + ext

        if self.run_config.get('reuse') and isfile(output_fpath) and getsize(output_fpath) > 0:
            self.log_print(output_fpath + ' exists, reusing')
            return output_fpath

        self.log_print(cmdline)

        res = subprocess.call(
            cmdline.split(),
            stdout=open(output_fpath, 'w') if to_stdout else open(self.run_config['log'], 'a'),
            stderr=open(self.run_config['log'] + '_err', 'w'))
        if res != 0:
            with open(self.run_config['log'] + '_err') as err:
                self.log_print('')
                self.log_print(err.read())
                self.log_print('')
            self.log_print('Command returned status ' + str(res) + ('. Log in ' + self.run_config['log']))
            exit(1)
        else:
            with open(self.run_config['log'] + '_err') as err, open(self.run_config['log'], 'a') as log:
                log.write('')
                log.write(err.read())
                log.write('')
            self.log_print('Saved to ' + output_fpath)
            print('Log in ' + self.run_config['log'])

        if not self.run_config.get('save_intermediate'):
            os.remove(input_fpath)
        self.log_print('Now processing ' + output_fpath)
        return output_fpath


    def _get_java_tool_cmdline(self, name):
        cmdline_pattern = self._get_tool_cmdline('java', name)
        jvm_opts = self.system_config['resources'][name].get('jvm_opts', [])
        return cmdline_pattern % (' '.join(jvm_opts) + ' -jar')


    def _get_tool_cmdline(self, executable, name):
        check_executable(executable)
        if 'resources' not in self.system_config:
            error('System config yaml must contain resources section with ' + name + ' path.')
        if name not in self.system_config['resources']:
            error('System config resources section must contain ' + name + ' info (with a path to the tool).')
        tool_config = self.system_config['resources'][name]
        if 'path' not in tool_config:
            error(name + ' section in the system config must contain a path to the tool.')
        tool_path = tool_config['path']
        check_existence(tool_path)
        return executable + ' %s ' + tool_path


    def snpsift_annotate(self, dbname, conf, input_fpath):
        self.log_print('')
        self.log_print('*' * 70)

        db_path = conf.get('path')
        annotations = conf.get('annotations', [])
        anno_line = ','.join(annotations)

        cmdline = self._get_java_tool_cmdline('snpsift') + \
                  ' annotate -v -info %s %s %s' % \
                  (anno_line, db_path, input_fpath)
        return self._call_and_rename(cmdline, input_fpath, dbname, to_stdout=True)


    def snpsift_db_nsfp(self, input_fpath):
        if 'db_nsfp' not in self.run_config:
            return input_fpath

        self.log_print('')
        self.log_print('*' * 70)

        db_path = self.run_config['db_nsfp'].get('path')
        assert db_path, 'Please, provide a path to db nsfp file in run_config.'

        annotations = self.run_config['db_nsfp'].get('annotations', [])
        ann_line = ','.join(annotations)

        cmdline = self._get_java_tool_cmdline('snpsift') + ' dbnsfp -f %s -v %s %s' % (ann_line, db_path, input_fpath)
        return self._call_and_rename(cmdline, input_fpath, 'db_nsfp', to_stdout=True)


    def snpeff(self, input_fpath):
        if 'snpeff' not in self.run_config:
            return input_fpath

        self.log_print('')
        self.log_print('*' * 70)

        ref_name = self.run_config['genome_build']

        db_path = self.run_config['snpeff'].get('path')
        assert db_path, 'Please, provide a path to db nsfp file in run_config.'

        cmdline = self._get_java_tool_cmdline('snpeff') + ' eff -dataDir %s -noStats -cancer ' \
                  '-noLog -1 -i vcf -o vcf %s %s' % (db_path, ref_name, input_fpath)
        return self._call_and_rename(cmdline, input_fpath, 'snpEff', to_stdout=True)


    def tracks(self, track_path, input_fpath):
        self.log_print('')
        self.log_print('*' * 70)

        check_executable('vcfannotate')

        field_name = splitext(basename(track_path))[0]
        cmdline = 'vcfannotate -b %s -k %s %s' % (track_path, field_name, input_fpath)
        return self._call_and_rename(cmdline, input_fpath, field_name, to_stdout=True)


    def gatk(self, input_fpath):
        if 'gatk' not in self.run_config:
            return input_fpath

        self.log_print('')
        self.log_print('*' * 70)

        base_name, ext = os.path.splitext(input_fpath)
        output_fpath = base_name + '.gatk' + ext

        ref_fpath = self.run_config['reference']

        cmdline = self._get_java_tool_cmdline('gatk') + ' -R %s -T VariantAnnotator -o %s --variant %s' % (ref_fpath, output_fpath, input_fpath)

        annotations = self.run_config['gatk'].get('annotations', [])
        for ann in annotations:
            cmdline += " -A " + ann

        return self._call_and_rename(cmdline, input_fpath, 'gatk', to_stdout=False)


    def extract_fields(self, input_fpath):
        if 'tsv_fields' not in self.run_config:
            return

        self.log_print('')
        self.log_print('*' * 70)

        anno_line = ' '.join(self.run_config['tsv_fields'])

        snpsift_cmline = self._get_java_tool_cmdline('snpsift')
        vcfoneperline_cmline = self._get_tool_cmdline('perl', 'vcfoneperline') % ''

        cmdline = vcfoneperline_cmline + ' | ' + snpsift_cmline + ' extractFields - ' + anno_line

        basepath, ext = os.path.splitext(input_fpath)
        tsv_fpath = basepath + '.extract.tsv'
        if isfile(tsv_fpath):
            os.remove(tsv_fpath)

        self.log_print(cmdline)
        res = subprocess.call(cmdline,
                              stdin=open(input_fpath),
                              stdout=open(tsv_fpath, 'w'),
                              stderr=open(self.run_config['log'], 'a'),
                              shell=True)
        self.log_print('')
        if res != 0:
            self.log_print('Command returned status ' + str(res) + ('. Log in ' + self.run_config['log']))
            exit(1)
            # return input_fpath
        else:
            self.log_print('Saved TSV file to ' + tsv_fpath)
            print('Log in ' + self.run_config['log'])


    def process_rna(self, sample_fpath):
        self.log_print('')
        self.log_print('*' * 70)

        sample_fname = os.path.basename(sample_fpath)
        sample_basename, ext = os.path.splitext(sample_fname)
        check_executable('vcf-subset')
        cmdline = 'vcf-subset -c %s -e %s' % (sample_basename.replace('-ensemble', ''), sample_fpath)

        return self._call_and_rename(cmdline, sample_fpath, '.ensm', to_stdout=True)


    def process_ensemble(self, sample_fpath):
        self.log_print('')
        self.log_print('*' * 70)
        self.log_print('Filtering ensemble reject lines.')

        sample_basepath, ext = os.path.splitext(sample_fpath)
        pass_sample_fpath = sample_basepath + '.pass' + ext
        with open(sample_fpath) as sample, open(pass_sample_fpath, 'w') as pass_sample:
            for line in sample.readlines():
                if 'REJECT' not in line:
                    pass_sample.write(line)
        if self.run_config.get('save_intermediate'):
            return pass_sample_fpath
        else:
            os.remove(sample_fpath)
            os.rename(pass_sample_fpath, sample_fpath)
            return sample_fpath


    def annotate(self):
        assert 'resources' in self.system_config

        sample_fpath = self.sample_fpath

        if self.run_config.get('split_genotypes'):
            sample_basepath, ext = os.path.splitext(sample_fpath)
            result_fpath = sample_basepath + '.split' + ext
            sample_fpath = self.split_genotypes(sample_fpath, result_fpath)
            self.log_print('Saved to ' + result_fpath)

        if self.run_config.get('rna'):
            sample_fpath = self.process_rna(sample_fpath)

        if self.run_config.get('ensemble'):
            sample_fpath = self.process_ensemble(sample_fpath)

        assert 'genome_build' in self.run_config, 'Please, provide genome build (genome_build).'
        assert 'reference' in self.run_config, 'Please, provide path to the reference file (reference).'
        check_existence(self.run_config['reference'])

        #{vcfs: {db_snp: {path: '', annotation: []}, cosmic: {path: '', ann:[]}}}

        if 'vcfs' in self.run_config:
            for dbname, conf in self.run_config['vcfs'].items():
                sample_fpath = self.snpsift_annotate(dbname, conf, sample_fpath)

        sample_fpath = self.snpsift_db_nsfp(sample_fpath)

        if 'tracks' in self.run_config:
            for track in self.run_config['tracks']:
                sample_fpath = self.tracks(track, sample_fpath)

        sample_fpath = self.gatk(sample_fpath)
        sample_fpath = self.snpeff(sample_fpath)
        self.extract_fields(sample_fpath)


    def split_genotypes(self, sample_fpath, result_fpath):
        self.log_print('')
        self.log_print('*' * 70)
        self.log_print('Splitting genotypes.')

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
                        for alt in set(alts):
                            line = '\t'.join(tokens[:2] + ['.'] + [tokens[3]] + [alt] + tokens[5:]) + '\n'
                            out.write(line)
                    else:
                        line = '\t'.join(tokens[:2] + ['.'] + tokens[3:]) + '\n'
                        out.write(line)

        if self.run_config.get('save_intermediate'):
            return result_fpath
        else:
            os.remove(sample_fpath)
            os.rename(result_fpath, sample_fpath)
            return sample_fpath


def remove_quotes(str):
    if str and str[0] == '"':
        str = str[1:]
    if str and str[-1] == '"':
        str = str[:-1]
    return str


def main(args):
    if len(args) < 2:
        sys.stderr.write('Usage: python ' + __file__ + ' system_info_local.yaml run_info.yaml\n')
        exit(1)

    assert os.path.isfile(args[0]), args[0] + ' does not exist of is a directory.'
    assert os.path.isfile(args[1]), args[1] + ' does not exist of is a directory.'
    system_config = load(open(args[0]), Loader=Loader)
    run_config = load(open(args[1]), Loader=Loader)
    print('Loaded system config ' + args[0])
    print('Loaded run config ' + args[1])

    annotator = Annotator(system_config, run_config)

    print('Note: please, load modules before start:')
    print('   source /etc/profile.d/modules.sh')
    print('   module load java')
    print('   module load perl')
    print('Use "module load bcbio-nextgen" if you want to annotate with bed tracks.')
    # print ''
    # print 'In Waltham, run this as well:'
    # print '   export PATH=$PATH:/group/ngs/src/snpEff/snpEff3.5/scripts'
    # print '   export PERL5LIB=$PERL5LIB:/opt/az/local/bcbio-nextgen/stable/0.7.6/tooldir/lib/perl5/site_perl'

    annotator.annotate()


if __name__ == '__main__':
    main(sys.argv[1:])