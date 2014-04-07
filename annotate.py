#!/usr/bin/env python

from distutils.version import LooseVersion
import os
from os.path import join, splitext, basename, realpath, isfile, getsize, dirname, exists
import subprocess
import sys
import shutil
from yaml import load
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader


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


def file_exists(fpath, description=''):
    if not fpath:
        sys.stderr.write((description + ': f' if description else 'F') + 'ile name is empty.\n')
        return False
    if not exists(fpath):
        sys.stderr.write((description + ': ' if description else '') + fpath + ' does not exist.\n')
        return False
    if not isfile(fpath):
        sys.stderr.write((description + ': ' if description else '') + fpath + ' not a file.\n')
        return False
    if getsize(fpath) <= 0:
        sys.stderr.write((description + ': ' if description else '') + fpath + ' is empty.\n')
        return False
    return True


class Annotator:
    def _set_up_dir(self):
        result_dir = realpath(self.run_config.get('output_dir', os.getcwd()))
        if not result_dir:
            exit('Result directory path is empty.')
        if isfile(result_dir):
            exit(result_dir + ' is a file.')
        if not exists(result_dir):
            try:
                os.mkdir(result_dir)
            except:
                exit(result_dir + ' does not exist.')
        self.run_config['output_dir'] = result_dir
        return result_dir


    def __init__(self, system_config_path, run_config_path):
        self.system_config = load(open(system_config_path), Loader=Loader)
        self.run_config = load(open(run_config_path), Loader=Loader)
        self.all_annotations = []

        if not 'file' in self.run_config:
            exit('Run config does not contain field "file".')
        sample_fpath = self.run_config['file']
        self.sample_fpath = realpath(sample_fpath)
        if not file_exists(self.sample_fpath, 'Sample file'):
            exit(1)

        sample_fname = os.path.basename(self.sample_fpath)
        sample_basename, ext = os.path.splitext(sample_fname)

        output_dir = self._set_up_dir()

        if output_dir != os.path.realpath(os.path.dirname(self.sample_fpath)):
            if self.run_config.get('save_intermediate'):
                new_sample_fname = sample_fname
            else:
                new_sample_fname = sample_basename + '.anno' + ext

            new_sample_fpath = os.path.join(output_dir, new_sample_fname)

            if os.path.exists(new_sample_fpath):
                os.remove(new_sample_fpath)
            shutil.copyfile(self.sample_fpath, new_sample_fpath)
            self.sample_fpath = new_sample_fpath
        else:
            if not self.run_config.get('save_intermediate'):
                new_sample_fpath = join(output_dir, sample_basename + '.anno' + ext)
                shutil.copyfile(self.sample_fpath, new_sample_fpath)
                self.sample_fpath = new_sample_fpath

        if 'log' in self.run_config:
            self.run_config['log'] = os.path.join(os.path.dirname(self.sample_fpath), sample_basename + '.log')
            if os.path.isfile(self.run_config['log']):
                os.remove(self.run_config['log'])

        self.log_print('Loaded system config ' + system_config_path)
        self.log_print('Loaded run config ' + run_config_path)
        self.log_print('')
        self.log_print('Writing into ' + output_dir)
        if 'log' in self.run_config:
            self.log_print('Logging to ' + self.run_config['log'])

        if not which('java'):
            sys.stderr.write('\nWARNING: Please, run "module load java"\n')
        if not which('perl'):
            sys.stderr.write('\nWARNING: Please, run "module load perl"\n')
        if not self._get_tool_cmdline('vcfannotate'):
            sys.stderr.write('\nWARNING: Please, run "module load bcbio-nextgen" '
                             'if you want to annotate with bed tracks.\n')


    def annotate(self):
        if not 'resources' in self.system_config:
            exit('"resources" section in system config required.')

        sample_fpath = self.sample_fpath

        if self.run_config.get('split_genotypes'):
            sample_basepath, ext = os.path.splitext(sample_fpath)
            result_fpath = sample_basepath + '.split' + ext
            sample_fpath = self.split_genotypes(sample_fpath, result_fpath)
            self.log_print('Saved to ' + sample_fpath)

        if self.run_config.get('split_samples'):
            #TODO
            pass

        if self.run_config.get('rna'):
            sample_fpath = self.process_rna(sample_fpath)

        if self.run_config.get('ensemble'):
            sample_fpath = self.process_ensemble(sample_fpath)

        if not 'genome_build' in self.run_config:
            exit('Please, provide genome build (genome_build).')
        if not 'reference' in self.run_config:
            exit('Please, provide path to the reference file (reference).')
        if not file_exists(self.run_config['reference'], 'Reference'):
            exit(1)

        sample_fpath = self.gatk(sample_fpath)
        if 'vcfs' in self.run_config:
            for dbname, conf in self.run_config['vcfs'].items():
                sample_fpath = self.snpsift_annotate(dbname, conf, sample_fpath)
        sample_fpath = self.snpsift_db_nsfp(sample_fpath)
        sample_fpath = self.snpeff(sample_fpath)
        if self.run_config.get('tracks'):
            for track in self.run_config['tracks']:
                sample_fpath = self.tracks(track, sample_fpath)
        self.extract_fields(sample_fpath)

        self.log_print('\nFinal VCF in ' + sample_fpath)
        if 'log' in self.run_config:
            print('Log in ' + self.run_config['log'])


    def log_print(self, msg=''):
        print(msg)
        if 'log' in self.run_config:
            open(self.run_config['log'], 'a').write(msg + '\n')


    def log_exit(self, msg=''):
        if 'log' in self.run_config:
            open(self.run_config['log'], 'a').write(msg + '\n')
        exit(msg)


    def log_err(self, msg=''):
        if 'log' in self.run_config:
            open(self.run_config['log'], 'a').write(msg + '\n')
        sys.stderr.write(msg + '\n')


    def _call_and_rename(self, cmdline, input_fpath, suffix, result_to_stdout=True, to_remove=None):
        to_remove = to_remove or []
        basepath, ext = splitext(input_fpath)
        output_fpath = basepath + '.' + suffix + ext

        if self.run_config.get('save_intermediate') and \
                (self.run_config.get('reuse_intermediate') or self.run_config.get('reuse')):
            if isfile(output_fpath) and getsize(output_fpath) > 0:
                self.log_print(output_fpath + ' exists, reusing')
                return output_fpath

        self.log_print(cmdline)

        err_fpath = os.path.join(os.path.dirname(self.sample_fpath)) + 'annotate_py_err.tmp'
        to_remove.append(err_fpath)

        if self.run_config.get('verbose', True):
            res = subprocess.call(
                cmdline.split(),
                stdout=open(output_fpath, 'w') if result_to_stdout else subprocess.PIPE,
                stderr=subprocess.PIPE)
            if res != 0:
                for fpath in to_remove:
                    if fpath and isfile(fpath):
                        os.remove(fpath)
                self.log_exit('Command returned status ' + str(res) +
                              ('. Log in ' + self.run_config['log'] if 'log' in self.run_config else '.'))
        else:
            res = subprocess.call(
                cmdline.split(),
                stdout=open(output_fpath, 'w') if result_to_stdout else open(self.run_config['log'], 'a'),
                stderr=open(err_fpath, 'w'))
            if res != 0:
                with open(err_fpath) as err:
                    self.log_print('')
                    self.log_print(err.read())
                    self.log_print('')
                for fpath in to_remove:
                    if fpath and isfile(fpath):
                        os.remove(fpath)
                self.log_exit('Command returned status ' + str(res) +
                              ('. Log in ' + self.run_config['log'] if 'log' in self.run_config else '.'))
            else:
                if 'log' in self.run_config:
                    with open(err_fpath) as err, open(self.run_config['log'], 'a') as log:
                        log.write('')
                        log.write(err.read())
                        log.write('')
        for fpath in to_remove:
            if fpath and isfile(fpath):
                os.remove(fpath)

        if not self.run_config.get('save_intermediate'):
            os.remove(input_fpath)
            os.rename(output_fpath, input_fpath)
            self.log_print('Saved to ' + input_fpath)
            return input_fpath
        else:
            self.log_print('Saved to ' + output_fpath)
            return output_fpath


    def _get_java_tool_cmdline(self, name):
        cmdline_template = self._get_script_cmdline_template('java', name)
        jvm_opts = self.system_config['resources'][name].get('jvm_opts', [])
        return cmdline_template % (' '.join(jvm_opts) + ' -jar')


    def _get_script_cmdline_template(self, executable, script_name):
        if not which(executable):
            exit(executable + ' executable required, maybe you need '
                 'to run "module load ' + executable + '"?')
        if 'resources' not in self.system_config:
            self.log_exit('System config yaml must contain resources section with '
                          + script_name + ' path.')
        if script_name not in self.system_config['resources']:
            self.log_exit('System config resources section must contain '
                          + script_name + ' info (with a path to the tool).')
        tool_config = self.system_config['resources'][script_name]
        if 'path' not in tool_config:
            self.log_exit(script_name + ' section in the system config must contain a path to the tool.')
        tool_path = tool_config['path']
        if not file_exists(tool_path, script_name):
            exit(1)
        return executable + ' %s ' + tool_path


    def _get_tool_cmdline(self, tool_name):
        tool_path = which(tool_name) or None

        if not 'resources' in self.system_config \
                or tool_name not in self.system_config['resources']\
                or 'path' not in self.system_config['resources'][tool_name]:
            if tool_path:
                return tool_path
            else:
                self.log_err(tool_name + ' executable was not found. '
                             'You can either specify path in the system config, or load into your '
                             'PATH environment variable.')
                return None

        tool_path = self.system_config['resources'][tool_name]['path']
        if file_exists(tool_path, tool_name):
            return tool_path
        else:
            self.log_err(tool_path + ' for ' + tool_name +
                         ' does not exist or is not a file.')
            return None


    def snpsift_annotate(self, dbname, conf, input_fpath):
        self.log_print('')
        self.log_print('*' * 70)

        executable = self._get_java_tool_cmdline('snpsift')
        db_path = conf.get('path')
        annotations = conf.get('annotations', [])
        self.all_annotations.extend(annotations)
        anno_line = ('-info ' + ','.join(annotations)) if annotations else ''

        cmdline = '{executable} annotate -v {anno_line} {db_path} {input_fpath}'.format(**locals())

        return self._call_and_rename(cmdline, input_fpath, dbname, result_to_stdout=True)


    def snpsift_db_nsfp(self, input_fpath):
        if 'db_nsfp' not in self.run_config:
            return input_fpath

        self.log_print('')
        self.log_print('*' * 70)

        executable = self._get_java_tool_cmdline('snpsift')

        db_path = self.run_config['db_nsfp'].get('path')
        if not db_path:
            exit('Please, provide a path to db nsfp file in run_config.')

        annotations = self.run_config['db_nsfp'].get('annotations', [])
        self.all_annotations.extend(annotations)
        ann_line = ('-f ' + ','.join(annotations)) if annotations else ''

        cmdline = '{executable} dbnsfp {ann_line} -v {db_path} {input_fpath}'.format(**locals())

        return self._call_and_rename(cmdline, input_fpath, 'db_nsfp', result_to_stdout=True)


    def snpeff(self, input_fpath):
        if 'snpeff' not in self.run_config:
            return input_fpath

        self.log_print('')
        self.log_print('*' * 70)

        executable = self._get_java_tool_cmdline('snpeff')
        ref_name = self.run_config['genome_build']
        db_path = self.run_config['snpeff'].get('path')
        assert db_path, 'Please, provide a path to db nsfp file in run_config.'

        cmdline = ('{executable} eff -dataDir {db_path} -noStats -noLog -1 '
                   '-i vcf -o vcf {ref_name} {input_fpath}').format(**locals())

        if self.run_config['snpeff'].get('clinical_reporting') or \
                self.run_config['snpeff'].get('canonical'):
            cmdline += ' -canon -hgvs '

        if self.run_config['snpeff'].get('cancer'):
            cmdline += ' -cancer '

        return self._call_and_rename(cmdline, input_fpath, 'snpEff', result_to_stdout=True)


    def tracks(self, track_path, input_fpath):
        self.log_print('')
        self.log_print('*' * 70)

        toolpath = self._get_tool_cmdline('vcfannotate')
        if not toolpath:
            self.log_err('WARNING: Skipping annotation with tracks: vcfannotate '
                         'executable not found, you probably need to run "module load bcbio-nextgen"')
            return

        field_name = splitext(basename(track_path))[0]
        self.all_annotations.append(field_name)

        cmdline = 'vcfannotate -b {track_path} -k {field_name} {input_fpath}'.format(**locals())

        return self._call_and_rename(cmdline, input_fpath, field_name, result_to_stdout=True)


    def _get_gatk_version(self):
        cmdline = self._get_java_tool_cmdline('gatk') + ' -version'

        version = None
        with subprocess.Popen(cmdline,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT,
                              shell=True).stdout as stdout:
            out = stdout.read().strip()
            print 'out =', out
            last_line = out.split('\n')[-1].strip()
            print 'out = stdout.read().strip() =', last_line
            # versions earlier than 2.4 do not have explicit version command,
            # parse from error output from GATK
            if out.find("ERROR") >= 0:
                flag = "The Genome Analysis Toolkit (GATK)"
                for line in last_line.split("\n"):
                    if line.startswith(flag):
                        version = line.split(flag)[-1].split(",")[0].strip()
            else:
                version = last_line
        if not version:
            self.log_print('Warning: could not determine Gatk version, using 1.0')
            return '1.0'
        if version.startswith("v"):
            version = version[1:]
        return version


    def _gatk_major_version(self):
        """Retrieve the GATK major version, handling multiple GATK distributions.

        Has special cases for GATK nightly builds, Appistry releases and
        GATK prior to 2.3.
        """
        full_version = self._get_gatk_version()
        # Working with a recent version if using nightlies
        if full_version.startswith("nightly-"):
            return "2.8"
        parts = full_version.split("-")
        if len(parts) == 4:
            appistry_release, version, subversion, githash = parts
        elif len(parts) == 3:
            version, subversion, githash = parts
        # version was not properly implemented in earlier GATKs
        else:
            version = "2.3"
        if version.startswith("v"):
            version = version[1:]
        return version


    def _gatk_type(self):
        """Retrieve type of GATK jar, allowing support for older GATK lite.
        Returns either `lite` (targeting GATK-lite 2.3.9) or `restricted`,
        the latest 2.4+ restricted version of GATK.
        """
        if LooseVersion(self._gatk_major_version()) > LooseVersion("2.3"):
            return "restricted"
        else:
            return "lite"

    def gatk(self, input_fpath):
        if 'gatk' not in self.run_config:
            return input_fpath

        self.log_print('')
        self.log_print('*' * 70)

        executable = self._get_java_tool_cmdline('gatk')

        base_name, ext = os.path.splitext(input_fpath)
        output_fpath = base_name + '.gatk' + ext

        ref_fpath = self.run_config['reference']

        cmdline = ('{executable} -R {ref_fpath} -T VariantAnnotator -o {output_fpath}'
                   ' --variant {input_fpath}').format(**locals())
        bam = self.run_config.get('bam')
        if bam:
            if not isfile(bam):
                self.log_exit('Error. Not such file: ' + bam)
            cmdline += ' -I ' + bam

        annotations = self.run_config['gatk'].get('annotations', [])
        self.all_annotations.extend(annotations)
        for ann in annotations:
            if ann == 'DepthOfCoverage' and self._gatk_type() == 'restricted':
                self.log_print('Notice: in the restricted Gatk version, DepthOfCoverage is renamed to Coverage. '
                               'Using the name Coverage.\n')
                ann = 'Coverage'
            if ann == 'Coverage' and self._gatk_type() == 'lite':
                self.log_print('Notice: in the lite Gatk version, the Coverage annotation goes by '
                               'name of DepthOfCoverage. '
                               'In the system config, the lite version of Gatk is specified; using DepthOfCoverage.\n')
                ann = 'DepthOfCoverage'
            cmdline += " -A " + ann

        res = self._call_and_rename(cmdline, input_fpath, 'gatk', result_to_stdout=False,
                                    to_remove=[output_fpath + '.idx', input_fpath + '.idx'])
        return res


    def extract_fields(self, input_fpath):
        fields = self.all_annotations + ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

        # if 'tsv_fields' not in self.run_config:
        #     return
        # fields = self.run_config.get('tsv_fields', [])
        if not fields:
            return
        anno_line = ' '.join(fields)

        self.log_print('')
        self.log_print('*' * 70)

        snpsift_cmline = self._get_java_tool_cmdline('snpsift')
        vcfoneperline_cmline = self._get_script_cmdline_template('perl', 'vcfoneperline') % ''

        cmdline = vcfoneperline_cmline + ' | ' + snpsift_cmline + ' extractFields - ' + anno_line

        basepath, ext = os.path.splitext(input_fpath)
        tsv_fpath = basepath + '.tsv'
        if isfile(tsv_fpath):
            os.remove(tsv_fpath)

        self.log_print(cmdline)
        res = subprocess.call(cmdline,
                              stdin=open(input_fpath),
                              stdout=open(tsv_fpath, 'w'),
                              shell=True)
        if res != 0:
            self.log_print('Command returned status ' + str(res) +
                           ('. Log in ' + self.run_config['log'] if 'log' in self.run_config else '.'))
            exit(1)
            # return input_fpath
        else:
            self.log_print('Saved TSV file to ' + tsv_fpath)


    def process_rna(self, sample_fpath):
        self.log_print('')
        self.log_print('*' * 70)

        cmdline = self._get_tool_cmdline('vcf-subset')
        if not cmdline:
            self.log_err('WARNING: vcf-subset executable not found in PATH system config, '
                         'skipping subset process for ensemble VCF.')
            return

        sample_fname = os.path.basename(sample_fpath)
        sample_basename, ext = os.path.splitext(sample_fname)
        name = sample_basename.replace('-ensemble', '')
        cmdline = 'vcf-subset -c {name} -e {sample_fpath}'.format(**locals())

        return self._call_and_rename(cmdline, sample_fpath, '.ensm', result_to_stdout=True)


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
                            line = '\t'.join(tokens[:2] + ['.'] + [tokens[3]] + [alt] + tokens[5:8]) + '\n'
                            out.write(line)
                    else:
                        line = '\t'.join(tokens[:2] + ['.'] + tokens[3:8]) + '\n'
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
    if len(args) < 1:
        exit('Usage: python ' + __file__ + ' system_info.yaml run_info.yaml\n'
             '    or python ' + __file__ + ' run_info.yaml')

    if len(args) == 1:
        run_config_path = args[0]
        system_config_path = join(dirname(realpath(__file__)), 'system_info_rask.yaml')
        sys.stderr.write('Warning: Using system_info_rask.yaml as a default tools configutation file.\n')
    else:
        system_config_path = args[0]
        run_config_path = args[1]

    if not os.path.isfile(system_config_path):
        exit(system_config_path + ' does not exist or is a directory.\n')
    if not os.path.isfile(run_config_path):
        exit(run_config_path + ' does not exist or is a directory.\n')

    if not system_config_path.endswith('.yaml'):
        sys.stderr.write(system_config_path + ' does not end with .yaml, maybe incorrect parameter?\n')
    if not run_config_path.endswith('.yaml'):
        sys.stderr.write(run_config_path + ' does not end with .yaml, maybe incorrect parameter?\n')

    annotator = Annotator(system_config_path, run_config_path)

    # print ''
    # print 'In Waltham, run this as well:'
    # print '   export PATH=$PATH:/group/ngs/src/snpEff/snpEff3.5/scripts'
    # print '   export PERL5LIB=$PERL5LIB:/opt/az/local/bcbio-nextgen/stable/0.7.6/tooldir/lib/perl5/site_perl'

    try:
        annotator.annotate()
    except KeyboardInterrupt:
        exit(1)


if __name__ == '__main__':
    main(sys.argv[1:])