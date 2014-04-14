#!/usr/bin/env python

import sys
if sys.version_info[:2] < (2, 5):
    exit('Python version 2.5 and higher is supported (you running ' +
         '.'.join(map(str, sys.version_info[:3])) + ')\n')

from distutils.version import LooseVersion
import os
from os.path import join, splitext, basename, realpath, isfile, getsize, dirname, exists
import subprocess
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


def remove_annotation(field_to_del, input_fpath):
    output_fpath = input_fpath + '_tmp'
    with open(input_fpath) as inp, open(output_fpath, 'w') as out:
        for l in inp:
            if field_to_del in l:
                l = l.lstrip()
                if l and l.startswith('##INFO='):
                    try:
                        if l.split('=', 1)[1].split(',', 1)[0].split('=')[1] == field_to_del:
                            continue
                    except:
                        self.log_exit('Incorrect VCF at line: ' + l)
                elif l.strip() and l.strip()[0] != '#':
                    fields = l.split('\t')
                    info_line = fields[7]
                    info_pairs = [attr.split('=') for attr in info_line.split(';')]
                    info_pairs = filter(lambda pair: pair[0] != field_to_del, info_pairs)
                    info_line = ';'.join('='.join(pair) if len(pair) == 2 else pair[0] for pair in info_pairs)
                    fields = fields[:7] + [info_line] + fields[8:]
                    l = '\t'.join(fields)
            out.write(l)
    os.rename(output_fpath, input_fpath)


class Annotator:
    def _set_up_dir(self):
        output_dir = realpath(self.run_config.get('output_dir', os.getcwd()))
        if not output_dir:
            exit('Result directory path is empty.')
        if isfile(output_dir):
            exit(output_dir + ' is a file.')
        if not exists(output_dir):
            try:
                os.mkdir(output_dir)
            except:
                exit(output_dir + ' does not exist.')
        self.output_dir = output_dir
        return output_dir


    def __init__(self, system_config_path, run_config_path):
        self.system_config = load(open(system_config_path), Loader=Loader)
        self.run_config = load(open(run_config_path), Loader=Loader)

        output_dir = self._set_up_dir()

        if 'log' in self.run_config:
            self.log = os.path.join(output_dir, 'log.txt')
            if os.path.isfile(self.log):
                os.remove(self.log)

            shutil.copy(run_config_path, output_dir)
        else:
            self.log = None

        self.log_print('Loaded system config ' + system_config_path)
        self.log_print('Loaded run config ' + run_config_path)
        self.log_print('')
        self.log_print('Writing into ' + output_dir)
        if self.log:
            self.log_print('Logging to ' + self.log)
        self.log_print('')

        if not which('java'):
            sys.stderr.write('* Warning: Java not found. You may want to run "module load java", '
                             'or better ". /group/ngs/bin/bcbio-prod.sh"\n\n')
        if not which('perl'):
            sys.stderr.write('* Warning: Perl not found. You may want to run "module load perl", '
                             'or better ". /group/ngs/bin/bcbio-prod.sh"\n\n')
        if not self._get_tool_cmdline('vcfannotate',
                                      extra_warn='You may want to load BCBio with ". /group/ngs/bin/bcbio-prod.sh"'):
            sys.stderr.write('* Warning: skipping annotation with bed tracks.\n')

        data = []
        if 'input' not in self.run_config:
            if 'file' not in self.run_config:
                self.log_exit('ERROR: Run config does not contain "input" section.')
            data.append({'vcf': realpath(self.run_config['file']),
                         'bam': self.run_config.get('bam'),
                         'bam_per_sample': None})
        else:
            for rec in self.run_config['input']:
                if 'vcf' not in rec:
                    self.log_exit('ERROR: Input section does not contain field "vcf".')
                data.append({'vcf': realpath(rec['vcf']),
                             'bam': rec.get('bam'),
                             'bam_per_sample': rec.get('bams')})

        self.data = []
        for rec in data:
            # Check VCF and BAMs existence
            inp_fpath = rec['vcf']
            if not file_exists(inp_fpath, 'Input file'):
                exit(1)
            bam_fpath = rec['bam']
            if bam_fpath:
                if not file_exists(bam_fpath, 'Bam file'):
                    exit(1)
            bams = rec['bam_per_sample']
            if bams:
                bam_fpaths = [fpath for sample, fpath in bams.items()]
                if bam_fpaths:
                    for bam_fpath in bam_fpaths:
                        if not file_exists(bam_fpath, 'Bam file'):
                            exit(1)

            # Move input VCF
            fname = os.path.basename(inp_fpath)
            base_name, ext = os.path.splitext(fname)

            if output_dir != os.path.realpath(os.path.dirname(inp_fpath)):
                if self.run_config.get('save_intermediate'):
                    new_fname = fname
                else:
                    new_fname = base_name + '.anno' + ext
                new_fpath = os.path.join(output_dir, new_fname)

                if os.path.exists(new_fpath):
                    os.remove(new_fpath)
                shutil.copyfile(inp_fpath, new_fpath)
                rec['vcf'] = new_fpath
            else:
                if not self.run_config.get('save_intermediate'):
                    new_fpath = join(output_dir, base_name + '.anno' + ext)
                    shutil.copyfile(inp_fpath, new_fpath)
                    rec['vcf'] = new_fpath
                else:
                    rec['vcf'] = inp_fpath
            inp_fpath = rec['vcf']
            basic_fields = next(l.strip()[1:].split() for l in open(inp_fpath)
                                if l.strip().startswith('#CHROM'))
            samples = basic_fields[9:]

            # Split by samples
            if bams or self.run_config.get('split_samples'):
                for sample in (bams.keys() if bams else []):
                    if sample not in samples:
                        self.log_exit('ERROR: sample ' + sample + ' is not in VCF ' + inp_fpath + '\n'
                                      'Available samples: ' + ', '.join(samples))
                for sample in samples:
                    bam_fpath = bams.get(sample) if bams else None
                    new_vcf = self.split_samples(inp_fpath, sample)
                    self.data.append({'name': sample, 'vcf': new_vcf, 'bam': bam_fpath})
            else:
                self.data.append({'vcf': inp_fpath, 'bam': bam_fpath,
                                  'name': samples[0] if len(samples) == 1 else None})


    def annotate(self):
        for rec in self.data:
            self.annotate_one(rec['vcf'], rec['bam'], rec.get('name'))


    def annotate_one(self, input_fpath, bam_fpath, sample_name=None):
        if not 'resources' in self.system_config:
            self.log_exit('"resources" section in system config required.')

        if sample_name:
            self.log_print('')
            self.log_print('')
            self.log_print('*' * 70)
            msg = '*' * 3 + ' Sample ' + sample_name + ' '
            self.log_print(msg + ('*' * (70 - len(msg)) if len(msg) < 70 else ''))
            self.log_print('VCF: ' + input_fpath)
            if bam_fpath:
                self.log_print('BAM: ' + bam_fpath)

        self.all_fields = []

        remove_annotation('EFF', input_fpath)

        if self.run_config.get('split_genotypes'):
            base_path, ext = os.path.splitext(input_fpath)
            result_fpath = base_path + '.split' + ext
            input_fpath = self.split_genotypes(input_fpath, result_fpath)
            self.log_print('Saved to ' + input_fpath)

        if self.run_config.get('rna'):
            input_fpath = self.process_rna(input_fpath)

        if self.run_config.get('ensemble'):
            input_fpath = self.process_ensemble(input_fpath)

        if not 'genome_build' in self.run_config:
            self.log_exit('Please, provide genome build (genome_build).')
        if not 'reference' in self.run_config:
            self.log_exit('Please, provide path to the reference file (reference).')
        if not file_exists(self.run_config['reference'], 'Reference'):
            exit(1)

        input_fpath = self.gatk(input_fpath, bam_fpath)
        if 'vcfs' in self.run_config:
            for dbname, conf in self.run_config['vcfs'].items():
                input_fpath = self.snpsift_annotate(dbname, conf, input_fpath)

        input_fpath = self.snpsift_db_nsfp(input_fpath)

        input_fpath = self.snpeff(input_fpath)

        if self.run_config.get('tracks'):
            for track in self.run_config['tracks']:
                input_fpath = self.tracks(track, input_fpath)

        input_fpath = self.filter_fields(input_fpath)

        tsv_fpath = self.extract_fields(input_fpath, sample_name)

        manual_tsv_fields = self.run_config.get('tsv_fields')
        if manual_tsv_fields:
            field_map = dict((rec.keys()[0], rec.values()[0]) for rec in manual_tsv_fields)
            tsv_fpath = self.rename_fields(tsv_fpath, field_map)
            self.log_print('Saved final TSV file with nice names to ' + tsv_fpath)

        if self.run_config.get('save_intermediate'):
            corr_tsv_fpath = correct_tabs(tsv_fpath)
            self.log_print('TSV file with dots to ' + corr_tsv_fpath)
            self.log_print('View with the commandline:  column -t ' + corr_tsv_fpath + ' | less -S')

        self.log_print('\nFinal VCF in ' + input_fpath)
        if self.log:
            print('Log in ' + self.log)


    def split_samples(self, input_fpath, sample):
        self.log_print('')
        self.log_print('-' * 70)
        self.log_print('Separating out sample ' + sample)
        self.log_print('-' * 70)
        executable = self._get_java_tool_cmdline('gatk')
        output_fpath = splitext(input_fpath)[0] + '.' + sample + '.vcf'
        ref_fpath = self.run_config['reference']
        cmd = '{executable} -nt 30 -R {ref_fpath} -T SelectVariants ' \
              '--variant {input_fpath} -o {output_fpath} -sn {sample}'.format(**locals())
        self._call_and_rename(cmd, input_fpath, suffix=sample, result_to_stdout=False,
                              rename=False)
        return output_fpath


    def log_print(self, msg=''):
        print(msg)
        if self.log:
            open(self.log, 'a').write(msg + '\n')


    def log_exit(self, msg=''):
        if self.log:
            open(self.log, 'a').write(msg + '\n')
        exit(msg)


    def log_err(self, msg=''):
        if self.log:
            open(self.log, 'a').write(msg + '\n')
        sys.stderr.write(msg + '\n')


    def _call_and_rename(self, cmdline, input_fpath, suffix, result_to_stdout=True, to_remove=None, rename=True):
        to_remove = to_remove or []
        basepath, ext = splitext(input_fpath)
        output_fpath = basepath + '.' + suffix + ext

        if self.run_config.get('save_intermediate') and \
                (self.run_config.get('reuse_intermediate') or self.run_config.get('reuse')):
            if isfile(output_fpath) and getsize(output_fpath) > 0:
                self.log_print(output_fpath + ' exists, reusing')
                return output_fpath

        self.log_print(cmdline)

        err_fpath = os.path.join(self.output_dir, 'annotate_py_err.tmp')
        to_remove.append(err_fpath)
        if err_fpath and isfile(err_fpath):
            os.remove(err_fpath)

        if self.run_config.get('verbose', True):
            proc = subprocess.Popen(
                cmdline.split(),
                stdout=open(output_fpath, 'w') if result_to_stdout else subprocess.PIPE,
                stderr=subprocess.STDOUT if not result_to_stdout else subprocess.PIPE)

            if proc.stdout:
                for line in iter(proc.stdout.readline, ''):
                    self.log_print('   ' + line.strip())
            elif proc.stderr:
                for line in iter(proc.stderr.readline, ''):
                    self.log_print('   ' + line.strip())

            ret_code = proc.wait()
            if ret_code != 0:
                for fpath in to_remove:
                    if fpath and isfile(fpath):
                        os.remove(fpath)
                self.log_exit('Command returned status ' + str(ret_code) +
                              ('. Log in ' + self.log if self.log else '.'))
        else:
            res = subprocess.call(
                cmdline.split(),
                stdout=open(output_fpath, 'w') if result_to_stdout else open(err_fpath, 'a'),
                stderr=open(err_fpath, 'a'))
            if res != 0:
                with open(err_fpath) as err:
                    self.log_print('')
                    self.log_print(err.read())
                    self.log_print('')
                for fpath in to_remove:
                    if fpath and isfile(fpath):
                        os.remove(fpath)
                self.log_exit('Command returned status ' + str(res) +
                              ('. Log in ' + self.log if self.log else '.'))
            else:
                if self.log:
                    with open(err_fpath) as err, open(self.log, 'a') as log:
                        log.write('')
                        log.write(err.read())
                        log.write('')
        for fpath in to_remove:
            if fpath and isfile(fpath):
                os.remove(fpath)

        if rename and not self.run_config.get('save_intermediate'):
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


    def _get_tool_cmdline(self, tool_name, extra_warn=''):
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
                if extra_warn:
                    self.log_err(extra_warn)
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
        self.log_print('-' * 70)
        self.log_print('Annotate with ' + dbname)
        self.log_print('-' * 70)

        executable = self._get_java_tool_cmdline('snpsift')
        db_path = conf.get('path')
        annotations = conf.get('annotations', [])
        self.all_fields.extend(annotations)
        anno_line = ('-info ' + ','.join(annotations)) if annotations else ''

        cmdline = '{executable} annotate -v {anno_line} {db_path} {input_fpath}'.format(**locals())

        output_fpath = self._call_and_rename(cmdline, input_fpath, dbname, result_to_stdout=True)

        corr_out_fpath = output_fpath + '_tmp'
        with open(output_fpath) as inp, open(corr_out_fpath, 'w') as out:
            for l in inp:
                if l and not l.strip().startswith('#'):
                    l = l.replace(' ', '_')
                out.write(l)
        os.rename(corr_out_fpath, output_fpath)
        return output_fpath


    def snpsift_db_nsfp(self, input_fpath):
        if 'db_nsfp' not in self.run_config:
            return input_fpath

        self.log_print('')
        self.log_print('-' * 70)
        self.log_print('DB SNFP')
        self.log_print('-' * 70)

        executable = self._get_java_tool_cmdline('snpsift')

        db_path = self.run_config['db_nsfp'].get('path')
        if not db_path:
            self.log_exit('Please, provide a path to db nsfp file in run_config.')

        annotations = self.run_config['db_nsfp'].get('annotations', [])
        self.all_fields.extend(['dbNSFP_' + ann for ann in annotations])
        ann_line = ('-f ' + ','.join(annotations)) if annotations else ''

        cmdline = '{executable} dbnsfp {ann_line} -v {db_path} {input_fpath}'.format(**locals())

        return self._call_and_rename(cmdline, input_fpath, 'db_nsfp', result_to_stdout=True)


    def snpeff(self, input_fpath):
        if 'snpeff' not in self.run_config:
            return input_fpath

        self.log_print('')
        self.log_print('-' * 70)
        self.log_print('SnpEff')
        self.log_print('-' * 70)

        self.all_fields.extend([
            "EFF[*].EFFECT", "EFF[*].IMPACT", "EFF[*].FUNCLASS", "EFF[*].CODON",
            "EFF[*].AA", "EFF[*].AA_LEN", "EFF[*].GENE", "EFF[*].CODING",
            "EFF[*].TRID", "EFF[*].RANK"])

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
        field_name = splitext(basename(track_path))[0]

        self.log_print('')
        self.log_print('-' * 70)
        self.log_print('Intersecting with ' + field_name)
        self.log_print('-' * 70)

        toolpath = self._get_tool_cmdline('vcfannotate')
        if not toolpath:
            self.log_err('WARNING: Skipping annotation with tracks: vcfannotate '
                         'executable not found, you probably need to '
                         'run the commandline:  . /group/ngs/bin/bcbio-prod.sh"')
            return

        self.all_fields.append(field_name)

        cmdline = 'vcfannotate -b {track_path} -k {field_name} {input_fpath}'.format(**locals())

        out_fpath = self._call_and_rename(cmdline, input_fpath, field_name, result_to_stdout=True)

        # Set TRUE or FALSE for tracks
        corr_out_fpath = out_fpath + '_tmp'
        with open(out_fpath) as inp, open(corr_out_fpath, 'w') as out:
            for l in inp:
                if field_name in l:
                    l = l.lstrip()
                    if l.strip() and l.strip()[0] != '#':
                        fields = l.split('\t')
                        info_line = fields[7]
                        info_pairs = [attr.split('=') for attr in info_line.split(';')]
                        info_pairs = [[pair[0], ('TRUE' if pair[1] else 'FALSE')]
                                      if pair[0] == field_name and len(pair) > 1
                                      else pair for pair in info_pairs]
                        info_line = ';'.join('='.join(pair) if len(pair) == 2 else pair[0] for pair in info_pairs)
                        fields = fields[:7] + [info_line] + fields[8:]
                        l = '\t'.join(fields)
                out.write(l)
        os.rename(corr_out_fpath, out_fpath)
        return out_fpath


    def _get_gatk_version(self):
        cmdline = self._get_java_tool_cmdline('gatk') + ' -version'

        version = None
        with subprocess.Popen(cmdline,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT,
                              shell=True).stdout as stdout:
            out = stdout.read().strip()
            last_line = out.split('\n')[-1].strip()
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
            self.log_print('WARNING: could not determine Gatk version, using 1.0')
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

    def gatk(self, input_fpath, bam_fpath):
        if 'gatk' not in self.run_config:
            return input_fpath

        self.log_print('')
        self.log_print('-' * 70)
        self.log_print('GATK')
        self.log_print('-' * 70)

        executable = self._get_java_tool_cmdline('gatk')

        base_name, ext = os.path.splitext(input_fpath)
        output_fpath = base_name + '.gatk' + ext

        ref_fpath = self.run_config['reference']

        cmdline = ('{executable} -nt 20 -R {ref_fpath} -T VariantAnnotator -o {output_fpath}'
                   ' --variant {input_fpath}').format(**locals())
        if bam_fpath:
            cmdline += ' -I ' + bam_fpath

        GATK_ANNOS_DICT = {
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

        annotations = self.run_config['gatk'].get('annotations', [])
        self.all_fields.extend(GATK_ANNOS_DICT.get(ann) for ann in annotations)
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


    def rename_fields(self, tsv_fpath, field_map):
        self.log_print('')
        self.log_print('-' * 70)
        self.log_print('Renaming fields.')
        self.log_print('-' * 70)

        with open(tsv_fpath) as f:
            first_line = f.readline()[1:]
        fields = first_line.split()
        new_fields = [field_map.get(f, f) for f in fields]
        new_line = '\t'.join(new_fields)

        out_fpath = splitext(tsv_fpath)[0] + '.renamed' + '.tsv'
        with open(out_fpath, 'w') as out:
            out.write(new_line + '\n')
            with open(tsv_fpath) as f:
                for i, l in enumerate(f):
                    if i != 0:
                        out.write(l)

        if self.run_config.get('save_intermediate'):
            return out_fpath
        else:
            os.remove(tsv_fpath)
            os.rename(out_fpath, tsv_fpath)
            return tsv_fpath


    def extract_fields(self, vcf_fpath, sample_name=None):
        tsv_fpath = os.path.splitext(vcf_fpath)[0] + '.tsv'

        if self.run_config.get('save_intermediate') and \
            (self.run_config.get('reuse_intermediate') or self.run_config.get('reuse')):
            if isfile(tsv_fpath) and getsize(tsv_fpath) > 0:
                self.log_print(tsv_fpath + ' exists, reusing')
                return tsv_fpath
        else:
            if isfile(tsv_fpath):
                os.remove(tsv_fpath)

        first_line = next(l.strip()[1:].split() for l in open(vcf_fpath) if l.strip().startswith('#CHROM'))
        basic_fields = [f for f in first_line[:9] if f != 'INFO' and f != 'FORMAT' and f != sample_name]
        manual_annots = filter(lambda f: f and f != 'ID', self.all_fields)

        manual_tsv_fields = self.run_config.get('tsv_fields')
        if manual_tsv_fields:
            fields = [rec.keys()[0] for rec in manual_tsv_fields]
        else:
            fields = (basic_fields + manual_annots + self.run_config.get('additional_tsv_fields', []))
        if not fields:
            return

        self.log_print('')
        self.log_print('-' * 70)
        self.log_print('Extracting fields')
        self.log_print('-' * 70)

        format_fields = set()
        tmp_vcf = None
        # Split FORMAT
        if sample_name:
            tmp_vcf = vcf_fpath + '_tmp'
            with open(vcf_fpath) as inp, open(tmp_vcf, 'w') as out:
                for l in inp:
                    if not l or l.startswith('#'):
                        out.write(l)
                        continue
                    vals = l.strip().split('\t')
                    if len(vals) <= 9:
                        out.write(l)
                        continue
                    info = vals[7]
                    format_fields = vals[8].split(':')
                    sample_fields = vals[9].split(':')
                    for f, s in zip(format_fields, sample_fields):
                        if f not in ['DP', 'MQ']:
                            f = 'gt_' + f
                            info += ';' + f + '=' + s
                            format_fields.append(f)
                    l = '\t'.join(vals[:7] + [info])
                    out.write(l + '\n')

        anno_line = ' '.join(fields + list(format_fields))
        snpsift_cmline = self._get_java_tool_cmdline('snpsift')
        vcfoneperline_cmline = self._get_script_cmdline_template('perl', 'vcfoneperline') % ''
        cmdline = vcfoneperline_cmline + ' | ' + snpsift_cmline + ' extractFields - ' + anno_line
        self.log_print(cmdline)
        res = subprocess.call(cmdline,
                              stdin=open(tmp_vcf or vcf_fpath),
                              stdout=open(tsv_fpath, 'w'), shell=True)
        # if tmp_vcf:
        #     os.remove(tmp_vcf)

        self.log_print('')
        if res != 0:
            self.log_exit('Command returned status ' + str(res) +
                          ('. Log in ' + self.log if self.log else '.'))
        self.log_print('Saved TSV file to ' + tsv_fpath)
        return tsv_fpath


    def process_rna(self, input_fpath):
        self.log_print('')
        self.log_print('-' * 70)

        cmdline = self._get_tool_cmdline('vcf-subset')
        if not cmdline:
            self.log_err('WARNING: vcf-subset executable not found in PATH system config, '
                         'skipping subset process for ensemble VCF.')
            return

        input_fname = os.path.basename(input_fpath)
        base_name, ext = os.path.splitext(input_fname)
        name = base_name.replace('-ensemble', '')
        cmdline = 'vcf-subset -c {name} -e {input_fpath}'.format(**locals())

        return self._call_and_rename(cmdline, input_fpath, '.ensm', result_to_stdout=True)


    def process_ensemble(self, input_fpath):
        self.log_print('')
        self.log_print('-' * 70)
        self.log_print('Filtering ensemble reject lines.')
        self.log_print('-' * 70)

        base_path, ext = os.path.splitext(input_fpath)
        output_fpath = base_path + '.pass' + ext

        if self.run_config.get('save_intermediate') and \
            (self.run_config.get('reuse_intermediate') or self.run_config.get('reuse')):
            if isfile(output_fpath) and getsize(output_fpath) > 0:
                self.log_print(output_fpath + ' exists, reusing')
                return output_fpath

        with open(input_fpath) as input_f, open(output_fpath, 'w') as pass_f:
            for line in input_f.readlines():
                if 'REJECT' not in line:
                    pass_f.write(line)
        if self.run_config.get('save_intermediate'):
            return output_fpath
        else:
            os.remove(input_fpath)
            os.rename(output_fpath, input_fpath)
            return input_fpath


    def split_genotypes(self, input_fpath, output_fpath):
        self.log_print('')
        self.log_print('-' * 70)
        self.log_print('Splitting genotypes.')
        self.log_print('-' * 70)

        if self.run_config.get('save_intermediate') and \
            (self.run_config.get('reuse_intermediate') or self.run_config.get('reuse')):
            if isfile(output_fpath) and getsize(output_fpath) > 0:
                self.log_print(output_fpath + ' exists, reusing')
                return output_fpath

        with open(input_fpath) as vcf, open(output_fpath, 'w') as out:
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
            return output_fpath
        else:
            os.remove(input_fpath)
            os.rename(output_fpath, input_fpath)
            return input_fpath


    def filter_fields(self, input_fpath):
        self.log_print('')
        self.log_print('-' * 70)
        self.log_print('Filtering incorrect fields.')
        self.log_print('-' * 70)

        output_fpath = splitext(input_fpath)[0] + '.filt.vcf'

        if self.run_config.get('save_intermediate') and \
            (self.run_config.get('reuse_intermediate') or self.run_config.get('reuse')):
            if isfile(output_fpath) and getsize(output_fpath) > 0:
                self.log_print(output_fpath + ' exists, reusing')
                return output_fpath

        with open(input_fpath) as inp, open(output_fpath, 'w') as out:
            for l in inp:
                if l.strip() and l.strip()[0] != '#':
                    if ',.' in l or '.,' in l:
                        fields = l.split('\t')
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
                        l = '\t'.join(fields)
                out.write(l)

        if self.run_config.get('save_intermediate'):
            return output_fpath
        else:
            os.remove(input_fpath)
            os.rename(output_fpath, input_fpath)
            return input_fpath


    # def split_format(self, tsv_fpath, sample_name):
    #     output_fpath = splitext(tsv_fpath)[0] + '_tmp'
    #
    #     with open(tsv_fpath) as inp, open(output_fpath, 'w') as out:
    #         for i, l in enumerate(inp):
    #             if i != 0 and l.strip() and l.strip()[0] != '#':
    #                 fields = l.split('\t')
    #                 info_line = fields[7]
    #                 info_pairs = [attr.split('=') for attr in info_line.split(';')]
    #                 new_info_pairs = []
    #                 for p in info_pairs:
    #                     if len(p) == 2:
    #                         if p[1].endswith(',.'):
    #                             p[1] = p[1][:-2]
    #                         if p[1].startswith('.,'):
    #                             p[1] = p[1][2:]
    #                         new_info_pairs.append('='.join(p))
    #                 info_line = ';'.join(new_info_pairs)
    #                 fields = fields[:7] + [info_line] + fields[8:]
    #                 l = '\t'.join(fields)
    #             out.write(l)
    #
    #     os.remove(tsv_fpath)
    #     os.rename(output_fpath, tsv_fpath)
    #     return tsv_fpath


def correct_tabs(tsv_fpath):
    out_fpath = splitext(tsv_fpath)[0] + '.tabs' + '.tsv'
    with open(out_fpath, 'w') as out, open(tsv_fpath) as inp:
        for l in inp:
            while '\t\t' in l:
                l = l.replace('\t\t', '\t.\t')
            out.write(l)

    return out_fpath


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
        sys.stderr.write('Notice: using system_info_rask.yaml as a default tools configutation file.\n\n')
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