#!/usr/bin/env python

import sys
if sys.version_info[:2] < (2, 5):
    exit('Python version 2.5 and higher is supported (you running ' +
         '.'.join(map(str, sys.version_info[:3])) + ')\n')

import subprocess
import shutil
import os

from os.path import join, splitext, basename, realpath, isfile, getsize, dirname, exists
from distutils.version import LooseVersion
from yaml import load
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

from src.bcbio_utils import which, splitext_plus, add_suffix, file_exists
from src.transaction import file_transaction


def verify_file(fpath, description=''):
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
    output_fpath = add_suffix(input_fpath, 'tmp')
    with open(input_fpath) as inp, file_transaction(output_fpath) as tx_out:
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
                    info_line = ';'.join('='.join(pair) if len(pair) == 2
                                         else pair[0] for pair in info_pairs)
                    fields = fields[:7] + [info_line] + fields[8:]
                    l = '\t'.join(fields)
            tx_out.write(l)
    os.rename(output_fpath, input_fpath)


class Annotator:
    def introduce_step(self, name):
        self.log_print('')
        self.log_print('-' * 70)
        self.log_print(name)
        self.log_print('-' * 70)


    def _set_up_dir(self):
        output_dir = realpath(self.run_cnf.get('output_dir', os.getcwd()))
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


    def _set_up_log(self):
        self.log = os.path.join(self.output_dir, 'log.txt')
        if os.path.isfile(self.log):
            os.remove(self.log)


    def _check_reference(self):
        if not 'genome_build' in self.run_cnf:
            self.log_exit('Please, provide genome build (genome_build).')
        if not 'reference' in self.run_cnf:
            self.log_exit('Please, provide path to the reference file (reference).')
        if not verify_file(self.run_cnf['reference'], 'Reference'):
            exit(1)


    def _check_executable(self):
        if not which('java'):
            sys.stderr.write('* Warning: Java not found. You may want to run "module load java", '
                             'or better ". /group/ngs/bin/bcbio-prod.sh"\n\n')
        if not which('perl'):
            sys.stderr.write('* Warning: Perl not found. You may want to run "module load perl", '
                             'or better ". /group/ngs/bin/bcbio-prod.sh"\n\n')
        if not self._get_tool_cmdline('vcfannotate',
                                      extra_warn='You may want to load BCBio '
                                                 'with ". /group/ngs/bin/bcbio-prod.sh"'):
            sys.stderr.write('* Warning: skipping annotation with bed tracks.\n')


    def _read_input(self):
        data = []

        if 'input' not in self.run_cnf:
            if 'file' not in self.run_cnf:
                self.log_exit('ERROR: Run config does not contain "input" section.')
            data.append({'vcf': realpath(self.run_cnf['file']),
                         'bam': self.run_cnf.get('bam'),
                         'bam_per_sample': None})
        else:
            for rec in self.run_cnf['input']:
                if 'vcf' not in rec:
                    self.log_exit('ERROR: Input section does not contain field "vcf".')
                data.append({'vcf': realpath(rec['vcf']),
                             'bam': rec.get('bam'),
                             'bam_per_sample': rec.get('bams')})
        return data


    def _check_file_for_rec(self, rec):
        inp_fpath = rec['vcf']
        if not verify_file(inp_fpath, 'Input file'):
            exit(1)
        bam_fpath = rec['bam']
        if bam_fpath:
            if not verify_file(bam_fpath, 'Bam file'):
                exit(1)
        bams = rec['bam_per_sample']
        if bams:
            bam_fpaths = [fpath for sample, fpath in bams.items()]
            if bam_fpaths:
                for bam_fpath in bam_fpaths:
                    if not verify_file(bam_fpath, 'Bam file'):
                        exit(1)
        return inp_fpath, bams, bam_fpath


    def _copy_to_final_dir(self, inp_fpath):
        fname = basename(inp_fpath)

        if self.output_dir != realpath(dirname(inp_fpath)):
            if self.run_cnf.get('save_intermediate'):
                new_fname = fname
            else:
                new_fname = add_suffix(fname, 'anno')
            new_fpath = join(self.output_dir, new_fname)

            if os.path.exists(new_fpath):
                os.remove(new_fpath)
            shutil.copyfile(inp_fpath, new_fpath)
            return new_fpath
        else:
            if not self.run_cnf.get('save_intermediate'):
                new_fpath = join(self.output_dir, add_suffix(fname, 'anno'))
                shutil.copyfile(inp_fpath, new_fpath)
                return new_fpath
            else:
                return inp_fpath


    def _read_samples(self, vcf_fpath):
        basic_fields = next(l.strip()[1:].split() for l in open(vcf_fpath)
                            if l.strip().startswith('#CHROM'))
        return basic_fields[9:]


    def _split_samples(self, vcf_fpath, bams):
        samples = self._read_samples(vcf_fpath)

        for sample in (bams.keys() if bams else []):
            if sample not in samples:
                self.log_exit('ERROR: sample ' + sample + ' is not in VCF ' + vcf_fpath + '\n'
                              'Available samples: ' + ', '.join(samples))
        for sample in samples:
            bam_fpath = bams.get(sample) if bams else None
            new_vcf = self.extract_sample(vcf_fpath, sample)
            self.samples[sample] = {'vcf': new_vcf, 'bam': bam_fpath}


    def __init__(self, system_config_path, run_config_path):
        self.sys_cnf = load(open(system_config_path), Loader=Loader)
        self.run_cnf = load(open(run_config_path), Loader=Loader)
        self.log = None
        self.samples = dict()

        if not 'resources' in self.sys_cnf:
            self.log_exit('"resources" section in system config required.')

        self._set_up_dir()

        if 'log' in self.run_cnf:
            self._set_up_log()
            shutil.copy(run_config_path, self.output_dir)

        self.log_print('Loaded system config ' + system_config_path)
        self.log_print('Loaded run config ' + run_config_path)
        self.log_print('')
        self.log_print('Writing into ' + self.output_dir)
        if self.log:
            self.log_print('Logging to ' + self.log)
        self.log_print('')

        self._check_executable()

        for rec in self._read_input():
            inp_fpath, bams, bam_fpath = self._check_file_for_rec(rec)

            inp_fpath = rec['vcf'] = self._copy_to_final_dir(inp_fpath)

            # Split by samples
            if bams or self.run_cnf.get('split_samples'):
                self._split_samples(inp_fpath, bams)
            else:
                self.samples[None] = {'vcf': inp_fpath, 'bam': bam_fpath}


    def annotate(self):
        if len(self.samples) == 1:
            name, sample = self.samples.items()[0]
            self.annotate_one(sample, name)
        else:
            try:
                from joblib import Parallel
            except ImportError:
                self.log_print('Joblib not found. You may want samples to be processed '
                               'in parallel, in this case, make sure python joblib intalled. '
                               '(pip install joblib).')
                for name, sample in self.samples:
                    self.annotate_one(sample, name)
            else:
                Parallel(n_jobs=len(self.samples))(
                    self.annotate_one(sample, name)
                    for name, sample in self.samples)


    def annotate_one(self, sample, sample_name):
        if sample_name:
            self.log_print('')
            self.log_print('')
            self.log_print('*' * 70)
            msg = '*' * 3 + ' Sample ' + sample_name + ' '
            self.log_print(msg + ('*' * (70 - len(msg)) if len(msg) < 70 else ''))
            self.log_print('VCF: ' + sample['vcf'])
            if sample.get('bam'):
                self.log_print('BAM: ' + sample['bam'])

        vcf_fpath = sample['vcf']
        # sample['fields'] = []

        remove_annotation('EFF', vcf_fpath)

        if self.run_cnf.get('split_genotypes'):
            vcf_fpath = self.split_genotypes(vcf_fpath)

        if self.run_cnf.get('ensemble'):
            vcf_fpath = self.process_ensemble(vcf_fpath)

        vcf_fpath = self.gatk(vcf_fpath, sample['bam'])

        if 'vcfs' in self.run_cnf:
            for dbname, conf in self.run_cnf['vcfs'].items():
                vcf_fpath = self.snpsift_annotate(dbname, conf, vcf_fpath)

        vcf_fpath = self.snpsift_db_nsfp(vcf_fpath)

        vcf_fpath = self.snpeff(vcf_fpath)

        if self.run_cnf.get('tracks'):
            for track in self.run_cnf['tracks']:
                vcf_fpath = self.tracks(track, vcf_fpath)

        vcf_fpath = self.filter_fields(vcf_fpath)

        tsv_fpath = self.extract_fields(vcf_fpath, sample_name)

        manual_tsv_fields = self.run_cnf.get('tsv_fields')
        if manual_tsv_fields:
            field_map = dict((rec.keys()[0], rec.values()[0]) for rec in manual_tsv_fields)
            tsv_fpath = self.rename_fields(tsv_fpath, field_map)
            self.log_print('Saved final TSV file with nice names to ' + tsv_fpath)

        if self.run_cnf.get('save_intermediate'):
            corr_tsv_fpath = dots_to_empty_cells(tsv_fpath)
            self.log_print('')
            self.log_print('TSV file with dots saved to ' + corr_tsv_fpath)
            self.log_print('View with the commandline:')
            self.log_print('    column -t ' + corr_tsv_fpath + ' | less -S')

        self.log_print('')
        self.log_print('Final VCF in ' + vcf_fpath)
        if self.log:
            print('Log in ' + self.log)


    def extract_sample(self, input_fpath, sample):
        self.introduce_step('Separating out sample ' + sample)

        executable = self._get_java_tool_cmdline('gatk')
        ref_fpath = self.run_cnf['reference']
        output_fpath = add_suffix(input_fpath, sample)
        with file_transaction(output_fpath) as tx_out_fpath:
            cmd = '{executable} -nt 30 -R {ref_fpath} -T SelectVariants ' \
                  '--variant {input_fpath} -o {output_fpath} -sn {sample}'.format(**locals())
            self._call_and_rename(cmd, input_fpath, suffix=sample,
                                  result_to_stdout=False, rename=False)
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


    def _call_and_rename(self, cmdline, input_fpath, suffix,
                         result_to_stdout=True, to_remove=None, rename=True):
        to_remove = to_remove or []

        output_fpath = add_suffix(input_fpath, suffix)

        if self.run_cnf.get('save_intermediate') and self.run_cnf.get('reuse_intermediate'):
            if file_exists(output_fpath):
                self.log_print(output_fpath + ' exists, reusing')
                return output_fpath

        self.log_print(cmdline)

        err_fpath = os.path.join(self.output_dir, 'annotate_py_err.tmp')
        to_remove.append(err_fpath)
        if err_fpath and isfile(err_fpath):
            os.remove(err_fpath)

        with file_transaction(output_fpath) as tx_out_fpath:
            if not result_to_stdout:
                cmdline += ' -o ' + tx_out_fpath

            if self.run_cnf.get('verbose', True):
                proc = subprocess.Popen(
                    cmdline.split(),
                    stdout=open(tx_out_fpath, 'w') if result_to_stdout else subprocess.PIPE,
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
                    stdout=open(tx_out_fpath, 'w') if result_to_stdout else open(err_fpath, 'a'),
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

        if rename and not self.run_cnf.get('save_intermediate'):
            os.rename(output_fpath, input_fpath)
            self.log_print('Saved to ' + input_fpath)
            return input_fpath
        else:
            self.log_print('Saved to ' + output_fpath)
            return output_fpath


    def _get_java_tool_cmdline(self, name):
        cmdline_template = self._get_script_cmdline_template('java', name)
        jvm_opts = self.sys_cnf['resources'][name].get('jvm_opts', [])
        return cmdline_template % (' '.join(jvm_opts) + ' -jar')


    def _get_script_cmdline_template(self, executable, script_name):
        if not which(executable):
            exit(executable + ' executable required, maybe you need '
                 'to run "module load ' + executable + '"?')
        if 'resources' not in self.sys_cnf:
            self.log_exit('System config yaml must contain resources section with '
                          + script_name + ' path.')
        if script_name not in self.sys_cnf['resources']:
            self.log_exit('System config resources section must contain '
                          + script_name + ' info (with a path to the tool).')
        tool_config = self.sys_cnf['resources'][script_name]
        if 'path' not in tool_config:
            self.log_exit(script_name + ' section in the system config must contain a path to the tool.')
        tool_path = tool_config['path']
        if not verify_file(tool_path, script_name):
            exit(1)
        return executable + ' %s ' + tool_path


    def _get_tool_cmdline(self, tool_name, extra_warn=''):
        tool_path = which(tool_name) or None

        if not 'resources' in self.sys_cnf \
                or tool_name not in self.sys_cnf['resources'] \
                or 'path' not in self.sys_cnf['resources'][tool_name]:
            if tool_path:
                return tool_path
            else:
                self.log_err(tool_name + ' executable was not found. '
                             'You can either specify path in the system config, or load into your '
                             'PATH environment variable.')
                if extra_warn:
                    self.log_err(extra_warn)
                return None

        tool_path = self.sys_cnf['resources'][tool_name]['path']
        if verify_file(tool_path, tool_name):
            return tool_path
        else:
            self.log_err(tool_path + ' for ' + tool_name +
                         ' does not exist or is not a file.')
            return None


    def snpsift_annotate(self, dbname, conf, input_fpath):
        self.introduce_step('Annotate with ' + dbname)

        executable = self._get_java_tool_cmdline('snpsift')
        db_path = conf.get('path')
        annotations = conf.get('annotations', [])
        # self.all_fields.extend(annotations)
        anno_line = ('-info ' + ','.join(annotations)) if annotations else ''
        cmdline = '{executable} annotate -v {anno_line} {db_path} {input_fpath}'.format(**locals())
        output_fpath = self._call_and_rename(cmdline, input_fpath, dbname, result_to_stdout=True)

        def proc_line(line):
            if not line.startswith('#'):
                line = line.replace(' ', '_')
                assert ' ' not in line
            return line
        output_fpath = self.iterate_file(output_fpath, proc_line, 'spaces')
        return output_fpath


    def snpsift_db_nsfp(self, input_fpath):
        if 'db_nsfp' not in self.run_cnf:
            return input_fpath

        self.introduce_step('DB SNFP')

        executable = self._get_java_tool_cmdline('snpsift')

        db_path = self.run_cnf['db_nsfp'].get('path')
        if not db_path:
            self.log_exit('Please, provide a path to db nsfp file in run_config.')

        annotations = self.run_cnf['db_nsfp'].get('annotations', [])
        # self.all_fields.extend(['dbNSFP_' + ann for ann in annotations])
        ann_line = ('-f ' + ','.join(annotations)) if annotations else ''

        cmdline = '{executable} dbnsfp {ann_line} -v {db_path} {input_fpath}'.format(**locals())
        return self._call_and_rename(cmdline, input_fpath, 'db_nsfp', result_to_stdout=True)


    def snpeff(self, input_fpath):
        if 'snpeff' not in self.run_cnf:
            return input_fpath

        self.introduce_step('SnpEff')

        # self.all_fields.extend([
        #     "EFF[*].EFFECT", "EFF[*].IMPACT", "EFF[*].FUNCLASS", "EFF[*].CODON",
        #     "EFF[*].AA", "EFF[*].AA_LEN", "EFF[*].GENE", "EFF[*].CODING",
        #     "EFF[*].TRID", "EFF[*].RANK"])

        executable = self._get_java_tool_cmdline('snpeff')
        ref_name = self.run_cnf['genome_build']
        db_path = self.run_cnf['snpeff'].get('path')
        assert db_path, 'Please, provide a path to db nsfp file in run_config.'

        cmdline = ('{executable} eff -dataDir {db_path} -noStats -noLog -1 '
                   '-i vcf -o vcf {ref_name} {input_fpath}').format(**locals())

        if self.run_cnf['snpeff'].get('clinical_reporting') or \
                self.run_cnf['snpeff'].get('canonical'):
            cmdline += ' -canon -hgvs '

        if self.run_cnf['snpeff'].get('cancer'):
            cmdline += ' -cancer '

        return self._call_and_rename(cmdline, input_fpath, 'snpEff', result_to_stdout=True)


    def tracks(self, track_path, input_fpath):
        field_name = splitext(basename(track_path))[0]

        self.introduce_step('Intersecting with ' + field_name)

        toolpath = self._get_tool_cmdline('vcfannotate')
        if not toolpath:
            self.log_err('WARNING: Skipping annotation with tracks: vcfannotate '
                         'executable not found, you probably need to '
                         'run the commandline:  . /group/ngs/bin/bcbio-prod.sh"')
            return

        # self.all_fields.append(field_name)

        cmdline = 'vcfannotate -b {track_path} -k {field_name} {input_fpath}'.format(**locals())

        out_fpath = self._call_and_rename(cmdline, input_fpath, field_name, result_to_stdout=True)

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
                    info_line = ';'.join('='.join(pair) if len(pair) == 2 else pair[0] for pair in info_pairs)
                    fields = fields[:7] + [info_line] + fields[8:]
                    return '\t'.join(fields)
            return line
        return self.iterate_file(out_fpath, proc_line)


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
        if 'gatk' not in self.run_cnf:
            return input_fpath

        self.introduce_step('GATK')

        executable = self._get_java_tool_cmdline('gatk')

        output_fpath = add_suffix(input_fpath, 'gatk')

        ref_fpath = self.run_cnf['reference']

        cmdline = ('{executable} -nt 20 -R {ref_fpath} -T VariantAnnotator'
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
        annotations = self.run_cnf['gatk'].get('annotations', [])

        # self.all_fields.extend(GATK_ANNOS_DICT.get(ann) for ann in annotations)
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

        return self._call_and_rename(cmdline, input_fpath, 'gatk',
                                     result_to_stdout=False,
                                     to_remove=[output_fpath + '.idx',
                                                input_fpath + '.idx'])


    def rename_fields(self, tsv_fpath, field_map):
        self.introduce_step('Renaming fields.')

        with open(tsv_fpath) as f:
            first_line = f.readline()[1:]
        fields = first_line.split()
        new_fields = [field_map.get(f, f) for f in fields]
        new_first_line = '\t'.join(new_fields)

        out_fpath = add_suffix(tsv_fpath, 'renamed')
        with file_transaction(out_fpath) as tx_out_fpath:
            with open(tx_out_fpath, 'w') as out:
                out.write(new_first_line + '\n')
                with open(tsv_fpath) as f:
                    for i, l in enumerate(f):
                        if i >= 1:
                            out.write(l)

        if self.run_cnf.get('save_intermediate'):
            return out_fpath
        else:
            os.rename(out_fpath, tsv_fpath)
            return tsv_fpath


    def extract_fields(self, vcf_fpath, sample_name=None):
        name, _ = splitext_plus(vcf_fpath)
        tsv_fpath = name + '.tsv'

        if self.run_cnf.get('save_intermediate') and self.run_cnf.get('reuse_intermediate'):
            if file_exists(tsv_fpath):
                self.log_print(tsv_fpath + ' exists, reusing')
                return tsv_fpath

        self.introduce_step('Extracting fields')

        all_format_fields = set()
        # Split FORMAT field and sample fields
        if sample_name:
            def proc_line(l):
                if l.startswith('#'):
                    return l
                vals = l.strip().split('\t')
                if len(vals) <= 9:
                    return l
                info = vals[7]
                format_fields = vals[8].split(':')
                sample_fields = vals[9].split(':')
                for f, s in zip(format_fields, sample_fields):
                    if f not in ['DP', 'MQ']:
                        if f == 'GT':
                            s = '"' + s + '"'
                        f = 'gt_' + f
                        info += ';' + f + '=' + s
                        all_format_fields.add(f)
                l = '\t'.join(vals[:7] + [info])
                return l
            tmp_vcf = self.iterate_file(vcf_fpath, proc_line)
        else:
            tmp_vcf = None

        manual_tsv_fields = self.run_cnf.get('tsv_fields')
        if manual_tsv_fields:
            fields = [rec.keys()[0] for rec in manual_tsv_fields]
        # else:
            # first_line = next(l.strip()[1:].split() for l in open(vcf_fpath) if l.strip().startswith('#CHROM'))
            # basic_fields = [f for f in first_line[:9] if f != 'INFO' and f != 'FORMAT' and f != sample_name]
            # manual_annots = filter(lambda f: f and f != 'ID', all_fields)
            # fields = (basic_fields + all_format_fields + manual_annots + self.run_cnf.get('additional_tsv_fields', []))
        else:
            return None

        anno_line = ' '.join(fields)
        snpsift_cmline = self._get_java_tool_cmdline('snpsift')
        vcfoneperline_cmline = self._get_script_cmdline_template('perl', 'vcfoneperline') % ''
        cmdline = vcfoneperline_cmline + ' | ' + snpsift_cmline + ' extractFields - ' + anno_line
        self.log_print(cmdline)
        with file_transaction(tsv_fpath) as tx_tsv_fpath:
            res = subprocess.call(cmdline,
                                  stdin=open(tmp_vcf or vcf_fpath),
                                  stdout=open(tx_tsv_fpath, 'w'), shell=True)
        if tmp_vcf:
            os.remove(tmp_vcf)

        self.log_print('')
        if res != 0:
            self.log_exit('Command returned status ' + str(res) +
                          ('. Log in ' + self.log if self.log else '.'))
        self.log_print('Saved TSV file to ' + tsv_fpath)
        return tsv_fpath


    def process_ensemble(self, input_fpath):
        self.introduce_step('Extracting dataset by filename, filtering ensemble reject line.')
        return self.iterate_file(input_fpath, lambda l: 'REJECT' not in l, 'pass')


    def iterate_file(self, input_fpath, proc_line_fun, suffix=None):
        output_fpath = add_suffix(input_fpath, suffix or 'tmp')

        if suffix and self.run_cnf.get('save_intermediate') \
                and self.run_cnf.get('reuse_intermediate'):
            if file_exists(output_fpath):
                self.log_print(output_fpath + ' exists, reusing')
                return output_fpath

        with file_transaction(output_fpath) as tx_fpath:
            with open(input_fpath) as vcf, open(tx_fpath, 'w') as out:
                for i, line in enumerate(vcf):
                    clean_line = line.strip()
                    if clean_line:
                        out.write(proc_line_fun(clean_line) + '\n')
                    else:
                        out.write(line)

        if suffix and self.run_cnf.get('save_intermediate'):
            return output_fpath
        else:
            os.rename(output_fpath, input_fpath)
            return input_fpath


    def split_genotypes(self, input_fpath):
        self.introduce_step('Splitting genotypes.')

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

        return self.iterate_file(input_fpath, proc_line, 'split_gt')


    def filter_fields(self, input_fpath):
        self.introduce_step('Filtering incorrect fields.')

        def proc_line(line):
            if line.startswith('#'):
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

        return self.iterate_file(input_fpath, proc_line, 'filt')


def dots_to_empty_cells(tsv_fpath):
    """Put dots instead of empty cells in order to view TSV with column -t
    """
    def proc_line(l):
        while '\t\t' in l:
            l = l.replace('\t\t', '\t.\t')
        return l
    return self.iterate_file(tsv_fpath, proc_line, 'tabs')


def remove_quotes(s):
    if s and s[0] == '"':
        s = s[1:]
    if s and s[-1] == '"':
        s = s[:-1]
    return s


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