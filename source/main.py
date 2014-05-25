#!/usr/bin/env python
import sys
import shutil
import os
from os.path import join, realpath, isdir, isfile, dirname, basename, join, realpath, expanduser
from optparse import OptionParser
from collections import OrderedDict

from source.utils import err, critical, verify_file,\
    join_parent_conf, info, get_java_tool_cmdline, call, safe_mkdir, verify_dir, verify_module

if verify_module('yaml'):
    from yaml import dump, load
    try:
        from yaml import CDumper as Dumper, CLoader as Loader
    except ImportError:
        from yaml import Dumper, Loader
else:
    critical('Cannot import module yaml.')

from source.bcbio_utils import which, open_gzipsafe, file_exists, splitext_plus


def common_main(name, opts):
    # options
    run_config_name = 'run_info_' + name + '.yaml'

    parser = OptionParser(
        usage='python ' + __file__ +
              ' [system_info.yaml] [' + run_config_name + '] ' +
              ' '.join(args[0] + example for args, example, _ in opts) +
              ' [--output_dir dir]\n'

              '    or python ' + __file__ +
              ' [system_info.yaml] ' + run_config_name)

    for args, _, kwargs in opts:
        parser.add_option(*args, **kwargs)
    parser.add_option('-o', '--output_dir', dest='output_dir', metavar='DIR')
    parser.add_option('--nt', '-t', dest='threads', help='number of threads to run GATK')
    (options, args) = parser.parse_args()

    # configs
    run_config_name = 'run_info_' + name + '.yaml'

    system_config_path = join(dirname(dirname(realpath(__file__))), 'system_info_rask.yaml')
    run_config_path = join(dirname(dirname(realpath(__file__))), run_config_name)
    if len(args) < 1:
        sys.stderr.write('Notice: using ' + run_config_name + ' as a default'
                         ' run configutation file.\n\n')
    else:
        run_config_path = args[0]

    if len(args) < 2:
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
        sys.stderr.write(system_config_path + ' does not end with .yaml,'
                                              ' maybe incorrect parameter?\n')
        to_exit = True
    if not run_config_path.endswith('.yaml'):
        sys.stderr.write(run_config_path + ' does not end with .yaml,'
                                           ' maybe incorrect parameter?\n')
        to_exit = True

    if to_exit:
        exit(1)

    _check_system_tools()

    sys_cnf = load(open(system_config_path), Loader=Loader)
    run_cnf = load(open(run_config_path), Loader=Loader)
    info('Loaded system config ' + system_config_path)
    info('Loaded run config ' + run_config_path)
    config = dict(run_cnf.items() + sys_cnf.items())

    if options.threads:
        if 'gatk' in config:
            gatk_opts = config['gatk'].get('options', [])
            new_opts = []
            for opt in gatk_opts:
                if opt.startswith('-nt '):
                    new_opts.append('-nt ' + options.threads)
                else:
                    new_opts.append(opt)
            config['gatk']['options'] = new_opts

    return config, options.__dict__


def read_samples_info_and_split(common_cnf, options):
    info('')
    info('Processing input details...')

    common_cnf['output_dir'] = common_cnf.get('output_dir', os.getcwd())
    common_cnf['filter_reject'] = common_cnf.get('filter_reject', False)
    common_cnf['split_samples'] = common_cnf.get('split_samples', False)
    common_cnf['log'] = None

    details = common_cnf.get('details')
    if not details and not options.get('vcf'):
        critical('ERROR: Provide input VCF by --var '
                 'or specify "details" section in run config.')
    if options.get('vcf'):
        vcf_conf = {'vcf': expanduser(options['vcf'])}
        if options.get('bam'):
            vcf_conf['bam'] = expanduser(options['bam'])
        if options.get('output_dir'):
            vcf_conf['output_dir'] = expanduser(options['output_dir'])
        details = [vcf_conf]

    all_samples = OrderedDict()

    for vcf_conf in details:
        if 'var' in vcf_conf:
            vcf_conf['vcf'] = vcf_conf['var']
            del vcf_conf['var']
        if 'vcf' not in vcf_conf:
            critical('ERROR: A section in details does not contain field "var".')
        vcf_conf['vcf'] = expanduser(vcf_conf['vcf'])
        if not verify_file(vcf_conf['vcf'], 'Input file'):
            exit(1)

        join_parent_conf(vcf_conf, common_cnf)

        vcf_header_samples = _read_sample_names_from_vcf(vcf_conf['vcf'])

        # MULTIPLE SAMPELS
        if ('samples' in vcf_conf or vcf_conf.get('split_samples')) and len(vcf_header_samples) == 0:
            sample_cnfs = _verify_sample_info(vcf_conf, vcf_header_samples)

            for header_sample_name in vcf_header_samples:
                if header_sample_name not in sample_cnfs:
                    sample_cnfs[header_sample_name] = vcf_conf.copy()

                if header_sample_name in all_samples:
                    critical('ERROR: duplicated sample name: ' + header_sample_name)

                cnf = all_samples[header_sample_name] = sample_cnfs[header_sample_name]
                cnf['name'] = header_sample_name
                _set_up_work_dir(cnf)

                cnf['vcf'] = extract_sample(cnf, vcf_conf['vcf'], cnf['name'], cnf['work_dir'])
                info(cnf.get('log'), '')

        # SINGLE SAMPLE
        else:
            if 'bam' in vcf_conf:
                vcf_conf['bam'] = expanduser(vcf_conf['bam'])
                if not verify_file(vcf_conf['bam']):
                    exit()

            cnf = vcf_conf

            cnf['name'] = splitext_plus(basename(cnf['vcf']))[0]

            _set_up_work_dir(cnf)

            work_vcf = join(cnf['work_dir'], cnf['name'] + '.vcf')
            if file_exists(work_vcf) and cnf.get('reuse_intermediate'):
                pass
            else:
                with open_gzipsafe(cnf['vcf']) as inp, open(work_vcf, 'w') as out:
                    out.write(inp.read())
            cnf['vcf'] = work_vcf

            all_samples[cnf['name']] = cnf

    if not all_samples:
        info('No samples.')
    else:
        info('Using samples: ' + ', '.join(all_samples) + '.')

    return all_samples


def _read_sample_names_from_vcf(vcf_fpath):
    basic_fields = next((l.strip()[1:].split() for l in open_gzipsafe(vcf_fpath)
                        if l.strip().startswith('#CHROM')), None)
    if not basic_fields:
        critical('Error: no VCF header in ' + vcf_fpath)
    if len(basic_fields) < 9:
        return []
    return basic_fields[9:]


def check_system_resources(cnf, required=list()):
    to_exit = False

    for program in required:
        if not which(program):
            resources = cnf.get('resources', None)

            if not resources:
                critical(cnf.get('log'), 'No "resources" section in system config.')

            data = resources.get(program)
            if data is None:
                err(program + ' is required. Specify path in system config or in your environment.')
                to_exit = True
            else:
                data['path'] = expanduser(data['path'])
                if not isdir(data['path']) and not file_exists(data['path']):
                    err(data['path'] + ' does not exist.')
                    to_exit = True
    if to_exit:
        exit()


def load_genome_resources(cnf, required=list()):
    if 'genome' not in cnf:
        critical('"genome" is not specified in run config.')
    if 'genomes' not in cnf:
        critical('"genomes" section is not specified in system config.')
    genome_name = cnf['genome']
    if genome_name not in cnf['genomes']:
        critical(genome_name + ' is not in "genomes section" of system config.')
    genome_cnf = cnf['genomes'][genome_name].copy()

    to_exit = False

    if 'seq' in required:
        required.remove('seq')

        if 'seq' not in genome_cnf:
            err('Please, provide path to the reference file (seq).')
            to_exit = True

        genome_cnf['seq'] = expanduser(genome_cnf['seq'])
        if not verify_file(genome_cnf['seq'], 'Reference seq'):
            to_exit = True

    if 'snpeff' in required:
        required.remove('snpeff')
        if 'snpeff' in genome_cnf:
            genome_cnf['snpeff'] = expanduser(genome_cnf['snpeff'])
            if not verify_dir(genome_cnf['snpeff'], 'snpeff'):
                to_exit = True

    for f in required:  #'dbsnp', 'cosmic', 'dbsnfp', '1000genomes':
        if f not in genome_cnf:
            err('Please, provide path to ' + f  + ' in system config genome section.')
            to_exit = True
        else:
            genome_cnf[f] = expanduser(genome_cnf[f])
            if not verify_file(genome_cnf[f], f):
                to_exit = True

    if to_exit:
        exit(1)

    cnf['genome'] = genome_cnf
    del cnf['genomes']
    genome_cnf['name'] = genome_name

    info('Loaded resources for ' + genome_cnf['name'])


def _check_system_tools():
    to_exit = False
    if not which('java'):
        err('\n* Warning: Java not found. You may want to run "module load java", '
            'or better ". /group/ngs/bin/bcbio-prod.sh"\n')
        to_exit = True

    if not which('perl'):
        err('\n* Warning: Perl not found. You may want to run "module load perl", '
            'or better ". /group/ngs/bin/bcbio-prod.sh"\n')

    # print ''
    # print 'In Waltham, run this as well:'
    # print '   export PATH=$PATH:/group/ngs/source/snpEff/snpEff3.5/scripts'
    # print '   export PERL5LIB=$PERL5LIB:/opt/az/local/bcbio-nextgen/'
    #       'stable/0.7.6/tooldir/lib/perl5/site_perl'
    if to_exit:
        exit()


def extract_sample(cnf, input_fpath, samplename, work_dir):
    info(cnf.get('log'), '-' * 70)
    info(cnf.get('log'), 'Separating out sample ' + samplename)
    info(cnf.get('log'), '-' * 70)

    executable = get_java_tool_cmdline(cnf, 'gatk')
    ref_fpath = cnf['genome']['seq']

    corr_samplename = ''.join([c if c.isalnum() else '_' for c in samplename])

    output_fname = splitext_plus(input_fpath)[0] + '.' + corr_samplename + '.vcf'
    output_fpath = join(work_dir, output_fname)

    cmd = '{executable} -nt 30 -R {ref_fpath} -T SelectVariants ' \
          '--variant {input_fpath} -sn {samplename} ' \
          '-o {output_fpath}'.format(**locals())
    call(cnf, cmd, None, output_fpath, stdout_to_outputfile=False,
         to_remove=[input_fpath + '.idx'])
    return output_fpath


def _set_up_work_dir(cnf):
    cnf['output_dir'] = expanduser(cnf['output_dir'])
    output_dirpath = realpath(cnf['output_dir'])
    safe_mkdir(output_dirpath, 'output_dir for ' + cnf['name'])
    work_dirpath = join(output_dirpath, 'work')
    safe_mkdir(work_dirpath, 'working directory')
    cnf['work_dir'] = work_dirpath
    if cnf.get('keep_intermediate'):
        cnf['log'] = join(work_dirpath, cnf['name'] + '.log')


def _verify_sample_info(vcf_conf, vcf_header_samples):
    # check bams if exist
    if 'samples' in vcf_conf:
        for header_sample_name, sample_conf in vcf_conf['samples'].items():
            join_parent_conf(sample_conf, vcf_conf)

            bam = sample_conf.get('bam')
            if bam and not verify_file(expanduser(bam), 'Bam file'):
                exit()
            sample_conf['bam'] = expanduser(bam)

    sample_cnfs = vcf_conf.get('samples') or OrderedDict()

    # compare input sample names to vcf header
    if sample_cnfs:
        for input_sample_name, sample_conf in sample_cnfs.items():
            if input_sample_name not in vcf_header_samples:
                critical('ERROR: sample ' + input_sample_name +
                         ' is not in VCF header ' + vcf_header_samples + '\n'
                         'Available samples: ' + ', '.join(vcf_header_samples))
    return sample_cnfs