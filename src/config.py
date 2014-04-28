import os
from os.path import basename, join, realpath
from collections import OrderedDict

from yaml import load
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

from src.quality_control import _check_quality_control_config
from src.utils import which, open_gzipsafe, file_exists, splitext_plus
from src.my_utils import get_tool_cmdline, err, critical, verify_file,\
    join_parent_conf, info, get_java_tool_cmdline, call, safe_mkdir, verify_dir


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


def _read_sample_names_from_vcf(vcf_fpath):
    basic_fields = next((l.strip()[1:].split() for l in open_gzipsafe(vcf_fpath)
                        if l.strip().startswith('#CHROM')), None)
    if not basic_fields:
        critical('Error: no VCF header in ' + vcf_fpath)
    if len(basic_fields) < 9:
        return []
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


def _check_system_resources(cnf):
    to_exit = False
    if not which('java'):
        err(cnf['log'], '\n* Warning: Java not found. You may want to run "module load java", '
            'or better ". /group/ngs/bin/bcbio-prod.sh"\n')
        to_exit = True

    if not which('perl'):
        err(cnf['log'], '\n* Warning: Perl not found. You may want to run "module load perl", '
            'or better ". /group/ngs/bin/bcbio-prod.sh"\n')
    if not get_tool_cmdline(cnf, 'vcfannotate',
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


def _read_samples_info(common_cnf):
    info('')
    info('Processing input details...')

    if 'details' not in common_cnf:
        critical('ERROR: Run config does not contain "details" section.')

    common_cnf['output_dir'] = common_cnf.get('output_dir', os.getcwd())
    common_cnf['filter_reject'] = common_cnf.get('filter_reject', False)
    common_cnf['split_samples'] = common_cnf.get('split_samples', False)
    common_cnf['log'] = None

    all_samples = OrderedDict()

    details = common_cnf['details']
    del common_cnf['details']

    for vcf_conf in details:
        if 'vcf' not in vcf_conf:
            critical('ERROR: A section in details does not contain field "vcf".')
        if not verify_file(vcf_conf['vcf'], 'Input file'):
            exit(1)

        join_parent_conf(vcf_conf, common_cnf)

        vcf_header_samples = _read_sample_names_from_vcf(vcf_conf['vcf'])

        # MULTIPLE SAMPELS
        if 'samples' in vcf_conf \
                or vcf_conf.get('split_samples') \
                or len(vcf_header_samples) == 0:
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

    return all_samples


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
            if bam and not verify_file(bam, 'Bam file'):
                exit()

    sample_cnfs = vcf_conf.get('samples') or OrderedDict()

    # compare input sample names to vcf header
    if sample_cnfs:
        for input_sample_name, sample_conf in sample_cnfs.items():
            if input_sample_name not in vcf_header_samples:
                critical('ERROR: sample ' + input_sample_name +
                         ' is not in VCF header ' + vcf_header_samples + '\n'
                         'Available samples: ' + ', '.join(vcf_header_samples))
    return sample_cnfs