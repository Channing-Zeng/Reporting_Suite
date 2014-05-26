#!/usr/bin/env python
import sys
import os
from os.path import basename, join, realpath, expanduser
from collections import OrderedDict

from source.bcbio_utils import open_gzipsafe, splitext_plus
from source.utils import critical, verify_file,\
    join_parent_conf, info, get_java_tool_cmdline, call, safe_mkdir, check_file_changed


def _set_up_work_dir(cnf):
    cnf['output_dir'] = expanduser(cnf['output_dir'])
    output_dirpath = realpath(cnf['output_dir'])
    safe_mkdir(output_dirpath, 'output_dir')
    work_dirpath = join(output_dirpath, 'work')
    safe_mkdir(work_dirpath, 'working directory')
    cnf['work_dir'] = work_dirpath


def read_samples_info_and_split(common_cnf, options, inputs):
    info('')
    info('Processing input details...')

    details = common_cnf.get('details')
    for key in inputs:
        if options.get(key):
            common_cnf[key] = expanduser(options[key])
            details = [common_cnf]

    all_samples = OrderedDict()

    for one_item_cnf in details:
        if 'var' in one_item_cnf:
            one_item_cnf['vcf'] = one_item_cnf['var']
            del one_item_cnf['var']
        if 'vcf' not in one_item_cnf:
            critical('ERROR: A section in details does not contain field "var".')
        one_item_cnf['vcf'] = expanduser(one_item_cnf['vcf'])
        if not verify_file(one_item_cnf['vcf'], 'Input file'):
            exit(1)

        join_parent_conf(one_item_cnf, common_cnf)

        work_vcf = join(one_item_cnf['work_dir'], basename(one_item_cnf['vcf']))
        check_file_changed(one_item_cnf, one_item_cnf['vcf'], work_vcf)
        if not one_item_cnf.get('reuse_intermediate'):
            with open_gzipsafe(one_item_cnf['vcf']) as inp, open(work_vcf, 'w') as out:
                out.write(inp.read())
        one_item_cnf['vcf'] = work_vcf

        vcf_header_samples = _read_sample_names_from_vcf(one_item_cnf['vcf'])

        # MULTIPLE SAMPELS
        if ('samples' in one_item_cnf or one_item_cnf.get('split_samples')) and len(vcf_header_samples) == 0:
            sample_cnfs = _verify_sample_info(one_item_cnf, vcf_header_samples)

            for header_sample_name in vcf_header_samples:
                if header_sample_name not in sample_cnfs:
                    sample_cnfs[header_sample_name] = one_item_cnf.copy()

                if header_sample_name in all_samples:
                    critical('ERROR: duplicated sample name: ' + header_sample_name)

                cnf = all_samples[header_sample_name] = sample_cnfs[header_sample_name]
                cnf['name'] = header_sample_name
                if cnf.get('keep_intermediate'):
                    cnf['log'] = join(cnf['work_dir'], cnf['name'] + '.log')

                cnf['vcf'] = extract_sample(cnf, one_item_cnf['vcf'], cnf['name'], cnf['work_dir'])
                info(cnf.get('log'), '')

        # SINGLE SAMPLE
        else:
            cnf = one_item_cnf

            if 'bam' in cnf:
                cnf['bam'] = expanduser(cnf['bam'])
                if not verify_file(cnf['bam']):
                    sys.exit(1)

            cnf['name'] = splitext_plus(basename(cnf['vcf']))[0]

            if cnf.get('keep_intermediate'):
                cnf['log'] = join(cnf['work_dir'], cnf['name'] + '.log')

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