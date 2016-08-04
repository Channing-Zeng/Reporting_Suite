#!/usr/bin/env python
# noinspection PyUnresolvedReferences
import bcbio_postproc

from collections import defaultdict, OrderedDict
import os
import sys
from os.path import join, abspath
from optparse import OptionParser

from source.bcbio.bcbio_structure import BCBioStructure, process_post_bcbio_args
from source.clinical_reporting.clinical_parser import capitalize_keep_uppercase, CombinedSampleInfo, Parameter
from source.clinical_reporting.combine_reports import run_combine_clinical_reports
from source.config import defaults
from source.logger import info, critical, err
from source.prepare_args_and_cnf import add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug, set_up_log
from source.file_utils import safe_mkdir, adjust_path
from source.utils import OrderedDefaultDict


def main():
    info(' '.join(sys.argv))
    info()
    description = 'This script makes clinical reports based on multiple bcbio projects.'

    parser = OptionParser(description=description)
    add_cnf_t_reuse_prjname_donemarker_workdir_genome_debug(parser)

    parser.add_option('--log-dir', dest='log_dir')
    parser.add_option('--email', dest='email', help='E-mail address to send notifications on errors and finished jobs.')
    parser.add_option('--metadata', dest='metadata_csv', help='CSV file with parameters of each sample.')
    parser.add_option('-o', dest='output_dir', help='Output directory for report combining.')

    cnf, bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths, tags, is_wgs_in_bcbio, is_rnaseq \
        = process_post_bcbio_args(parser)
    is_wgs = cnf.is_wgs = cnf.is_wgs or is_wgs_in_bcbio

    if not cnf.metadata_csv:
        critical('Provide the path to CSV file with information of each sample')
        critical('Usage: ' + __file__ + '  project_bcbio_path [project_bcbio_path] --metadata metadata_path [-o output_dir]')

    cnf.sample_names = []
    parameters_info, samples_data = parse_samples_metadata(cnf, cnf.metadata_csv)

    info()
    info('*' * 70)
    bcbio_structures = []
    for bcbio_project_dirpath, bcbio_cnf, final_dirpath in zip(
            bcbio_project_dirpaths, bcbio_cnfs, final_dirpaths):
        bs = BCBioStructure(cnf, bcbio_project_dirpath, bcbio_cnf, final_dirpath,
                            is_wgs=is_wgs, is_rnaseq=is_rnaseq)
        bcbio_structures.append(bs)

    if cnf.output_dir is None and cnf.project_name is None:
        critical('Either -o (output dir) or --project (project name) has to be specified')

    if not cnf.output_dir:
        cnf.output_dir = join(os.getcwd(), cnf.project_name)
    if not cnf.project_name:
        cnf.project_name = 'Combined_project'

    safe_mkdir(cnf.output_dir)

    cnf.log_dir = join(cnf.output_dir, 'log')
    info('log_dirpath: ' + cnf.log_dir)
    safe_mkdir(cnf.log_dir)
    set_up_log(cnf, 'combine_clin_reports', cnf.project_name, cnf.output_dir)

    cnf.work_dir = cnf.work_dir or adjust_path(join(cnf.output_dir, 'work'))
    safe_mkdir(cnf.work_dir)

    # shared_sample_names = set(s.name for bs in bcbio_structures for s in bs.samples)
    # if not shared_sample_names:
    #    sample_names = [bs.project_name + ': ' + ', '.join(s.name for s in bs.samples) for bs in bcbio_structures]
    #    critical('Not shared samples in projects.\n' + '\n'.join(sample_names))

    # info('Shared samples: ' + ', '.join(shared_sample_names))

    info('')
    info('*' * 70)
    run_combine_clinical_reports(cnf, bcbio_structures, parameters_info, samples_data)


def parse_samples_metadata(cnf, csv_fpath):
    sample_col = None
    project_col = None
    sample_num_col = None
    group_col = None
    additional_cols = []

    headers = []
    parameters = []

    parameters_info = OrderedDefaultDict(Parameter)
    samples_data = defaultdict(lambda : defaultdict(OrderedDict))
    with open(csv_fpath) as input_f:
        delim = None
        for i, l in enumerate(input_f):
            l = l.replace('\n', '')
            if not l:
                critical('Line ' + str(i) + ' is empty')
            if i == 0:
                if ',' in l:
                    delim = ','
                    info('Interpeting input as comma-separated')
                elif '\t' in l:
                    delim = '\t'
                    info('Interpeting input as tab-separated')
                else:
                    critical('Header is not separated by comma or tab: ' + l)
                headers = l.split(delim)
                sample_col = headers.index('sample') if 'sample' in headers else None
                project_col = headers.index('project_path') if 'project_path' in headers else None
                sample_num_col = headers.index('sample_num') if 'sample_num' in headers else None
                group_col = headers.index('report_num') if 'report_num' in headers else None
                if any(col is None for col in [sample_col, project_col, sample_num_col, sample_col, group_col]):
                    critical('Error: incorrect csv file format. File must have the following fields: sample,project_path,sample_num,report_num')

                additional_cols = [headers.index(col) for col in headers if headers.index(col) not in
                                   [sample_col, project_col, sample_num_col, group_col]]
                parameters = [capitalize_keep_uppercase(col) for col in headers[group_col + 1:]]
                continue
            fields = l.split(delim)
            if len(fields) < len(headers):
                critical('Error: len of line ' + str(i) + ' is ' + str(len(fields)) + ', which is less than the len of header (' + str(len(headers)) + ')')
            sample_name = fields[sample_col]
            cnf.sample_names.append(sample_name)
            project_path = abspath(fields[project_col])
            sample_num = fields[sample_num_col]
            group = fields[group_col]
            sample_info = CombinedSampleInfo(group=group, sample_num=sample_num)
            sample_data = sample_info.data
            for index, col in enumerate(additional_cols):
                parameter_name = parameters[index]
                parameter_value = fields[col]
                if parameter_value not in parameters_info[parameter_name].values:
                    parameters_info[parameter_name].values.add(parameter_value)
                sample_data[parameter_name] = parameter_value
            samples_data[project_path][sample_name] = sample_info
    for parameter_name, parameter in parameters_info.iteritems():
        values_list = list(parameter.values)
        values_list = [''] + sorted(values_list) + ['']
        def get_prefix_len(x):
            return len(os.path.commonprefix(x))
        for i, value in enumerate(values_list[:-1]):
            if i == 0:
                continue
            prefix_len = 1 + max(get_prefix_len(values_list[i-1:i+1]), get_prefix_len(values_list[i:i+2]))
            parameters_info[parameter_name].prefixes[value.lower()] = value[:prefix_len].capitalize()

    return parameters_info, samples_data


if __name__ == '__main__':
    main()
