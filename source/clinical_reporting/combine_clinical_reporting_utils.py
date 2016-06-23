import itertools
import json
import operator
from collections import defaultdict, OrderedDict
from os.path import dirname, basename, join

from source.clinical_reporting.clinical_parser import capitalize_keep_uppercase, get_sample_num, parse_vcf_record, get_record_from_vcf, \
    get_group_num, get_sample_info
from source.file_utils import open_gzipsafe
from source.reporting.reporting import Metric
from source.targetcov.bam_and_bed_utils import sambamba_depth
from source.targetcov.cov import get_mean_cov
from source.utils import gray
from source.variants import vcf_parser as vcf


def find_sample_parameters(row, mut_by_experiment, samples_data, cur_group_num=None):
    all_parameters = defaultdict(set)
    for e, m in mut_by_experiment.items():
        if cur_group_num and get_group_num(e.key) != cur_group_num:
            continue
        project_dirpath = dirname(dirname(e.sample.dirpath))
        sample_info = samples_data[project_dirpath][e.sample.name].data
        for parameter, value in sample_info.iteritems():
            all_parameters[parameter].add(value)
    for parameter, values in all_parameters.iteritems():
        row.add_record(parameter, ', '.join(sorted(values)))
    return row


def find_other_occurences(row, mut_by_experiment, cur_group_num, samples_data, parameters_info):
    num_by_samples = defaultdict(set)
    tooltips = []
    if cur_group_num:
        for e, m in mut_by_experiment.items():
            if get_group_num(e.key) == cur_group_num:
                continue
            sample_parameters = get_sample_info(e.sample.name, e.sample.dirpath, samples_data)
            sample_parameters = remove_parameters_to_combine(sample_parameters)
            short_parameters = [parameters_info.items()[i][1].prefixes[p.lower()] for i, p in enumerate(sample_parameters)]
            num_by_samples[tuple(short_parameters)].add(get_group_num(e.key))
            report_link = '<a href="' + basename(e.sample.clinical_html) + '" target="_blank">' + e.sample.name + '</a>'
            freq = Metric.format_value(m.freq, is_html=True, unit='%')
            tooltip = report_link + ':  ' + str(freq) + '  ' + str(m.depth) + '<br>'
            tooltips.append((e.sample.name, tooltip))
        tooltips = [tooltip[1] for tooltip in sorted(tooltips)]
        other_occurences = ', '.join([str(len(v)) + ''.join(k) for k, v in num_by_samples.iteritems()])
        other_occurences = add_tooltip(other_occurences, ''.join(tooltips))
        row.add_record('Other occurrences', other_occurences)
    return row


def get_depth_and_freq(e, mut, formatted_name, vcf_readers, filt_vcf_readers):
    depth = None
    freq = None
    tooltip = ''
    depth_str = None
    freq_record = None
    if e.rejected_mutations and (mut.gene.name, mut.pos) in e.rejected_mutations:
        rejected_mut = e.rejected_mutations[(mut.gene.name, mut.pos)]
        depth = rejected_mut.depth
        depth_str = gray(str(depth))
        tooltip = str(rejected_mut.reason)
        if rejected_mut.alt != mut.alt or (rejected_mut.aa_change and mut.aa_change != rejected_mut.aa_change):
            tooltip += '<br> Mutation: ' + str(rejected_mut.gene) + ' ' + str(rejected_mut.ref) + '>' + str(rejected_mut.alt) + \
                       ' ' + str(rejected_mut.aa_change)
        freq = rejected_mut.freq
        freq_record = add_tooltip(str(freq * 100), tooltip)
    else:
        record, vcf_reader = None, None
        if e in vcf_readers:
            record = get_record_from_vcf(vcf_readers[e], mut)
            vcf_reader = vcf_readers[e]
        if not record and e in filt_vcf_readers:
            record = get_record_from_vcf(filt_vcf_readers[e], mut)
            vcf_reader = filt_vcf_readers[e]
        if record:
            depth, freq, tooltip = parse_vcf_record(record, mut, e.sample.name, vcf_reader)
            if depth and freq:
                depth_str = gray(str(depth))
                freq_record = add_tooltip(str(freq * 100), tooltip)

    if not depth_str and not e.sample.bam:
        return None, None, None, None

    if not depth:
        depth = e.mutations_depth[mut.pos] if e.mutations_depth else 0
        depth_str = gray('%.0f' % depth)
        tooltip = 'No mutation'
    depth_record = add_tooltip(depth_str, tooltip)
    return depth_record, int(depth), freq_record, freq


def add_freq_depth_records(row, mut, mut_by_experiment, full_names, cur_group_num, samples_data, parameters_info,
                           vcf_readers, filt_vcf_readers):
    records_by_samples = defaultdict(list)
    if not mut_by_experiment.keys()[0].is_target2wgs_comparison:
        find_sample_parameters(row, mut_by_experiment, samples_data, cur_group_num)
        if cur_group_num:
            find_other_occurences(row, mut_by_experiment, cur_group_num, samples_data, parameters_info)
        else:
            used_names = set()
            for e, formatted_name in full_names.items():
                if e in mut_by_experiment:
                    used_names.add(formatted_name)
            row.add_record('Samples', len(used_names))
    for e, formatted_name in full_names.items():
        if e in mut_by_experiment:
            m = mut_by_experiment[e]
            if cur_group_num:
                row.add_record(formatted_name + ' Freq', m.freq if m else None, show_content=mut.is_canonical)
                row.add_record(formatted_name + ' Depth', m.depth if m else None, show_content=mut.is_canonical)
            else:
                sample_parameters = get_sample_info(e.sample.name, e.sample.dirpath, samples_data)
                records_by_samples[formatted_name].append((sample_parameters[-1], m.freq if m else None, m.depth if m else None, True))
        else:
            depth_record, depth, freq_record, freq = get_depth_and_freq(e, mut, formatted_name, vcf_readers, filt_vcf_readers)
            if cur_group_num:
                if depth_record:
                    row.add_record(formatted_name + ' Depth', depth_record, show_content=mut.is_canonical)
                    if freq_record:
                        row.add_record(formatted_name + ' Freq', freq_record, num=freq, show_content=mut.is_canonical, text_color='gray')
            else:
                sample_parameters = get_sample_info(e.sample.name, e.sample.dirpath, samples_data)
                records_by_samples[formatted_name].append((sample_parameters[-1], freq, depth, False))
    for formatted_name, records in records_by_samples.iteritems():
        tooltip = ''
        passed_records = [r for r in records if r[3]] or records
        max_index, freq = max(enumerate([r[1] for r in passed_records]), key=operator.itemgetter(1))
        depth = passed_records[max_index][2]
        is_passed = passed_records[max_index][3]
        for record in records:
            rec_freq = str(record[1] * 100) + '% ' if record[1] else '- '
            rec_depth = str(record[2]) if record[2] is not None else '-'
            rec_tooltip = record[0] + ': ' + rec_freq + rec_depth
            rec_is_passed = record[3]
            if not rec_is_passed:
                rec_tooltip = gray(rec_tooltip)
            tooltip += rec_tooltip + '<br>'

        if freq:
            freq_record = add_tooltip(str(freq * 100), tooltip)
            row.add_record(formatted_name + ' Freq', freq_record, num=freq, show_content=mut.is_canonical,
                           text_color='gray' if not is_passed else 'black')
        if depth:
            depth_record = add_tooltip(str(depth), tooltip)
            row.add_record(formatted_name + ' Depth', depth_record, num=depth, show_content=mut.is_canonical,
                           text_color='gray' if not is_passed else 'black')


def format_experiment_names(mutations_by_experiment, samples_data, cur_group_num=None):
    formatted_names = []
    short_names = OrderedDict()
    full_names = OrderedDict()
    used_full_names = defaultdict(int)

    for e in mutations_by_experiment.keys():
        formatted_name = ''
        if get_group_num(e.key) == cur_group_num or not cur_group_num:
            parameters = get_sample_info(e.sample.name, e.sample.dirpath, samples_data)
            if cur_group_num:
                formatted_name = ' '.join(parameters)
            else:
                formatted_name = ' '.join(remove_parameters_to_combine(parameters))
                formatted_name += '_' + get_sample_num(e.sample.name, e.sample.dirpath, samples_data)
                if formatted_name not in used_full_names:
                    used_full_names[formatted_name] = 0
        elif mutations_by_experiment.keys()[0].is_target2wgs_comparison:
            formatted_name = e.key
        formatted_names.append(formatted_name)

    if not mutations_by_experiment.keys()[0].is_target2wgs_comparison:
        prev_parameters = []
        for index in range(len(formatted_names)):
            formatted_name = formatted_names[index]
            if not formatted_name:
                continue
            if formatted_name in used_full_names and cur_group_num:
                used_full_names[formatted_name] += 1
                formatted_name += '_' + str(used_full_names[formatted_name])
            full_names[mutations_by_experiment.keys()[index]] = formatted_name
            parameters = formatted_name.split()
            formatted_name = parameters[-1]
            if len(parameters) > 1:
                for i, p in enumerate(parameters[-2: -len(parameters) - 1: -1]):
                    if p not in prev_parameters:
                        first_parameter = -(i + 2)
                        formatted_name = ' '.join(parameters[first_parameter:])
            prev_parameters = parameters
            short_names[mutations_by_experiment.keys()[index]] = formatted_name
    else:
        for index, e in enumerate(mutations_by_experiment.keys()):
            short_names[e] = formatted_names[index]
    return short_names, full_names


def group_for_venn_diagram(mutations_by_experiment, full_names, parameters_info, samples_data):
    samples_by_index = dict()
    set_labels = dict()
    parameters_sets = [remove_parameters_to_combine(parameter.values) for k, parameter in parameters_info.iteritems()]
    parameters_sets = [p_set for p_set in parameters_sets if p_set]
    base_groups = list(itertools.product(*parameters_sets))
    used_samples = set()
    set_index = 0
    for e in mutations_by_experiment.keys():
        if e not in full_names:
            continue
        sample_name = full_names[e]
        for index, g in enumerate(base_groups):
            if sample_name not in used_samples and \
                    all(capitalize_keep_uppercase(parameter) in get_sample_info(e.sample.name, e.sample.dirpath,
                                                                                samples_data) for parameter in g):
                used_samples.add(sample_name)
                samples_by_index[sample_name] = index
                set_labels[index] = ' '.join(g)
                set_index = max(set_index, index)
                set_index += 1
                break

    return samples_by_index, set_labels


def update_venn_diagram_data(sets, mut_by_experiment, samples_by_index, full_names):
    indexes = sorted(set([samples_by_index[full_names[e]] for e in mut_by_experiment.keys() if e in full_names]))
    for index in indexes:
        sets[index] += 1
    for len_set in range(2, len(indexes) + 1):
        for indexes_subset in itertools.combinations(indexes, len_set):
            sets[indexes_subset] += 1


def save_venn_diagram_data(venn_sets, set_labels):
    data = []
    for venn_set, size in venn_sets.iteritems():
        set_info = dict()
        set_info['size'] = size
        if isinstance(venn_set, tuple):
            set_info['sets'] = list(venn_set)
        else:
            set_info['sets'] = [venn_set]
        if isinstance(venn_set, int):
            set_info['label'] = set_labels[venn_set]
        data.append(set_info)
    return json.dumps(sorted(data, key=lambda x: x['sets']))


def get_vcf_readers(mutations_by_experiment, cur_group_num):
    vcf_readers, filt_vcf_readers = dict(), dict()
    for e, muts in mutations_by_experiment.items():
        if not cur_group_num or get_group_num(e.key) == cur_group_num:
            variant_caller = 'vardict' if 'vardict' in e.sample.variantcallers else 'vardict-java'
            if e.sample.vcf_by_callername.get(variant_caller):
                vcf_fpath = e.sample.vcf_by_callername.get(variant_caller)
                filt_vcf_fpath = e.sample.find_filt_vcf_by_callername(variant_caller)
                if vcf_fpath:
                    vcf_readers[e] = vcf.Reader(open_gzipsafe(vcf_fpath, 'r'))
                if filt_vcf_fpath:
                    filt_vcf_readers[e] = vcf.Reader(open_gzipsafe(filt_vcf_fpath, 'r'))
    return vcf_readers, filt_vcf_readers


def remove_parameters_to_combine(parameters):
    parameters_to_combine = ['WGS', 'AZ300', 'AZ50', 'Exome', 'WES']
    for parameter in parameters_to_combine:
        if parameter in parameters:
            parameters.remove(parameter)
    return parameters


def add_tooltip(text, tooltip):
    return ' <span class="my_hover"><div class="my_tooltip">' + tooltip + '</div> ' + text + ' </span>'

