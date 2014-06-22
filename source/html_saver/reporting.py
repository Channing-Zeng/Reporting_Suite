############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


# Here you can modify content and order of metrics in QUAST reports and names of metrcis as well
from collections import namedtuple
from os.path import basename
from source.logger import info, err


class Field:
    def __init__(self, name, unit='', thresholds=None):
        self.name = name
        self.unit = unit
        self.thresholds = thresholds

    def __str__(self, threshold=None):
        if threshold:
            return self.name % threshold
        else:
            return self.name

F = Field

# noinspection PyClassHasNoInit
class Fields:
####################################################################################
###########################  CONFIGURABLE PARAMETERS  ##############################
####################################################################################
    ### List of available fields for reports. Values (strings) should be unique! ###

    # Header
    NAME = F('SAMPLE')

    # Coverage
    READS = F('Reads')
    MAPPED_READS = F('Mapped reads')
    UNMAPPED_READS = F('Unmapped reads')
    PERCENT_MAPPED_READS = F('Percentage of mapped reads', '%')
    BASES_IN_TARGET = F('Bases in target', 'bp')
    COVERED_BASES_IN_TARGET = F('Covered bases in target', 'bp')
    PERCENT_TARGET_COVERED = F('Percentage of target covered by at least 1 read', '%')
    READS_MAPPED_ON_TARGET = F('Reads mapped on target')
    PERCENT_READS_MAPPED_ON_TARGET = F('Percentage of reads mapped on target', '%')
    READS_MAPPED_ON_PADDED_TARGET = F('Reads mapped on padded target')
    PERCENT_READS_MAPPED_ON_PADDED_TARGET = F('Percentage of reads mapped on padded target', '%')
    READ_BASES_MAPPED_ON_TARGET = F('Read bases mapped on target', 'bp')
    AVG_TARGET_COV_DEPTH = F('Average target coverage depth')
    STD_DEV_TARGET_COV_DEPTH = F('Std. dev. of target coverage depth')
    MAX_TARGET_COV_DEPTH = F('Maximum target coverage depth')
    PERCENT_TARGET_AROUND_MEAN = F('Percentage of target within 20% of mean depth', '%')

    PERCENT_TARGET_COVERED__AT_LEAST = F('Part of target covered at least by %dx', '%')

    order = [READS, MAPPED_READS, UNMAPPED_READS, PERCENT_MAPPED_READS, BASES_IN_TARGET,
             COVERED_BASES_IN_TARGET, PERCENT_TARGET_COVERED, READS_MAPPED_ON_TARGET,
             PERCENT_READS_MAPPED_ON_TARGET, READS_MAPPED_ON_PADDED_TARGET,
             PERCENT_READS_MAPPED_ON_PADDED_TARGET, READ_BASES_MAPPED_ON_TARGET,
             AVG_TARGET_COV_DEPTH, STD_DEV_TARGET_COV_DEPTH, MAX_TARGET_COV_DEPTH,
             PERCENT_TARGET_AROUND_MEAN]

####################################################################################
########################  END OF CONFIGURABLE PARAMETERS  ##########################
####################################################################################

    # noinspection PyClassHasNoInit
    class Quality:
        MORE_IS_BETTER = 'More is better'
        LESS_IS_BETTER = 'Less is better'
        EQUAL = 'Equal'

    quality_dict = {
        Quality.MORE_IS_BETTER: [READS, MAPPED_READS, PERCENT_MAPPED_READS, BASES_IN_TARGET,
                                 COVERED_BASES_IN_TARGET, PERCENT_TARGET_COVERED,
                                 READ_BASES_MAPPED_ON_TARGET, PERCENT_READS_MAPPED_ON_TARGET,
                                 READS_MAPPED_ON_PADDED_TARGET, PERCENT_READS_MAPPED_ON_PADDED_TARGET,
                                 READ_BASES_MAPPED_ON_TARGET, AVG_TARGET_COV_DEPTH, MAX_TARGET_COV_DEPTH,
                                 PERCENT_TARGET_AROUND_MEAN],

        Quality.LESS_IS_BETTER: [UNMAPPED_READS, STD_DEV_TARGET_COV_DEPTH],

        Quality.EQUAL: [],
    }

#################################################

import os
####################################################################################
# Reporting module (singleton) for QUAST
#
# See class Fields to available fields for report.
# Usage from QUAST modules:
#  from libs import reporting
#  report = reporting.get(fasta_filename)
#  report.add_field(reporting.Field.N50, n50)
#
####################################################################################

reports = {}  # basefilename -> Report
sample_fpaths = []  # for printing in appropriate order

#################################################

# def take_tuple_metric_apart(field):
#     metrics = []
#
#     if isinstance(field, tuple):
#         thresholds = map(int, ''.join(field[1]).split(','))
#         for i, feature in enumerate(thresholds):
#             metrics.append(field[0] % feature)
#     else:
#         metrics = [field]
#
#     return metrics


def get_quality(metric):
    for quality, metrics in Fields.quality_dict.iteritems():
        if metric in Fields.quality_dict[quality]:
            return quality
    return Fields.Quality.EQUAL


# Report for one filename, dict: field -> value
class Report(object):
    def __init__(self, name):
        self.d = {}
        self.add_field(Fields.NAME, name)

    def add_field(self, field, value):
        assert field in Fields.__dict__.itervalues(), 'Unknown field: %s' % field
        self.d[field] = value

    def append_field(self, field, value):
        assert field in Fields.__dict__.itervalues(), 'Unknown field: %s' % field
        self.d.setdefault(field, []).append(value)

    def get_field(self, field):
        assert field in Fields.__dict__.itervalues(), 'Unknown field: %s' % field
        return self.d.get(field, None)


def get(assembly_fpath):
    if assembly_fpath not in assembly_fpaths:
        assembly_fpaths.append(assembly_fpath)
    return reports.setdefault(assembly_fpath, Report(qutils.label_from_fpath(assembly_fpath)))


def delete(assembly_fpath):
    if assembly_fpath in assembly_fpaths:
        assembly_fpaths.remove(assembly_fpath)
    if assembly_fpath in reports.keys():
        reports.pop(assembly_fpath)


# ATTENTION! Contents numeric values, needed to be converted into strings
def table(order=Fields.order):
    if not isinstance(order[0], tuple):  # is not a groupped metrics order
        order = [('', order)]

    table = []

    def append_line(rows, field, are_multiple_tresholds=False, pattern=None, feature=None, i=None):
        quality = get_quality(field)
        values = []

        for assembly_fpath in assembly_fpaths:
            report = get(assembly_fpath)
            value = report.get_field(field)

            if are_multiple_tresholds:
                values.append(value[i] if (value and i < len(value)) else None)
            else:
                values.append(value)

        if filter(lambda v: v is not None, values):
            metric_name = field if (feature is None) else pattern % feature
            # ATTENTION! Contents numeric values, needed to be converted to strings.
            rows.append({
                'metricName': metric_name,
                'quality': quality,
                'values': values,
            })

    for group_name, metrics in order:
        rows = []
        table.append((group_name, rows))

        for field in metrics:
            if isinstance(field, tuple):  # TODO: rewrite it nicer
                for i, feature in enumerate(field[1]):
                    append_line(rows, field, True, field[0], feature, i)
            else:
                append_line(rows, field)

    if not isinstance(order[0], tuple):  # is not a groupped metrics order
        group_name, rows = table[0]
        return rows
    else:
        return table


def is_groupped_table(table):
    return isinstance(table[0], tuple)


def get_all_rows_out_of_table(table):
    all_rows = []
    if is_groupped_table(table):
        for group_name, rows in table:
            all_rows += rows
    else:
        all_rows = table

    return all_rows


def val_to_str(val):
    if val is None:
        return '-'
    else:
        return str(val)


def save_txt(fpath, table):
    all_rows = get_all_rows_out_of_table(table)

    # determine width of columns for nice spaces
    colwidths = [0] * (len(all_rows[0]['values']) + 1)
    for row in all_rows:
        for i, cell in enumerate([row['metricName']] + map(val_to_str, row['values'])):
            colwidths[i] = max(colwidths[i], len(cell))
            # output it

    txt_file = open(fpath, 'w')

    for row in all_rows:
        print >>txt_file, '  '.join('%-*s' % (colwidth, cell) for colwidth, cell
            in zip(colwidths, [row['metricName']] + map(val_to_str, row['values'])))

    txt_file.close()


def save_tsv(fpath, table):
    all_rows = get_all_rows_out_of_table(table)

    tsv_file = open(fpath, 'w')

    for row in all_rows:
        print >>tsv_file, '\t'.join([row['metricName']] + map(val_to_str, row['values']))

    tsv_file.close()


def parse_number(val):
    # Float?
    try:
        num = int(val)
    except ValueError:
        # Int?
        try:
            num = float(val)
        except ValueError:
            num = None

    return num


def get_num_from_table_value(val):
    if isinstance(val, int) or isinstance(val, float):
        num = val
    elif isinstance(val, basestring):
        num = parse_number(val)
    else:
        num = val

    return num


def save_tex(fpath, table, is_transposed=False):
    all_rows = get_all_rows_out_of_table(table)

    tex_file = open(fpath, 'w')
    # Header
    print >>tex_file, '\\documentclass[12pt,a4paper]{article}'
    print >>tex_file, '\\begin{document}'
    print >>tex_file, '\\begin{table}[ht]'
    print >>tex_file, '\\begin{center}'
    print >>tex_file, '\\caption{All statistics are based on contigs of size $\geq$ %d bp, unless otherwise noted ' % qconfig.min_contig + \
                      '(e.g., "\# contigs ($\geq$ 0 bp)" and "Total length ($\geq$ 0 bp)" include all contigs).}'

    rows_n = len(all_rows[0]['values'])
    print >>tex_file, '\\begin{tabular}{|l*{' + val_to_str(rows_n) + '}{|r}|}'
    print >>tex_file, '\\hline'

    # Body
    for row in all_rows:
        values = row['values']
        quality = row['quality'] if ('quality' in row) else Fields.Quality.EQUAL

        if is_transposed or quality not in [Fields.Quality.MORE_IS_BETTER, Fields.Quality.LESS_IS_BETTER]:
            cells = map(val_to_str, values)
        else:
            # Checking the first value, assuming the others are the same type and format
            num = get_num_from_table_value(values[0])
            if num is None:  # Not a number
                cells = map(val_to_str, values)
            else:
                nums = map(get_num_from_table_value, values)
                best = None
                if quality == Fields.Quality.MORE_IS_BETTER:
                    best = max(nums)
                if quality == Fields.Quality.LESS_IS_BETTER:
                    best = min(nums)

                if len([num for num in nums if num != best]) == 0:
                    cells = map(val_to_str, values)
                else:
                    cells = ['HIGHLIGHTEDSTART' + val_to_str(v) + 'HIGHLIGHTEDEND'
                             if get_num_from_table_value(v) == best
                             else val_to_str(v)
                             for v in values]

        row = ' & '.join([row['metricName']] + cells)
        # escape characters
        for esc_char in "\\ % $ # _ { } ~ ^".split():
            row = row.replace(esc_char, '\\' + esc_char)
        # more pretty '>=' and '<=', '>'
        row = row.replace('>=', '$\\geq$')
        row = row.replace('<=', '$\\leq$')
        row = row.replace('>', '$>$')
        # pretty indent
        if row.startswith(Fields.TAB):
            row = "\hspace{5mm}" + row.lstrip()
        # pretty highlight
        row = row.replace('HIGHLIGHTEDSTART', '{\\bf ')
        row = row.replace('HIGHLIGHTEDEND', '}')
        row += ' \\\\ \\hline'
        print >>tex_file, row

    # Footer
    print >>tex_file, '\\end{tabular}'
    print >>tex_file, '\\end{center}'
    print >>tex_file, '\\end{table}'
    print >>tex_file, '\\end{document}'
    tex_file.close()

    if os.path.basename(fpath) == 'report.tex':
        pass


def save(output_dirpath, report_name, order):
    # Where total report will be saved
    tab = table(order)

    info('  Creating total report...')
    report_txt_fpath = os.path.join(output_dirpath, report_name) + '.txt'
    report_tsv_fpath = os.path.join(output_dirpath, report_name) + '.tsv'
    report_tex_fpath = os.path.join(output_dirpath, report_name) + '.tex'
    save_txt(report_txt_fpath, tab)
    save_tsv(report_tsv_fpath, tab)
    save_tex(report_tex_fpath, tab)
    reports_fpaths = report_txt_fpath + ', ' + basename(report_tsv_fpath) + ', and ' + basename(report_tex_fpath)
    info('    saved to ' + reports_fpaths)

    return reports_fpaths


def save_total(output_dirpath, report_name):
    info('Summarizing...')
    return save(output_dirpath, report_name, Fields.order)
