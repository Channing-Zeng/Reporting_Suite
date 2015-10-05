# coding: utf-8

import os
import shutil
import re
import sys
import datetime
from os.path import join, abspath, dirname, isdir, splitext
from json import dumps, JSONEncoder
import traceback
from math import floor
import math

from source.bcbio_structure import VariantCaller, BCBioSample
from source.file_utils import verify_file, file_transaction
from source.file_utils import file_exists
from ext_modules.jsontemplate import jsontemplate
from source.logger import err
from source.utils import mean


def get_real_path(path_in_html_saver):
    return join(dirname(abspath(__file__)), path_in_html_saver)

scripts_inserted = False

template_fpath = get_real_path('template.html')
static_template_fpath = get_real_path('static_template.html')

static_dirname = 'static'
static_dirpath = get_real_path(static_dirname)

aux_dirname = 'html_aux'

css_files = [
    'bootstrap/bootstrap.min.css',
    'common.css',
    'report.css',
    'table_sorter/style.css'
]
js_files = [
    'jquery-1.8.2.min.js',
    # 'flot/jquery.flot.min.js',
    # 'flot/jquery.flot.dashes.js',
    'scripts/hsvToRgb.js',
    'scripts/utils.js',
    # 'scripts/build_total_report.js',
    # 'scripts/build_report.js',
    # 'flot/excanvas.min.js',
    # 'ie_html5.js',
    'bootstrap/bootstrap.min.js',
    'bootstrap/bootstrap-tooltip-vlad.js',
    # 'bootstrap/bootstrap-tooltip-5px-lower.min.js',
    'dragtable.js',
    'table_sorter/tsort.js',
]
image_files = [
    'table_sorter/arrow_asc.png',
    'table_sorter/arrow_desc.png',
    # 'img/draggable.png',
]


# def preprocess_report(report):


def _build_report(report, type_):
    report_html = ''
    report_html += _build_common_records(report.get_common_records())

    for section in report.metric_storage.sections:
        column_names = [m.name for m in section.metrics]
        # columnOrder = (recoverOrderFromCookies section.name) or report.order or [0...columnNames.length]
        column_order = range(len(column_names))  # TODO: read column order

        report_html += _build_total_report(report, section, column_order)

    return report_html


def _build_common_records(general_records):
    for rec in general_records:
        _calc_record_cell_contents(rec)

    table = '<table cellspacing="0" class="common_table" id="common_table">'
    for rec in general_records:
        table += '\n<tr><td>'
        if rec.html_fpath:
            table += '<a href="' + rec.html_fpath + '>' + rec.cell_contents + '</a>'
        else:
            table += '<span class="metric_name">' + __get_metric_name_html(rec.metric, use_full_name=True) + ': </span>' + rec.cell_contents
        table += '</td></tr>'
    table += '\n</table>\n'
    table += '<div class="space_30px"></div>'
    return table


def __get_metric_name_html(metric, use_full_name=False):
    if metric.short_name and not use_full_name:
        metric_name = metric.short_name
        description = metric.description or metric.name
        return '<a class="metric_name" rel="tooltip" title="' + description + '">' + metric_name + '</a>'
    else:
        return metric.name


def _build_total_report(report, section, column_order):
    html = ''

    if section.title:
        html += '\n<h3 class="table_name">' + section.title + '</h3>'
    _calc_cell_contents(report, section)

    table = '\n<table cellspacing="0" class="report_table tableSorter fix-align-char" id="report_table_' + section.name + '">'
    table += '\n<thead>\n<tr class="top_row_tr">'
    table += '\n<th class="top_left_td left_column_td" data-sortBy="numeric"><span>Sample</span></th>'

    for col_num in range(len(section.metrics)):
        pos = column_order[col_num]
        metric = section.metrics[pos]
        if not metric.is_hidden:
            if metric.numbers:
                sort_by = 'nosort' if metric.all_values_equal else 'numeric'
                direction = 'ascending' if metric.quality == 'Less is better' else 'descending'
                table += ('\n<th class="second_through_last_col_headers_td" data-sortBy=' + sort_by +
                          ' data-direction=' + direction + ' position=' + str(pos) + '>' +
                          '<span class="metricName">' + __get_metric_name_html(metric) + '</span></th>')
            elif metric.values:
                table += ('\n<th class="second_through_last_col_headers_td">' +
                          '<span class="metricName">' + __get_metric_name_html(metric) + '</span></th>')

    table += '\n</tr>\n</thead>\n<tbody>'

    i = 0
    if 'sample_reports' not in report.__dict__:
        sample_reports = [report]
    else:
        sample_reports = report.sample_reports
    sample_reports_length = len(sample_reports)

    for sample_report in sample_reports:
        line_caption = sample_report.display_name  # sample name
        max_sample_name_len = 50
        if len(line_caption) > max_sample_name_len:
            line_caption = '<span title="' + line_caption + '">' + line_caption[:max_sample_name_len] + '...</span>'

        second_row_tr = 'second_row_tr' if i == 0 else ''
        table += '\n<tr class="' + second_row_tr + '">\n<td class="left_column_td td" data-sortAs=' + str(sample_reports_length - i) + '>'
        if sample_reports_length == 1:
            table += '<span class="sample_name">' + line_caption + '</span>'
        else:
            if sample_report.html_fpath:
                if isinstance(sample_report.html_fpath, basestring):
                    table += '<a class="sample_name" href="' + sample_report.html_fpath + '">' + line_caption + '</a>'
                else:  # several links for one sample are possible multi-reports (e.g. TargQC)
                    if len(sample_report.html_fpath) == 0:
                        table += '<span class="sample_name">' + line_caption + '</span>'
                    else:
                        links = ', '.join('<a href="' + html_fpath + '">' + report_name + '</a>'
                                          for report_name, html_fpath in sample_report.html_fpath.items()
                                          if html_fpath)
                        table += '<span class="sample_name">' + line_caption + '(' + links + ')</span>'
            else:
                table += '<span class="sample_name">' + line_caption + '</span>'

        table += '\n</td>'

        for col_num in range(len(section.metrics)):
            pos = column_order[col_num]
            metric = section.metrics[pos]
            if '100x' in metric.name:
                pass
            if not metric.values:
                continue
            if metric.is_hidden:
                continue
            rec = sample_report.find_record(sample_report.records, metric.name)
            if not rec:
                table += "\n<td>-</td>"
                continue

            table += '\n<td metric="' + metric.name + '" style="background-color: ' + rec.color + '; color: ' + rec.text_color + \
                     '" quality="' + str(metric.quality) + '" class="td '
            if rec.num:
                table += ' number" number="' + str(rec.value) + '" data-sortAs="' + str(rec.value) + '">'
            else:
                table += '">'

            if rec.right_shift:
                padding_style = 'margin-left: ' + str(rec.right_shift) + 'px; margin-right: -' + str(rec.right_shift) + 'px;'
            else:
                padding_style = ''

            if rec.html_fpath:
                if isinstance(rec.html_fpath, basestring):
                    table += '<a href="' + rec.html_fpath + '">' + rec.cell_contents + '</a></td>'
                else:  # varQC -- several variant callers for one sample are possible
                    if len(rec.html_fpath) == 0:
                        rec.value = None
                        _calc_record_cell_contents(rec)
                        table += rec.cell_contents + '</td>'
                    else:
                        caller_links = ', '.join(
                            '<a href="' + html_fpath + '">' + caller + '</a>'
                            for caller, html_fpath in rec.html_fpath.items()
                            if html_fpath)
                        table += rec.cell_contents + '(' + caller_links + ')</td>'

            else:
                table += '<a style="' + str(padding_style) + '" ' + __get_meta_tag_contents(rec) + '>' + rec.cell_contents + '</a></td>'

        table += '\n</tr>'
        i += 1
    table += '\n</tbody>\n</table>\n'

    html += table
    return html


def __get_meta_tag_contents(rec):
    # meta = rec.meta
    #
    # if meta? and (a for a of meta).length != 0
    #     if typeof meta is 'string'
    #         return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\" title=\"#{meta}\""
    #
    #     else  # qc
    #         (k for own k of meta).length isnt 0
    #         meta_table = '<table class=\'qc_meta_table\'>\n<tr><td></td>'
    #         for novelty, values of meta when novelty isnt 'all'
    #             meta_table += "<td>#{novelty}</td>"
    #         meta_table += '</tr>\n'
    #
    #         for novelty, val_by_db of meta
    #             dbs = (db for db, val of val_by_db when db isnt 'average')
    #             dbs.push 'average'
    #             break
    #
    #         short_table = true
    #         for novelty, val_by_db of meta
    #             if not check_all_values_equal (val for db, val of val_by_db when db isnt 'average')
    #                 short_table = false
    #
    #         if short_table  # Values are the same for each database
    #             meta_table += '<tr><td></td>'
    #             for novelty, val_by_db of meta when novelty isnt 'all'
    #                 pretty_str = toPrettyString(val_by_db[dbs[0]], rec.metric.unit)
    #                 meta_table += "<td>#{pretty_str}</td>"
    #             meta_table += '</tr>\n'
    #         else
    #             for db in dbs
    #                 meta_table += "<tr><td>#{db}</td>"
    #                 for novelty, val_by_db of meta when novelty isnt 'all'
    #                     meta_table += "<td>#{toPrettyString(val_by_db[db], rec.metric.unit)}</td>"
    #                 meta_table += '</tr>\n'
    #
    #         meta_table += '</table>\n'
    #
    #         return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\" title=\"#{meta_table}\""
    # else
    return "class=\"meta_info_span tooltip-meta\" rel=\"tooltip\""


# Color heat map
BLUE_HUE = 240
BLUE_OUTER_BRT = 55
BLUE_INNER_BRT = 65

GREEN_HUE = 120
GREEN_OUTER_BRT = 50
GREEN_INNER_BRT = 60

RED_HUE = 0
RED_OUTER_BRT = 50
RED_INNER_BRT = 60

MIN_NORMAL_BRT = 80
MEDIAN_BRT = 100  # just white.


def hsl2rgb(h, s, l):
    r, g, b = None, None, None

    if s == 0:
        r = g = b = l  # achromatic
    else:
        q = l * (1 + s) if l < 0.5 else l + s - l * s
        p = 2 * l - q
        r = hue2rgb(p, q, h + 1./3)
        g = hue2rgb(p, q, h)
        b = hue2rgb(p, q, h - 1./3)

    return map(int, [round(r * 255), round(g * 255), round(b * 255)])

def hue2rgb(p, q, t):
    if t < 0: t += 1
    if t > 1: t -= 1
    if t < 1./6: return p + (q - p) * 6 * t
    if t < 1./2: return q
    if t < 2./3: return p + (q - p) * (2./3 - t) * 6
    return p


def get_color(hue, lightness):
    lightness = lightness or 92
    # lightness = Math.round (Math.pow hue - 75, 2) / 350 + 35
    rgb = hsl2rgb(float(hue) / 360, 0.8, float(lightness) / 100)
    hex_rgb = [hex(c)[2:] for c in rgb]
    return '#' + ''.join(hex_rgb)


def _calc_record_cell_contents(rec):
    rec.cell_contents = rec.format_html()
    if rec.value is not None and (isinstance(rec.value, int) or isinstance(rec.value, float)):
        rec.num = rec.value

    #TODO: intPartTextWidth
    return rec


def _calc_cell_contents(report, section):
    max_frac_widths_by_metric = dict()

    # First round: calculatings max/min integral/fractional widths (for decimal alingment) and max/min values (for heatmaps)
    if 'sample_reports' in report.__dict__:
        sample_reports = report.sample_reports
    else:
        sample_reports = [report]

    for sample_report in sample_reports:
        for rec in sample_report.records:
            if rec.metric.name in section.metrics_by_name:
                _calc_record_cell_contents(rec)

        for rec in sample_report.records:
            if rec.metric.name in section.metrics_by_name:
                if rec.metric.name not in max_frac_widths_by_metric or \
                                rec.frac_width > max_frac_widths_by_metric[rec.metric.name]:
                    max_frac_widths_by_metric[rec.metric.name] = rec.frac_width
                if rec.num:
                    rec.metric.numbers.append(rec.num)
                if rec.value is not None:
                    rec.metric.values.append(rec.value)
    # else:
    #     for rec in report.records:
    #         if rec.metric.name in section.metrics_by_name:
    #             if rec.num:
    #                 if not rec.metric.values:
    #                     rec.metric.values = []
    #                 rec.metric.values.append(rec.num)

    for metric in section.metrics:
        if metric.numbers:
            def _cmp(a, b):  # None is always less than anything
                if a and b:
                    return cmp(a, b)
                elif a:
                    return 1
                else:
                    return -1

            numbers = sorted(
                [v for v in metric.numbers],
                cmp=_cmp)
            l = len(numbers)

            metric.min = numbers[0]
            metric.max = numbers[l - 1]
            metric.all_values_equal = metric.min == metric.max
            metric.med = numbers[(l - 1) / 2] if l % 2 != 0 else mean([numbers[l / 2], numbers[(l / 2) - 1]])
            q1 = numbers[int(floor((l - 1) / 4))]
            q3 = numbers[int(floor((l - 1) * 3 / 4))]

            d = q3 - q1
            metric.low_outer_fence = q1 - 3   * d
            metric.low_inner_fence = q1 - 1.5 * d
            metric.top_inner_fence = q3 + 1.5 * d
            metric.top_outer_fence = q3 + 3   * d

    if 'sample_reports' in report.__dict__:
        # Second round: setting shift and color properties based on max/min widths and vals
        for sample_report in report.sample_reports:
            for rec in sample_report.records:
                if rec.metric and rec.metric.name in section.metrics_by_name:
                    # Padding based on frac width
                    if rec.frac_width:
                        rec.right_shift = max_frac_widths_by_metric[rec.metric.name] - rec.frac_width

                    metric = rec.metric

                    # For metrics where we know the "normal value" - we want to color everything above normal white,
                    #   and everything below - red, starting from normal, finishing with bottom
                    if metric.ok_threshold is not None:
                        if isinstance(metric.ok_threshold, int) or isinstance(metric.ok_threshold, float):
                            metric.med = metric.ok_threshold
                            if metric.bottom is not None:
                                metric.low_outer_fence = metric.bottom

                    # Color heatmap
                    if rec.num:
                        [top_hue, inner_top_brt, outer_top_brt] = [BLUE_HUE, BLUE_INNER_BRT, BLUE_OUTER_BRT]
                        [low_hue, inner_low_brt, outer_low_brt] = [RED_HUE, RED_INNER_BRT, RED_OUTER_BRT]

                        if metric.quality == 'Less is better':  # then swap colors
                            [top_hue, low_hue] = [low_hue, top_hue]
                            [inner_top_brt, inner_low_brt] = [inner_low_brt, inner_top_brt]
                            [outer_top_brt, outer_low_brt] = [outer_low_brt, outer_top_brt]

                        if metric.ok_threshold is not None:
                            if isinstance(metric.ok_threshold, int) or isinstance(metric.ok_threshold, float):
                                if rec.num >= metric.ok_threshold:
                                    continue  # white on blak

                                # rec_to_align_with = sample_report.find_record(sample_report.records, metric.threshold)
                                # if rec_to_align_with:
                                #     rec.text_color = lambda: rec_to_align_with.text_color()
                                #     continue

                        if not metric.all_values_equal:
                            rec.text_color = 'black'

                            # Low outliers
                            if rec.num < rec.metric.low_outer_fence:
                                rec.color = get_color(low_hue, outer_low_brt)
                                rec.text_color = 'white'

                            elif rec.num < rec.metric.low_inner_fence:
                                rec.color = get_color(low_hue, inner_low_brt)

                            # Normal values
                            elif rec.num < metric.med:
                                k = (MEDIAN_BRT - MIN_NORMAL_BRT) / (metric.med - rec.metric.low_inner_fence)
                                brt = round(MEDIAN_BRT - (metric.med - rec.num) * k)
                                rec.color = get_color(low_hue, brt)

                            # High outliers
                            elif rec.num > rec.metric.top_inner_fence:
                                rec.color = get_color(top_hue, inner_top_brt)

                            elif rec.num > rec.metric.top_outer_fence:
                                rec.color = get_color(top_hue, outer_top_brt)
                                rec.text_color = 'white'

                            elif rec.num > metric.med:
                                k = (MEDIAN_BRT - MIN_NORMAL_BRT) / (rec.metric.top_inner_fence - metric.med)
                                brt = round(MEDIAN_BRT - (rec.num - metric.med) * k)
                                rec.color = get_color(top_hue, brt)

            for rec in sample_report.records:
                if rec.metric and rec.metric.name in section.metrics_by_name:
                    metric = rec.metric

                    if metric.ok_threshold is not None:
                        if isinstance(metric.ok_threshold, basestring):
                            rec_to_align_with = sample_report.find_record(sample_report.records, metric.ok_threshold)
                            if rec_to_align_with:
                                rec.text_color = rec_to_align_with.text_color
                                rec.color = rec_to_align_with.color
    return report


def write_html_report(report, type_, html_fpath, caption=''):
    with open(template_fpath) as template_file:
        html = template_file.read()

    html = _insert_into_html(html, caption, 'caption')
    html = _insert_into_html(html, datetime.datetime.now().strftime('%d %B %Y, %A, %H:%M:%S'), 'report_date')

    report_html = _build_report(report, type_)
    html = _insert_into_html(html, report_html, 'report')

    plots_html = ''.join('<img src="' + plot + '"/>' for plot in report.plots)
    html = _insert_into_html(html, plots_html, 'plots')

    html = _embed_css_and_scripts(html)

    with open(html_fpath, 'w') as f:
        f.write(html)
    return html_fpath


def _copy_aux_files(results_dirpath):
    aux_dirpath = join(results_dirpath, aux_dirname)

    if isdir(aux_dirpath):
        shutil.rmtree(aux_dirpath)

    if not isdir(aux_dirpath):
        try:
            os.mkdir(aux_dirpath)
        except OSError:
            pass

    def copy_aux_file(fname):
        src_fpath = join(static_dirpath, fname)

        if file_exists(src_fpath):
            dst_fpath = join(aux_dirpath, fname)

            if not file_exists(dirname(dst_fpath)):
                try:
                    os.makedirs(dirname(dst_fpath))
                except OSError:
                    pass
            try:
                shutil.copyfile(src_fpath, dst_fpath)
            except OSError:
                pass

    for aux_f_relpath in js_files + css_files + image_files:
        if aux_f_relpath.endswith('.js'):
            for ext in ['.js', '.coffee', '.map']:
                copy_aux_file(splitext(aux_f_relpath)[0] + ext)

        elif aux_f_relpath.endswith('.css'):
            for ext in ['.css', '.sass']:
                copy_aux_file(splitext(aux_f_relpath)[0] + ext)

        else:
            copy_aux_file(aux_f_relpath)


def _embed_css_and_scripts(html):
    js_line_tmpl = '<script type="text/javascript" src="/static/{file_rel_path}"></script>'
    js_l_tag = '<script type="text/javascript" name="{}">'
    js_r_tag = '    </script>'

    css_line_tmpl = '<link rel="stylesheet" type="text/css" href="/static/{file_rel_path}" />'
    css_l_tag = '<style type="text/css" rel="stylesheet" name="{}">'
    css_r_tag = '    </style>'

    for line_tmpl, files, l_tag, r_tag in [
            (js_line_tmpl, js_files, js_l_tag, js_r_tag),
            (css_line_tmpl, css_files, css_l_tag, css_r_tag)]:
        for rel_fpath in files:
            with open(join(static_dirpath, join(*rel_fpath.split('/')))) as f:
                contents = f.read()
                contents = '\n'.join(' ' * 8 + l for l in contents.split('\n'))

                line = line_tmpl.format(file_rel_path=rel_fpath)
                l_tag_fmt = l_tag.format(rel_fpath)
                # line_ascii = line.encode('ascii', 'ignore')
                # html = html.decode('utf-8', 'ignore').encode('ascii')
                # try:
                #     html_ascii = html.decode('utf-8', 'ignore').encode('ascii')  #.encode('ascii', 'ignore')
                # except:
                #     pass
                try:
                    html = ''.join(i for i in html if ord(i) < 128)
                    contents = ''.join(i for i in contents if ord(i) < 128)
                    html = html.replace(line, l_tag_fmt + '\n' + contents + '\n' + r_tag)
                except:
                    err()
                    # print contents
                    err()
                    err(traceback.format_exc())
                    err()
                    err('Encoding problem embeding this file into html: ' + rel_fpath)
                    err()

                # try:
                #     html = html_ascii.encode('ascii').replace(line, l_tag_fmt + '\n' + contents + '\n' + r_tag)
                # except:
                #     pass
                    # err('line - cannot replace')  # + str(l_tag_fmt + '\n' + contents + '\n' + r_tag))

    return html


def _insert_into_html(html, text, keyword):
    # substituting template text with json
    # html_text = re.sub('{{ ' + keyword + ' }}', text, html_text)
    html = ''.join(i for i in html if ord(i) < 128)
    html = html.replace('{{ ' + keyword + ' }}', text)

    return html


def write_static_html_report(work_dir, data_dict, html_fpath, tmpl_fpath=None):
    tmpl_fpath = tmpl_fpath or static_template_fpath
    with open(tmpl_fpath) as f:
        html = f.read()

    html = jsontemplate.expand(html, data_dict)

    html = _embed_css_and_scripts(html)

    with file_transaction(work_dir, html_fpath) as tx:
        with open(tx, 'w') as f:
            f.write(html)

    return html_fpath
