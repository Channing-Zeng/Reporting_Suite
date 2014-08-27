from __future__ import with_statement

import os
import shutil
import re
from os.path import join, abspath, dirname, isdir, splitext
import sys
from source.file_utils import verify_file

from source.quast_reporting import json_saver
from source.file_utils import file_exists


def get_real_path(path_in_html_saver):
    return join(dirname(abspath(__file__)), path_in_html_saver)

scripts_inserted = False

template_fpath = get_real_path('template.html')

static_dirname = 'static'
static_dirpath = get_real_path(static_dirname)

aux_dirname = 'html_aux'
aux_files = [
    'jquery-1.8.2.min.js',
    # 'flot/jquery.flot.min.js',
    # 'flot/excanvas.min.js',
    # 'flot/jquery.flot.dashes.js',
    'scripts/build_total_report.js',
    # 'scripts/draw_cumulative_plot.js',
    # 'scripts/draw_nx_plot.js',
    # 'scripts/draw_gc_plot.js',
    'scripts/utils.js',
    'scripts/hsvToRgb.js',
    # 'scripts/draw_genes_plot.js',
    'scripts/build_report.js',
    'dragtable.js',
    'ie_html5.js',
    'img/draggable.png',
    'bootstrap/bootstrap-tooltip-5px-lower.min.js',
    'bootstrap/bootstrap.min.css',
    'bootstrap/bootstrap.min.js',
    'bootstrap/bootstrap-tooltip-vlad.js',
    'report.css',
    'common.css',
    'table_sorter/tsort.js',
    'table_sorter/style.css',
    'table_sorter/arrow_asc.png',
    'table_sorter/arrow_desc.png',
]


def write_html_reports(output_dirpath, work_dirpath, full_reports, report_base_name, caption):
    common_records = get_common_records(full_reports)
    json_fpath = json_saver.save_total_report(work_dirpath, report_base_name, full_reports, common_records)

    if not verify_file(json_fpath):
        sys.exit(1)

    html_fpath = init_html(output_dirpath, report_base_name + '.html', caption)
    append(html_fpath, json_fpath, 'totalReport')
    return html_fpath


def get_common_records(full_reports):
    if not isinstance(full_reports, list):
        full_reports = [full_reports]

    common_records = list()
    for full_report in full_reports:
        if full_report.sample_reports:
            sample_report = full_report.sample_reports[0]
            for record in sample_report.records:
                if record.metric.common:
                    common_records.append(record)
            for sample_report in full_report.sample_reports:
                sample_report.records = [record for record in sample_report.records if not record.metric.common]
    return common_records


def init_html(results_dirpath, report_fname, caption=''):
#    shutil.copy(template_fpath, os.path.join(results_dirpath, report_fname))
    aux_dirpath = join(results_dirpath, aux_dirname)
    if isdir(aux_dirpath):
        shutil.rmtree(aux_dirpath)
    os.mkdir(aux_dirpath)

    def copy_aux_file(fname):
        src_fpath = join(static_dirpath, fname)

        if file_exists(src_fpath):
            dst_fpath = join(aux_dirpath, fname)

            if not file_exists(dirname(dst_fpath)):
                os.makedirs(dirname(dst_fpath))

            shutil.copyfile(src_fpath, dst_fpath)

    for aux_f_relpath in aux_files:
        if aux_f_relpath.endswith('.js'):
            for ext in ['.js', '.coffee', '.map']:
                copy_aux_file(splitext(aux_f_relpath)[0] + ext)

        elif aux_f_relpath.endswith('.css'):
            for ext in ['.css', '.sass']:
                copy_aux_file(splitext(aux_f_relpath)[0] + ext)

        else:
            copy_aux_file(aux_f_relpath)

    with open(template_fpath) as template_file:
        html = template_file.read()
        html = html.replace("/" + static_dirname, aux_dirname)
        html = html.replace('{{ caption }}', caption)

        html_fpath = os.path.join(results_dirpath, report_fname)
        if os.path.exists(html_fpath):
            os.remove(html_fpath)
        with open(html_fpath, 'w') as f_html:
            f_html.write(html)

    return html_fpath


def append(html_fpath, json_fpath, keyword):
    # reading JSON file
    with open(json_fpath) as f_json:
        json_text = f_json.read()
    os.remove(json_fpath)

    # reading html template file
    with open(html_fpath) as f_html:
        html_text = f_html.read()

    # substituting template text with json
    html_text = re.sub('{{ ' + keyword + ' }}', json_text, html_text)

    # writing substituted html to final file
    with open(html_fpath, 'w') as f_html:
        f_html.write(html_text)

    return html_fpath
