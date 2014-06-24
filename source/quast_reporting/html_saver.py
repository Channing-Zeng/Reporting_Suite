from __future__ import with_statement

import os
import shutil
import re
from os.path import join, abspath, dirname, isdir
import sys
from source.file_utils import verify_file

from source.quast_reporting import json_saver
from source.utils_from_bcbio import file_exists


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
]


def write_html_report(output_dirpath, work_dirpath, report,
                      sample_names, report_base_name, caption):

    json_fpath = json_saver.save_total_report(
        work_dirpath, report_base_name, sample_names, report)

    if not verify_file(json_fpath):
        sys.exit(1)

    html_fpath = init_html(output_dirpath, report_base_name + '.html', caption)
    append(html_fpath, json_fpath, 'totalReport')
    return html_fpath


def init_html(results_dirpath, report_fname, caption=''):
#    shutil.copy(template_fpath, os.path.join(results_dirpath, report_fname))
    aux_dirpath = join(results_dirpath, aux_dirname)
    if isdir(aux_dirpath):
        shutil.rmtree(aux_dirpath)
    os.mkdir(aux_dirpath)

    for aux_f_relpath in aux_files:
        src_fpath = join(static_dirpath, aux_f_relpath)
        dst_fpath = join(aux_dirpath, aux_f_relpath)

        if not file_exists(dirname(dst_fpath)):
            os.makedirs(dirname(dst_fpath))

        if not file_exists(dst_fpath):
            shutil.copyfile(src_fpath, dst_fpath)

    with open(template_fpath) as template_file:
        html = template_file.read()
        html = html.replace("/" + static_dirname, aux_dirname)
        html = html.replace('{{ caption }}', caption)
        with open(get_real_path('glossary.json')) as glos_f:
            html = html.replace('{{ glossary }}', glos_f.read())

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
