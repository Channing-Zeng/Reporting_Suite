from __future__ import with_statement

import os
import shutil
import re
from os.path import join, abspath, dirname, isdir
import sys

from source.quast_reporting import json_saver
from source.logger import info, critical
from source.utils import verify_file


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


def save_total_report(results_dirpath, sample_names, report_base_name, report):
    json_fpath = json_saver.save_total_report(results_dirpath, report_base_name,
                                              sample_names, report)
    if not verify_file(json_fpath):
        sys.exit(1)

    html_fpath = append(results_dirpath, report_base_name + '.html', json_fpath, 'totalReport')
    info('  HTML version saved to ' + html_fpath)
    return html_fpath


def _init_html(results_dirpath, report_fname):
#    shutil.copy(template_fpath,     os.path.join(results_dirpath, report_fname))
    aux_dirpath = os.path.join(results_dirpath, aux_dirname)
    if not isdir(aux_dirpath):
        os.mkdir(aux_dirpath)

    for aux_f_relpath in aux_files:
        src_fpath = os.path.join(static_dirpath, aux_f_relpath)
        dst_fpath = os.path.join(aux_dirpath, aux_f_relpath)

        if not os.path.exists(os.path.dirname(dst_fpath)):
            os.makedirs(os.path.dirname(dst_fpath))

        if not os.path.exists(dst_fpath):
            shutil.copyfile(src_fpath, dst_fpath)

    with open(template_fpath) as template_file:
        html = template_file.read()
        html = html.replace("/" + static_dirname, aux_dirname)
        html = html.replace('{{ glossary }}', open(get_real_path('glossary.json')).read())

        html_fpath = os.path.join(results_dirpath, report_fname)
        if os.path.exists(html_fpath):
            os.remove(html_fpath)
        with open(html_fpath, 'w') as f_html:
            f_html.write(html)


def append(results_dirpath, report_fname, json_fpath, keyword):
    html_fpath = join(results_dirpath, report_fname)

    if not os.path.isfile(html_fpath):
        _init_html(results_dirpath, report_fname)

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


