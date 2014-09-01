import os
import shutil
import re
import sys
import datetime
from os.path import join, abspath, dirname, isdir, splitext
from json import dumps, JSONEncoder

from source.bcbio_structure import VariantCaller, Sample
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


def write_html_report(json, output_dirpath, report_base_name, caption):
    html_fpath = _init_html(output_dirpath, report_base_name + '.html', caption)
    _append(html_fpath, json, 'totalReport')
    return html_fpath


def _init_html(results_dirpath, report_fname, caption=''):
    # shutil.copy(template_fpath, os.path.join(results_dirpath, report_fname))
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


def _append(html_fpath, json, keyword):
    # reading html template file
    with open(html_fpath) as f_html:
        html_text = f_html.read()

    # substituting template text with json
    html_text = re.sub('{{ ' + keyword + ' }}', json, html_text)

    # writing substituted html to final file
    with open(html_fpath, 'w') as f_html:
        f_html.write(html_text)

    return html_fpath
