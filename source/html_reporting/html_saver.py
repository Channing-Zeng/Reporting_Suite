import os
import shutil
import re
import sys
import datetime
from os.path import join, abspath, dirname, isdir, splitext
from json import dumps, JSONEncoder
import traceback

from source.bcbio_structure import VariantCaller, BCBioSample
from source.file_utils import verify_file, file_transaction
from source.html_reporting import json_saver
from source.file_utils import file_exists
from ext_modules.jsontemplate import jsontemplate
from source.logger import err


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
    'scripts/build_total_report.js',
    'scripts/build_report.js',
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


def write_html_report(json, output_dirpath, report_base_name, caption=''):
    html_fpath = _init_html(output_dirpath, report_base_name + '.html', caption)
    _append(html_fpath, json, 'totalReport')
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


def _init_html(results_dirpath, report_fname, caption=''):
    # Temporary:
    aux_dirpath = join(results_dirpath, aux_dirname)
    if isdir(aux_dirpath):
        shutil.rmtree(aux_dirpath)
    # Temporary.

    with open(template_fpath) as template_file:
        html = template_file.read()

    html = html.replace('{{ caption }}', caption)

    html = _embed_css_and_scripts(html)

    html_fpath = os.path.join(results_dirpath, report_fname)
    if os.path.exists(html_fpath):
        os.remove(html_fpath)

    with open(html_fpath, 'w') as f:
        f.write(html)

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


def _write_static_html_report(work_dir, data_dict, html_fpath):
    with open(static_template_fpath) as f:
        html = f.read()

    html = jsontemplate.expand(html, data_dict)

    html = _embed_css_and_scripts(html)

    with file_transaction(work_dir, html_fpath) as tx:
        with open(tx, 'w') as f:
            f.write(html)

    return html_fpath
