import datetime
import json
import os
from os.path import join


def save_total_report(cnf, report_base_name, sample_names, report):
    t = datetime.datetime.now()

    return save(join(cnf['work_dir'], report_base_name + '.json'), {
        'date': t.strftime('%d %B %Y, %A, %H:%M:%S'),
        'report': report,
        'sampleNames': sample_names,
    })


def save(fpath, contents):
    if os.path.exists(fpath):
        os.remove(fpath)

    with open(fpath, 'w') as f:
        json.dump(contents, f, separators=(',', ':'))

    return fpath














