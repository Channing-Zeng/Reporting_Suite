import datetime
from json import dump, JSONEncoder
import os
from os.path import join


class Encoder(JSONEncoder):
    def default(self, o):
        return o.__dict__


def save_total_report(work_dir, report_base_name, report):
    t = datetime.datetime.now()

    return save(join(work_dir, report_base_name + '.json'), dict(
        date=t.strftime('%d %B %Y, %A, %H:%M:%S'),
        report=report,
    ))


def save(fpath, contents):
    if os.path.exists(fpath):
        os.remove(fpath)

    with open(fpath, 'w') as f:
        dump(contents, f, separators=(',', ':'), cls=Encoder)

    return fpath














