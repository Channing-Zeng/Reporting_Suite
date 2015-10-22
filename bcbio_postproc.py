#!/usr/bin/env python

import sys
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are using %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from site import addsitedir
from os.path import dirname, join, splitext, realpath
this_py_fpath = splitext(__file__)[0] + '.py'
this_py_real_fpath = realpath(this_py_fpath)
project_dir = dirname(this_py_real_fpath)
addsitedir(join(project_dir))
addsitedir(join(project_dir, 'ext_modules'))


if __name__ == '__main__':
    from source.bcbio.main import main
    main()

