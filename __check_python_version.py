import sys
if not ((2, 7) <= sys.version_info[:2] < (3, 0)):
    sys.exit('Python 2, versions 2.7 and higher is supported '
             '(you are using %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from site import addsitedir
from os.path import abspath, dirname, join, splitext
import subprocess
this_py_fpath = splitext(__file__)[0] + '.py'
link_contents = subprocess.check_output('readlink ' + this_py_fpath, shell=True).strip()
abs_this_py_fpath = abspath(join(dirname(this_py_fpath), link_contents))
project_dir = dirname(abs_this_py_fpath)
addsitedir(join(project_dir))
addsitedir(join(project_dir, 'ext_modules'))